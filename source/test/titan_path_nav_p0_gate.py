#!/usr/bin/env python3
"""No-MD NAV-P0 fixed-probe, Recovery, and thread-determinism gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "test" / "titan_path_g3_example1"
PROBE_LADDER = [1.0, 0.5, 0.25]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_status(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        key, value = line.split("=", 1)
        if key in values:
            raise RuntimeError(f"duplicate status key {key}")
        values[key] = value
    return values


def finite(status: dict[str, str], key: str) -> float:
    try:
        value = float(status[key])
    except (KeyError, ValueError) as error:
        raise RuntimeError(f"invalid or missing {key}") from error
    if not math.isfinite(value):
        raise RuntimeError(f"non-finite {key}")
    return value


def run_case(titan: Path, argument: str, output: Path) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=False)
    completed = subprocess.run(
        [str(titan), "-a", str(FIXTURE / argument)],
        cwd=output,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    (output / "console.txt").write_text(
        completed.stdout, encoding="utf-8", newline="\n"
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"TITAN returned {completed.returncode} for {argument}; see {output}"
        )
    status_path = output / "pull_direction_status.txt"
    if not status_path.is_file():
        raise RuntimeError(f"missing {status_path}")
    status = parse_status(status_path)
    if status.get("STATUS") != "CONTINUE" or status.get("REASON") != "OK":
        raise RuntimeError(f"unexpected status for {argument}: {status}")
    if status.get("FIXED_PROBE_ENABLED") != "1":
        raise RuntimeError(f"fixed probe was not enabled for {argument}")
    observed_ladder = [
        float(value) for value in status["CONFIGURED_PROBE_LENGTHS_A"].split(",")
    ]
    if observed_ladder != PROBE_LADDER:
        raise RuntimeError(f"probe ladder drift for {argument}: {observed_ladder}")
    displacement = finite(status, "NEXT_REFERENCE_DISPLACEMENT_A")
    vector = [finite(status, f"NEXT_VECTOR_{axis}_A") for axis in "XYZ"]
    vector_length = math.sqrt(sum(value * value for value in vector))
    if abs(vector_length - displacement) > 2.0e-8:
        raise RuntimeError(f"probe/vector mismatch for {argument}")
    if not any(abs(displacement - value) <= 1.0e-8 for value in PROBE_LADDER):
        raise RuntimeError(f"selected probe is outside the frozen ladder: {displacement}")
    path_progress = finite(status, "PATH_PROGRESS_A")
    applied_progress = finite(status, "APPLIED_PROGRESS_A")
    target_progress = finite(status, "TARGET_PROGRESS_A")
    requested_target_progress = finite(status, "REQUESTED_TARGET_PROGRESS_A")
    applied_delta = finite(status, "APPLIED_DELTA_PROGRESS_A")
    if abs(requested_target_progress - target_progress) > 1.0e-8:
        raise RuntimeError("TARGET_PROGRESS_A lost its V1 requested-target meaning")
    if abs(applied_progress - path_progress - applied_delta) > 2.0e-8:
        raise RuntimeError("applied progress accounting mismatch")
    forbidden = []
    for suffix in ("*.xtc", "*.trr", "*.edr", "*.cpt", "*.tpr"):
        forbidden.extend(path.name for path in output.rglob(suffix))
    if forbidden:
        raise RuntimeError(f"NAV-P0 gate unexpectedly produced MD artifacts: {forbidden}")
    return {
        "argument": argument,
        "status": status,
        "status_sha256": sha256(status_path),
        "vector_sha256": sha256(output / "DE_out_vector.dat"),
        "report_sha256": sha256(output / "DE_report.dat"),
    }


def comparable(case: dict[str, object]) -> dict[str, str]:
    status = case["status"]
    assert isinstance(status, dict)
    excluded = {
        "CLEARANCE_PAIR_CHECKS",
    }
    return {key: value for key, value in status.items() if key not in excluded}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--titan", type=Path, default=ROOT / "titan.exe")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    titan = args.titan.resolve(strict=True)
    output = args.output.resolve()
    if output.exists():
        raise RuntimeError(f"refusing to reuse output path: {output}")
    output.mkdir(parents=True)

    advance_1 = run_case(titan, "arg_nav_p0_nt1.in", output / "advance_nt1")
    advance_4 = run_case(titan, "arg_nav_p0_nt4.in", output / "advance_nt4")
    recovery_1 = run_case(
        titan, "arg_nav_p0_recovery_nt1.in", output / "recovery_nt1"
    )
    recovery_4 = run_case(
        titan, "arg_nav_p0_recovery_nt4.in", output / "recovery_nt4"
    )
    recovery_short_1 = run_case(
        titan,
        "arg_nav_p0_recovery_short_nt1.in",
        output / "recovery_short_nt1",
    )
    recovery_short_4 = run_case(
        titan,
        "arg_nav_p0_recovery_short_nt4.in",
        output / "recovery_short_nt4",
    )

    for name, left, right in (
        ("advance", advance_1, advance_4),
        ("recovery", recovery_1, recovery_4),
        ("recovery-short", recovery_short_1, recovery_short_4),
    ):
        if comparable(left) != comparable(right):
            raise RuntimeError(f"{name} status differs between nt=1 and nt=4")
        if left["vector_sha256"] != right["vector_sha256"]:
            raise RuntimeError(f"{name} vector differs between nt=1 and nt=4")

    advance_status = advance_1["status"]
    recovery_status = recovery_1["status"]
    recovery_short_status = recovery_short_1["status"]
    assert (
        isinstance(advance_status, dict)
        and isinstance(recovery_status, dict)
        and isinstance(recovery_short_status, dict)
    )
    if advance_status.get("MODE") != "Advance" or finite(
        advance_status, "APPLIED_DELTA_PROGRESS_A"
    ) <= 1.0e-6:
        raise RuntimeError("Advance did not make strict applied progress")
    if advance_status.get("RECOVERY_DIRECT_PROBE_USED") != "0":
        raise RuntimeError("Advance unexpectedly used the Recovery probe")
    if recovery_status.get("MODE") != "Recovery":
        raise RuntimeError("synthetic Recovery case did not enter Recovery")
    if recovery_status.get("RECOVERY_DIRECT_PROBE_USED") != "1":
        raise RuntimeError("empty-forward Recovery did not use its direct probe")
    if finite(recovery_status, "APPLIED_PATH_DISTANCE_A") >= finite(
        recovery_status, "PATH_DISTANCE_A"
    ) - 1.0e-6:
        raise RuntimeError("Recovery direct probe did not approach the pathway")
    if int(recovery_status["CANDIDATE_PROPOSALS"]) != 1:
        raise RuntimeError("Recovery direct probe must be a single deterministic proposal")
    if abs(
        finite(recovery_short_status, "NEXT_REFERENCE_DISPLACEMENT_A") - 0.5
    ) > 1.0e-8:
        raise RuntimeError("Recovery trust region did not fall back from 1.0 to 0.5 A")
    if int(recovery_short_status["PROBE_MODE_REJECTIONS"]) < 1:
        raise RuntimeError("short-probe fixture did not exercise mode rejection")
    if finite(recovery_short_status, "APPLIED_PATH_DISTANCE_A") >= finite(
        recovery_short_status, "PATH_DISTANCE_A"
    ) - 1.0e-6:
        raise RuntimeError("short Recovery probe did not approach the pathway")

    report = {
        "schema": "titan-path-nav-p0-gate-v1",
        "status": "PASS_NO_MD_NAV_P0",
        "titan_sha256": sha256(titan),
        "probe_ladder_A": PROBE_LADDER,
        "advance": advance_1,
        "advance_nt4": advance_4,
        "recovery": recovery_1,
        "recovery_nt4": recovery_4,
        "recovery_short": recovery_short_1,
        "recovery_short_nt4": recovery_short_4,
        "thread_deterministic": True,
        "md_executed": False,
    }
    report_path = output / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    (output / "report.json.sha256").write_text(
        f"{sha256(report_path)}  report.json\n", encoding="ascii", newline="\n"
    )
    print(
        "PASS TITAN_PATH_NAV_P0 "
        f"advance_probe_A={finite(advance_status, 'NEXT_REFERENCE_DISPLACEMENT_A'):.2f} "
        f"recovery_distance_A={finite(recovery_status, 'PATH_DISTANCE_A'):.6f}->"
        f"{finite(recovery_status, 'APPLIED_PATH_DISTANCE_A'):.6f} "
        f"recovery_fallback_A={finite(recovery_short_status, 'NEXT_REFERENCE_DISPLACEMENT_A'):.2f} "
        "nt=1/4 deterministic no_md=true"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
