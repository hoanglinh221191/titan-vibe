#!/usr/bin/env python3
"""No-MD NAV-P1 navigation-cost and guided-seed mechanics gate."""

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
    if status.get("MODE") != "Advance":
        raise RuntimeError("NAV-P1 offset fixture must remain in Advance")
    if status.get("NAVIGATION_ENABLED") != "1":
        raise RuntimeError("navigation was not enabled")
    if status.get("NAVIGATION_COST_KIND") != "ALGORITHMIC_DIMENSIONLESS":
        raise RuntimeError("navigation cost kind is not audit-safe")
    if status.get("NAVIGATION_SELECTION_SCOPE") != "ADVANCE_ONLY":
        raise RuntimeError("navigation selection scope drifted")
    if status.get("NAVIGATION_DIRECTIONALITY") != (
        "GUIDED_SEEDS_PLUS_PROGRESS_SELECTION"
    ):
        raise RuntimeError("navigation directionality contract drifted")
    if status.get("NAVIGATION_OLD_DIRECTION_BLEND_FRACTIONS") != (
        "0.25000000,0.50000000,0.75000000"
    ):
        raise RuntimeError("navigation blend schedule drifted")
    if abs(finite(status, "NAVIGATION_RECOVERY_TANGENT_WEIGHT") - 0.25) > 1e-12:
        raise RuntimeError("Recovery tangent weight drifted")
    proposed = int(status["GUIDED_SEED_PROPOSALS"])
    evaluated = int(status["GUIDED_SEED_EVALUATIONS"])
    budget_granted = int(status["GUIDED_SEED_BUDGET_GRANTED"])
    valid = int(status["GUIDED_SEED_VALID"])
    if proposed < 9 or budget_granted > proposed or evaluated > budget_granted or valid < 1:
        raise RuntimeError(
            "guided seed coverage failed: "
            f"proposed={proposed} granted={budget_granted} "
            f"evaluated={evaluated} valid={valid}"
        )
    guide = [finite(status, f"NAVIGATION_GUIDE_{axis}") for axis in "XYZ"]
    if abs(math.sqrt(sum(value * value for value in guide)) - 1.0) > 2e-8:
        raise RuntimeError("navigation guide is not unit length")
    if guide[1] <= 0.0 or abs(guide[0]) + abs(guide[2]) <= 0.5:
        raise RuntimeError("offset fixture did not produce tangent-plus-return guidance")
    if finite(status, "PATH_DISTANCE_A") <= finite(
        status, "TARGET_APPLIED_CORRIDOR_RADIUS_A"
    ):
        raise RuntimeError("offset fixture did not start outside the flat-bottom tube")
    tube = finite(status, "TARGET_NAVIGATION_TUBE_COST")
    back = finite(status, "TARGET_NAVIGATION_BACKTRACK_COST")
    total = finite(status, "TARGET_NAVIGATION_COST")
    if min(tube, back, total) < 0.0 or abs(total - tube - back) > 2e-8:
        raise RuntimeError("navigation cost decomposition is inconsistent")
    probe = finite(status, "NEXT_REFERENCE_DISPLACEMENT_A")
    path_distance = finite(status, "APPLIED_PATH_DISTANCE_A")
    corridor_radius = finite(status, "TARGET_APPLIED_CORRIDOR_RADIUS_A")
    applied_delta = finite(status, "APPLIED_DELTA_PROGRESS_A")
    tube_excess = max(0.0, path_distance - corridor_radius)
    backtrack = max(0.0, -applied_delta)
    expected_tube = 0.5 * finite(status, "NAVIGATION_TUBE_WEIGHT") * (
        tube_excess / probe
    ) ** 2
    expected_back = 0.5 * finite(
        status, "NAVIGATION_BACKTRACK_WEIGHT"
    ) * (backtrack / probe) ** 2
    if (
        abs(finite(status, "TARGET_NAVIGATION_TUBE_EXCESS_A") - tube_excess)
        > 2e-8
        or abs(finite(status, "TARGET_NAVIGATION_BACKTRACK_A") - backtrack)
        > 2e-8
        or abs(tube - expected_tube) > 2e-8
        or abs(back - expected_back) > 2e-8
    ):
        raise RuntimeError("navigation cost formula does not match applied geometry")
    if int(status["NAVIGATION_ELIGIBLE_SIZE"]) < 1:
        raise RuntimeError("navigation epsilon pool is empty")
    if int(status["NAVIGATION_ELIGIBLE_SIZE"]) > int(
        status["ENERGY_ELIGIBLE_SIZE"]
    ):
        raise RuntimeError("navigation pool exceeds the physical-energy pool")
    if finite(status, "TARGET_ENERGY_KEY") > finite(
        status, "ARCHIVE_MIN_ENERGY_KEY"
    ) + finite(status, "ENERGY_KEY_EPSILON") + 2e-8:
        raise RuntimeError("selected candidate is outside the energy epsilon pool")
    if finite(status, "TARGET_NAVIGATION_COST") > finite(
        status, "ENERGY_ELIGIBLE_MIN_NAVIGATION_COST"
    ) + finite(status, "NAVIGATION_COST_EPSILON") + 2e-8:
        raise RuntimeError("selected candidate is outside the navigation epsilon pool")
    if finite(status, "APPLIED_DELTA_PROGRESS_A") <= 1.0e-6:
        raise RuntimeError("selected NAV-P1 candidate did not advance")
    if finite(status, "TARGET_MIN_CLEARANCE_A") < finite(
        status, "MIN_CLEARANCE_THRESHOLD_A"
    ):
        raise RuntimeError("selected NAV-P1 candidate violates hard clearance")
    selected_probe = finite(status, "NEXT_REFERENCE_DISPLACEMENT_A")
    if not any(abs(selected_probe - value) <= 1e-8 for value in PROBE_LADDER):
        raise RuntimeError("selected probe is outside the frozen ladder")
    forbidden: list[str] = []
    for suffix in ("*.xtc", "*.trr", "*.edr", "*.cpt", "*.tpr"):
        forbidden.extend(path.name for path in output.rglob(suffix))
    if forbidden:
        raise RuntimeError(f"NAV-P1 gate produced MD artifacts: {forbidden}")
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
    return {
        key: value
        for key, value in status.items()
        if key not in {"CLEARANCE_PAIR_CHECKS"}
    }


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

    nt1 = run_case(titan, "arg_nav_p1_nt1.in", output / "nt1")
    nt4 = run_case(titan, "arg_nav_p1_nt4.in", output / "nt4")
    energy_nt1 = run_case(
        titan, "arg_nav_p1_energy_nt1.in", output / "energy_nt1"
    )
    energy_nt4 = run_case(
        titan, "arg_nav_p1_energy_nt4.in", output / "energy_nt4"
    )
    budget_nt1 = run_case(
        titan, "arg_nav_p1_budget9_nt1.in", output / "budget9_nt1"
    )
    budget_nt4 = run_case(
        titan, "arg_nav_p1_budget9_nt4.in", output / "budget9_nt4"
    )
    for label, left, right in (
        ("geometry", nt1, nt4),
        ("energy", energy_nt1, energy_nt4),
        ("budget9", budget_nt1, budget_nt4),
    ):
        if comparable(left) != comparable(right):
            raise RuntimeError(
                f"NAV-P1 {label} status differs between nt=1 and nt=4"
            )
        if left["vector_sha256"] != right["vector_sha256"]:
            raise RuntimeError(
                f"NAV-P1 {label} selected vector differs between nt=1 and nt=4"
            )

    energy_status = energy_nt1["status"]
    assert isinstance(energy_status, dict)
    if energy_status.get("GEOMETRY_ONLY") != "0" or int(
        energy_status["ENERGY_EVALUATIONS"]
    ) < 1:
        raise RuntimeError("NAV-P1 energy arm did not exercise physical scoring")
    budget_status = budget_nt1["status"]
    assert isinstance(budget_status, dict)
    if (
        int(budget_status["GUIDED_SEED_PROPOSALS"]) != 9
        or int(budget_status["GUIDED_SEED_BUDGET_GRANTED"]) != 9
        or int(budget_status["GUIDED_SEED_EVALUATIONS"]) != 9
        or int(budget_status["CANDIDATE_EVALUATIONS"]) != 9
        or budget_status.get("TARGET_GUIDED_SEED") != "1"
    ):
        raise RuntimeError(
            "finite budget was not consumed deterministically by the guided prefix"
        )

    status = nt1["status"]
    assert isinstance(status, dict)
    report = {
        "schema": "titan-path-nav-p1-mechanics-gate-v1",
        "status": "PASS_NO_MD_NAV_P1_MECHANICS",
        "titan_sha256": sha256(titan),
        "route_sha256": sha256(FIXTURE / "route_nav_p1_offset.tpath"),
        "probe_ladder_A": PROBE_LADDER,
        "nt1": nt1,
        "nt4": nt4,
        "energy_nt1": energy_nt1,
        "energy_nt4": energy_nt4,
        "budget9_nt1": budget_nt1,
        "budget9_nt4": budget_nt4,
        "thread_deterministic": True,
        "md_executed": False,
    }
    report_path = output / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    report_hash = sha256(report_path)
    (output / "report.json.sha256").write_text(
        f"{report_hash}  report.json\n", encoding="ascii", newline="\n"
    )
    print(
        "PASS TITAN_PATH_NAV_P1_MECHANICS "
        f"guided={status['GUIDED_SEED_VALID']}/{status['GUIDED_SEED_EVALUATIONS']} "
        f"nav_cost={finite(status, 'TARGET_NAVIGATION_COST'):.8f} "
        f"energy_nav_cost={finite(energy_status, 'TARGET_NAVIGATION_COST'):.8f} "
        f"applied_ds_A={finite(status, 'APPLIED_DELTA_PROGRESS_A'):.8f} "
        "nt=1/4 deterministic no_md=true"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
