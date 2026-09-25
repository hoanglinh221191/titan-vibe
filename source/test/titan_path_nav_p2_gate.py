#!/usr/bin/env python3
"""No-MD NAV-P2 exact swept-translation hard-clearance gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "test" / "titan_path_g3_example1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_manifest(paths: list[Path]) -> dict[str, str]:
    manifest: dict[str, str] = {}
    for path in sorted({item.resolve(strict=True) for item in paths}):
        try:
            label = path.relative_to(ROOT).as_posix()
        except ValueError:
            label = str(path)
        manifest[label] = sha256(path)
    return manifest


def source_manifest() -> dict[str, str]:
    source_files = [path for path in (ROOT / "src").rglob("*") if path.is_file()]
    source_files.extend(
        [
            ROOT / ".gitignore",
            ROOT / "Makefile",
            ROOT / "VERSION",
            ROOT / "tools" / "generate_build_info.py",
            ROOT / "tools" / "titan_path_controller.py",
            Path(__file__),
        ]
    )
    return file_manifest(source_files)


def fixture_manifest() -> dict[str, str]:
    files = [
        FIXTURE / "arg_nav_p2_nt1.in",
        FIXTURE / "arg_nav_p2_nt4.in",
        FIXTURE / "arg_nav_p2_crossing_nt1.in",
        FIXTURE / "arg_nav_p2_crossing_nt4.in",
        FIXTURE / "route_nav_p2_crossing.tpath",
        ROOT / "example" / "1" / "complex_2.pdb",
        ROOT / "example" / "1" / "MOL.top",
        ROOT / "example" / "1" / "index.ndx",
    ]
    files.extend(
        path
        for path in (
            ROOT / "example" / "manuscript_mm_pbsa_tabi" / "input" / "lib_ff"
        ).rglob("*")
        if path.is_file()
    )
    return file_manifest(files)


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


def run_case(
    titan: Path,
    argument: str,
    output: Path,
    *,
    require_interior_only_rejection: bool,
) -> dict[str, object]:
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
    status = parse_status(status_path)
    if status.get("STATUS") != "CONTINUE" or status.get("REASON") != "OK":
        raise RuntimeError(f"unexpected swept-clearance outcome: {status}")
    if status.get("CLEARANCE_MODE") != "swept-linear-center":
        raise RuntimeError("swept clearance was not active")
    if status.get("CLEARANCE_SCOPE") != "CONTINUOUS_LINEAR_TRANSLATION":
        raise RuntimeError("swept clearance scope drifted")
    if (
        status.get("CLEARANCE_COORDINATE_FRAME")
        != "UPSTREAM_PBC_RECONSTRUCTED_CARTESIAN"
        or status.get("CLEARANCE_PERIODIC_DISTANCE_POLICY")
        != "DIRECT_IN_RECONSTRUCTED_FRAME"
        or status.get("CLEARANCE_PBC_PREPROCESSOR") != "pbcn"
    ):
        raise RuntimeError("swept clearance PBC-frame contract drifted")
    probe_evaluations = int(status["PROBE_EVALUATIONS"])
    swept_probe_evaluations = int(status["SWEPT_CLEARANCE_PROBE_EVALUATIONS"])
    pair_checks = int(status["CLEARANCE_PAIR_CHECKS"])
    swept_pair_checks = int(status["SWEPT_CLEARANCE_PAIR_CHECKS"])
    swept_rejections = int(status["SWEPT_CLEARANCE_REJECTIONS"])
    interior_only_rejections = int(status["SWEPT_INTERIOR_ONLY_REJECTIONS"])
    if (
        probe_evaluations < 1
        or swept_probe_evaluations != probe_evaluations
        or swept_pair_checks != pair_checks
        or pair_checks < 1
        or not 0 < swept_rejections <= swept_probe_evaluations
        or not 0 <= interior_only_rejections <= swept_rejections
    ):
        raise RuntimeError("swept-clearance counters are inconsistent")
    if require_interior_only_rejection and interior_only_rejections == 0:
        raise RuntimeError(
            "integration fixture did not reject a sweep whose two endpoints were safe"
        )
    swept = finite(status, "TARGET_SWEPT_MIN_CLEARANCE_A")
    endpoint = finite(status, "TARGET_ENDPOINT_MIN_CLEARANCE_A")
    active = finite(status, "TARGET_MIN_CLEARANCE_A")
    threshold = finite(status, "MIN_CLEARANCE_THRESHOLD_A")
    fraction = finite(status, "TARGET_SWEPT_MIN_FRACTION")
    if status.get("TARGET_SWEPT_MIN_COMPLETE") != "1":
        raise RuntimeError("selected swept minimum is not globally complete")
    if (
        abs(active - swept) > 2e-8
        or swept > endpoint + 2e-8
        or swept < threshold - 2e-8
        or not 0.0 <= fraction <= 1.0
    ):
        raise RuntimeError("selected swept-clearance witness is inconsistent")
    for key in (
        "TARGET_SWEPT_MIN_RECEPTOR_INDEX",
        "TARGET_SWEPT_MIN_LIGAND_INDEX",
    ):
        if int(status[key]) < 0:
            raise RuntimeError(f"negative witness index {key}")
    if finite(status, "APPLIED_DELTA_PROGRESS_A") <= 1e-6:
        raise RuntimeError("NAV-P2 selected candidate did not advance")
    if abs(finite(status, "NEXT_REFERENCE_DISPLACEMENT_A") - 0.25) > 2e-8:
        raise RuntimeError(
            "crossing fixture did not exercise swept-clearance trust-region fallback"
        )
    forbidden: list[str] = []
    for suffix in ("*.xtc", "*.trr", "*.edr", "*.cpt", "*.tpr"):
        forbidden.extend(path.name for path in output.rglob(suffix))
    if forbidden:
        raise RuntimeError(f"NAV-P2 gate produced MD artifacts: {forbidden}")
    return {
        "argument": argument,
        "status": status,
        "status_sha256": sha256(status_path),
        "vector_sha256": sha256(output / "DE_out_vector.dat"),
        "report_sha256": sha256(output / "DE_report.dat"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--titan", type=Path, default=ROOT / "titan.exe")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build-toolchain-label", default="not-recorded",
                        help="Explicit label for the build invocation; not inferred from the platform")
    args = parser.parse_args()
    titan = args.titan.resolve(strict=True)
    output = args.output.resolve()
    if output.exists():
        raise RuntimeError(f"refusing to reuse output path: {output}")
    output.mkdir(parents=True)

    nt1 = run_case(
        titan,
        "arg_nav_p2_nt1.in",
        output / "standard_probe" / "nt1",
        require_interior_only_rejection=False,
    )
    nt4 = run_case(
        titan,
        "arg_nav_p2_nt4.in",
        output / "standard_probe" / "nt4",
        require_interior_only_rejection=False,
    )
    status_1 = nt1["status"]
    status_4 = nt4["status"]
    assert isinstance(status_1, dict) and isinstance(status_4, dict)
    if status_1 != status_4:
        raise RuntimeError("NAV-P2 status differs between nt=1 and nt=4")
    if nt1["vector_sha256"] != nt4["vector_sha256"]:
        raise RuntimeError("NAV-P2 selected vector differs between nt=1 and nt=4")
    crossing_nt1 = run_case(
        titan,
        "arg_nav_p2_crossing_nt1.in",
        output / "endpoint_safe_crossing" / "nt1",
        require_interior_only_rejection=True,
    )
    crossing_nt4 = run_case(
        titan,
        "arg_nav_p2_crossing_nt4.in",
        output / "endpoint_safe_crossing" / "nt4",
        require_interior_only_rejection=True,
    )
    crossing_status_1 = crossing_nt1["status"]
    crossing_status_4 = crossing_nt4["status"]
    assert isinstance(crossing_status_1, dict) and isinstance(
        crossing_status_4, dict
    )
    if crossing_status_1 != crossing_status_4:
        raise RuntimeError("crossing status differs between nt=1 and nt=4")
    if crossing_nt1["vector_sha256"] != crossing_nt4["vector_sha256"]:
        raise RuntimeError("crossing selected vector differs between nt=1 and nt=4")

    version = subprocess.run(
        [str(titan), "--version"],
        cwd=output,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if version.returncode != 0:
        raise RuntimeError("TITAN --version failed during provenance capture")
    build_manifest_path = (
        ROOT / "build" / "native" / "generated" / "titan_build_manifest.json"
    ).resolve(strict=True)
    build_manifest = json.loads(build_manifest_path.read_text(encoding="utf-8"))
    report = {
        "schema": "titan-path-nav-p2-swept-clearance-gate-v2",
        "status": "PASS_NO_MD_NAV_P2_SWEPT_CLEARANCE",
        "titan_sha256": sha256(titan),
        "build_identity": {
            "build_request_profile": "native",
            "build_toolchain_label": args.build_toolchain_label,
            "generated_manifest": build_manifest,
            "generated_manifest_sha256": sha256(build_manifest_path),
            "titan_version_output": version.stdout.strip(),
            "python": sys.version,
            "platform": platform.platform(),
        },
        "source_manifest": source_manifest(),
        "fixture_manifest": fixture_manifest(),
        "standard_probe": {"nt1": nt1, "nt4": nt4},
        "endpoint_safe_crossing": {
            "nt1": crossing_nt1,
            "nt4": crossing_nt4,
        },
        "thread_deterministic": True,
        "clearance_definition": "exact continuous atom-center distance over rigid linear translation in one upstream PBC-reconstructed Cartesian frame",
        "periodic_distance_policy": "direct distance in the upstream reconstructed frame; no implicit MIC inside the sweep primitive",
        "soft_steric_enabled": False,
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
        "PASS TITAN_PATH_NAV_P2_SWEPT_CLEARANCE "
        f"swept_A={finite(status_1, 'TARGET_SWEPT_MIN_CLEARANCE_A'):.8f} "
        f"endpoint_A={finite(status_1, 'TARGET_ENDPOINT_MIN_CLEARANCE_A'):.8f} "
        f"tau={finite(status_1, 'TARGET_SWEPT_MIN_FRACTION'):.8f} "
        f"rejections={status_1['SWEPT_CLEARANCE_REJECTIONS']}/"
        f"{status_1['SWEPT_CLEARANCE_PROBE_EVALUATIONS']} "
        "crossing_interior_only="
        f"{crossing_status_1['SWEPT_INTERIOR_ONLY_REJECTIONS']} "
        "nt=1/4 deterministic no_md=true"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
