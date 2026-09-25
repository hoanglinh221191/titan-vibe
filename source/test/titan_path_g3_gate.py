#!/usr/bin/env python3
"""G3 thread-reproducibility gate on the repository example1 snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "test" / "titan_path_g3_example1"


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


def run_case(titan: Path, arg_file: Path, output: Path) -> dict[str, Any]:
    output.mkdir(parents=True)
    completed = subprocess.run(
        [str(titan), "-a", str(arg_file)],
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
            f"TITAN returned {completed.returncode} for {arg_file}; see {output / 'console.txt'}"
        )
    required = (
        output / "pull_direction_status.txt",
        output / "DE_out_vector.dat",
        output / "DE_report.dat",
    )
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"missing required output {path}")
    status = parse_status(required[0])
    if status.get("STATUS") != "CONTINUE" or status.get("REASON") != "OK":
        raise RuntimeError(f"unexpected TITAN status: {status}")
    return {
        "status": status,
        "vector_text": required[1].read_text(encoding="utf-8"),
        "vector_sha256": sha256(required[1]),
        "report_sha256": sha256(required[2]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--titan", type=Path, default=ROOT / "titan.exe")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    nt1 = run_case(args.titan.resolve(), FIXTURE / "arg_nt1.in", output / "nt1")
    nt4 = run_case(args.titan.resolve(), FIXTURE / "arg_nt4.in", output / "nt4")
    comparable_keys = (
        "SCHEMA",
        "STATUS",
        "REASON",
        "MODE",
        "MAX_TURN_ANGLE_DEG",
        "PATH_PROGRESS_A",
        "PATH_TOTAL_A",
        "PATH_DISTANCE_A",
        "ROUTE_KIND",
        "ROUTE_ID",
        "ROUTE_FILE_SHA256",
        "REFERENCE_SHA256",
        "REGISTRATION_MODE",
        "REGISTRATION_RMSD_A",
        "REQUESTED_TARGET_PROGRESS_A",
        "REQUESTED_DELTA_PROGRESS_A",
        "TARGET_PROGRESS_A",
        "APPLIED_PROGRESS_A",
        "APPLIED_DELTA_PROGRESS_A",
        "APPLIED_PATH_DISTANCE_A",
        "REQUESTED_PULL_DISTANCE_A",
        "FIXED_PROBE_ENABLED",
        "CONFIGURED_PROBE_LENGTH_A",
        "CONFIGURED_PROBE_LENGTHS_A",
        "NEXT_REFERENCE_DISPLACEMENT_A",
        "NEXT_VECTOR_X_A",
        "NEXT_VECTOR_Y_A",
        "NEXT_VECTOR_Z_A",
        "CORRIDOR_OFFSET_U_A",
        "CORRIDOR_OFFSET_V_A",
        "CORRIDOR_RADIUS_A",
        "TARGET_MIN_CLEARANCE_A",
        "PROBE_EVALUATIONS",
        "PROBE_CLEARANCE_REJECTIONS",
        "PROBE_MODE_REJECTIONS",
    )
    status_equal = all(
        nt1["status"].get(key) == nt4["status"].get(key)
        for key in comparable_keys
    )
    vector_equal = nt1["vector_text"] == nt4["vector_text"]
    route_hash = sha256(FIXTURE / "route.tpath")
    route_bound = (
        nt1["status"].get("ROUTE_FILE_SHA256") == route_hash
        and nt4["status"].get("ROUTE_FILE_SHA256") == route_hash
    )
    clearance_pass = all(
        float(case["status"]["TARGET_MIN_CLEARANCE_A"]) >= 1.5
        for case in (nt1, nt4)
    )
    passed = status_equal and vector_equal and route_bound and clearance_pass
    report = {
        "schema": "titan-path-g3-thread-report-v1",
        "status": "PASS" if passed else "FAIL",
        "titan_sha256": sha256(args.titan.resolve()),
        "example_snapshot_sha256": sha256(ROOT / "example" / "1" / "complex_2.pdb"),
        "route_sha256": route_hash,
        "status_keys_equal": status_equal,
        "vector_output_equal": vector_equal,
        "route_hash_bound": route_bound,
        "clearance_pass": clearance_pass,
        "nt1": nt1,
        "nt4": nt4,
    }
    (output / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
