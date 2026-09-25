#!/usr/bin/env python3
"""G3 energy/epsilon integration and thread-reproducibility gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "test"))

from titan_path_g3_gate import FIXTURE, run_case, sha256  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--titan", type=Path, default=ROOT / "titan.exe")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    nt1 = run_case(
        args.titan.resolve(), FIXTURE / "arg_energy_nt1.in", output / "nt1"
    )
    nt4 = run_case(
        args.titan.resolve(), FIXTURE / "arg_energy_nt4.in", output / "nt4"
    )
    keys = (
        "STATUS",
        "MODE",
        "TARGET_PROGRESS_A",
        "NEXT_VECTOR_X_A",
        "NEXT_VECTOR_Y_A",
        "NEXT_VECTOR_Z_A",
        "CORRIDOR_OFFSET_U_A",
        "CORRIDOR_OFFSET_V_A",
        "TARGET_MIN_CLEARANCE_A",
        "TARGET_RAW_SCORE",
        "TARGET_ENERGY_KEY",
        "ARCHIVE_MIN_ENERGY_KEY",
        "VALID_ARCHIVE_SIZE",
        "EPSILON_ELIGIBLE_SIZE",
        "ENERGY_KEY_EPSILON",
        "GEOMETRY_ONLY",
    )
    values_equal = all(
        nt1["status"].get(key) == nt4["status"].get(key) for key in keys
    )
    epsilon = float(nt1["status"]["ENERGY_KEY_EPSILON"])
    selected = float(nt1["status"]["TARGET_ENERGY_KEY"])
    minimum = float(nt1["status"]["ARCHIVE_MIN_ENERGY_KEY"])
    epsilon_constraint_pass = selected <= minimum + epsilon + 1.0e-8
    energy_mode = nt1["status"].get("GEOMETRY_ONLY") == "0"
    eligible_nonzero = int(nt1["status"]["EPSILON_ELIGIBLE_SIZE"]) > 0
    vector_equal = nt1["vector_text"] == nt4["vector_text"]
    passed = (
        values_equal
        and epsilon_constraint_pass
        and energy_mode
        and eligible_nonzero
        and vector_equal
    )
    report = {
        "schema": "titan-path-g3-energy-epsilon-report-v1",
        "status": "PASS" if passed else "FAIL",
        "titan_sha256": sha256(args.titan.resolve()),
        "thread_values_equal": values_equal,
        "vector_output_equal": vector_equal,
        "energy_mode": energy_mode,
        "archive_min_energy_key": minimum,
        "selected_energy_key": selected,
        "epsilon": epsilon,
        "epsilon_constraint_pass": epsilon_constraint_pass,
        "epsilon_eligible_size": int(nt1["status"]["EPSILON_ELIGIBLE_SIZE"]),
        "selected_progress_A": float(nt1["status"]["TARGET_PROGRESS_A"]),
        "selected_clearance_A": float(nt1["status"]["TARGET_MIN_CLEARANCE_A"]),
        "nt1_status": nt1["status"],
        "nt4_status": nt4["status"],
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
