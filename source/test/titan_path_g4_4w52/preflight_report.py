#!/usr/bin/env python3
"""Aggregate the accepted H/J route and TITAN 4W52 preflight evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def status(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        key, value = line.split("=", 1)
        if key in values:
            raise RuntimeError(f"duplicate status key {key} in {path}")
        values[key] = value
    return values


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    gate = root / "docs" / "titan_path" / "gates" / "G4_4w52_hj_preflight_20260825"
    parser.add_argument("--gate-dir", type=Path, default=gate)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    try:
        gate = args.gate_dir.resolve()
        route = root / "test" / "titan_path_g4_4w52" / "hj_literature_route_a.tpath"
        derivation = json.loads((gate / "route_derivation_report.json").read_text(encoding="utf-8"))
        geometry_path = gate / "geometry_v2" / "pull_direction_status.txt"
        energy_path = gate / "energy_v2" / "pull_direction_status.txt"
        geometry = status(geometry_path)
        energy = status(energy_path)
        route_hash = sha256(route)
        common_ok = all(
            item.get("SCHEMA") == "TITAN_PATH_STATUS_V1"
            and item.get("STATUS") == "CONTINUE"
            and item.get("ROUTE_ID") == "4w52-hj-literature-route-a"
            and item.get("ROUTE_FILE_SHA256") == route_hash
            and float(item["REGISTRATION_RMSD_A"]) <= 2.0
            and float(item["TARGET_MIN_CLEARANCE_A"])
            >= float(item["MIN_CLEARANCE_THRESHOLD_A"])
            for item in (geometry, energy)
        )
        energy_ok = (
            energy.get("GEOMETRY_ONLY") == "0"
            and float(energy["TARGET_ENERGY_KEY"])
            <= float(energy["ARCHIVE_MIN_ENERGY_KEY"])
            + float(energy["ENERGY_KEY_EPSILON"])
            + 1.0e-8
        )
        geometry_ok = geometry.get("GEOMETRY_ONLY") == "1"
        passed = derivation.get("status") == "PASS" and common_ok and energy_ok and geometry_ok
        report = {
            "schema": "titan-path-g4-4w52-hj-preflight-report-v1",
            "status": "PASS" if passed else "FAIL",
            "route_id": "4w52-hj-literature-route-a",
            "route_sha256": route_hash,
            "route_derivation_status": derivation.get("status"),
            "route_derivation_max_error_A": derivation.get("max_waypoint_error_A"),
            "geometry": {
                "registration_rmsd_A": float(geometry["REGISTRATION_RMSD_A"]),
                "valid_archive_size": int(geometry["VALID_ARCHIVE_SIZE"]),
                "target_progress_A": float(geometry["TARGET_PROGRESS_A"]),
                "target_clearance_A": float(geometry["TARGET_MIN_CLEARANCE_A"]),
            },
            "energy": {
                "registration_rmsd_A": float(energy["REGISTRATION_RMSD_A"]),
                "valid_archive_size": int(energy["VALID_ARCHIVE_SIZE"]),
                "epsilon_eligible_size": int(energy["EPSILON_ELIGIBLE_SIZE"]),
                "archive_min_energy_key": float(energy["ARCHIVE_MIN_ENERGY_KEY"]),
                "selected_energy_key": float(energy["TARGET_ENERGY_KEY"]),
                "target_progress_A": float(energy["TARGET_PROGRESS_A"]),
                "target_clearance_A": float(energy["TARGET_MIN_CLEARANCE_A"]),
            },
            "sha256": {
                "geometry_status": sha256(geometry_path),
                "energy_status": sha256(energy_path),
                "route_derivation_report": sha256(gate / "route_derivation_report.json"),
            },
        }
        output = args.report.resolve() if args.report else gate / "report.json"
        output.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
        return 0 if passed else 1
    except (OSError, RuntimeError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"TITAN_PATH_G4_PREFLIGHT_REPORT_ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
