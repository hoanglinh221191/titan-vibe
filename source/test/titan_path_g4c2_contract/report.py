#!/usr/bin/env python3
"""Validate the local TITAN-Path direction-contract/performance smoke."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[2]
GATE = (
    ROOT
    / "docs"
    / "titan_path"
    / "gates"
    / "G4C2_contract_performance_smoke_20260825"
)
ATTEMPT = GATE / "F_corridor_one_segment_attempt02"


class GateError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def vector_distance(left: list[float], right: list[float]) -> float:
    return math.sqrt(sum((float(a) - float(b)) ** 2 for a, b in zip(left, right)))


def extract_search_wall(report_text: str) -> float:
    match = re.search(r"Wall time\s*:\s*([0-9.eE+-]+)\s+s", report_text)
    if match is None:
        raise GateError("DE_report.dat is missing the pathway search wall time")
    return float(match.group(1))


def main() -> int:
    try:
        state_path = ATTEMPT / "controller_state.json"
        record_path = ATTEMPT / "segments" / "segment_0000" / "segment_record.json"
        de_report_path = ATTEMPT / "segments" / "segment_0000" / "DE_report.dat"
        state = json.loads(state_path.read_text(encoding="utf-8"))
        record = json.loads(record_path.read_text(encoding="utf-8"))
        status = record["titan_status"]
        if state.get("status") != "INCOMPLETE_MAX_SEGMENTS":
            raise GateError("one-segment controller did not reach its prespecified cap")
        if status.get("STATUS") != "CONTINUE" or status.get("USED_OLD_VECTOR") != "1":
            raise GateError("TITAN did not continue using the active old direction")
        if status.get("TARGET_CLEARANCE_VECTOR") != "NEXT_VECTOR":
            raise GateError("clearance is not bound to the applied vector")
        direction_delta = vector_distance(
            record["requested_direction_unit"], record["applied_direction_unit"]
        )
        if direction_delta > 1.0e-12 or int(record["remaining_direction_increments"]) != 0:
            raise GateError("controller modified the TITAN-evaluated direction")
        applied_turn = float(status["NEXT_VECTOR_ANGLE_DEG"])
        if applied_turn > 20.0000005:
            raise GateError("applied direction exceeds the frozen 20-degree cap")
        clearance = float(status["TARGET_MIN_CLEARANCE_A"])
        if clearance < 1.5:
            raise GateError("applied target violates the 1.5-A clearance gate")

        timing = record["timing_wall_seconds"]
        overhead = sum(
            float(timing[key])
            for key in ("direction_provider", "grompp", "controller_nonprocess")
        )
        segment_total = float(timing["segment_total"])
        provider = record["direction_provider_runtime"]
        search_wall = extract_search_wall(de_report_path.read_text(encoding="utf-8"))
        for native_log in ("gmx_check_xtc.log", "gmx_check_edr.log"):
            log_path = ATTEMPT / native_log
            if not log_path.is_file() or log_path.stat().st_size == 0:
                raise GateError(f"native GROMACS check is missing: {native_log}")

        report = {
            "schema": "titan-path-g4c2-contract-performance-smoke-v1",
            "status": "PASS",
            "scope": "one 2-ps local 4W52 F-arm integration segment",
            "claim_boundary": (
                "mechanical direction-contract and instrumentation evidence only; "
                "local GROMACS 2026.3 timing is not TrytonPlus performance and is "
                "not a method-effect estimate"
            ),
            "controller_status": state["status"],
            "used_old_vector": int(status["USED_OLD_VECTOR"]),
            "route_target_angle_deg": float(status["ROUTE_TARGET_ANGLE_DEG"]),
            "applied_turn_angle_deg": applied_turn,
            "turn_capped": int(status["TURN_CAPPED"]),
            "requested_applied_direction_delta": direction_delta,
            "remaining_direction_increments": int(record["remaining_direction_increments"]),
            "applied_target_clearance_A": clearance,
            "candidate_proposals": int(
                status.get("CANDIDATE_PROPOSALS", status["CANDIDATE_EVALUATIONS"])
            ),
            "candidate_evaluations": int(status["CANDIDATE_EVALUATIONS"]),
            "scan_shell_exclusions": int(status.get("SCAN_SHELL_EXCLUSIONS", 0)),
            "energy_evaluations": int(status["ENERGY_EVALUATIONS"]),
            "clearance_pair_checks": int(status["CLEARANCE_PAIR_CHECKS"]),
            "clearance_rejections": int(status["CLEARANCE_REJECTIONS"]),
            "timing_wall_seconds": {
                "mdrun": float(timing["mdrun"]),
                "direction_provider": float(timing["direction_provider"]),
                "provider_conversion": float(provider["conversion_wall_seconds"]),
                "provider_titan_process": float(provider["titan_wall_seconds"]),
                "pathway_search_core": search_wall,
                "grompp": float(timing["grompp"]),
                "controller_nonprocess": float(timing["controller_nonprocess"]),
                "controller_overhead": overhead,
                "segment_total": segment_total,
                "controller_overhead_fraction": overhead / segment_total,
            },
            "sha256": {
                "controller_state": sha256(state_path),
                "segment_record": sha256(record_path),
                "xtc_check_log": sha256(ATTEMPT / "gmx_check_xtc.log"),
                "edr_check_log": sha256(ATTEMPT / "gmx_check_edr.log"),
            },
            "retained_negative_attempt": "F_corridor_one_segment",
            "negative_attempt_reason": (
                "controller correctly rejected a second interpolation caused by "
                "four-decimal direction serialization; fixed by higher status precision "
                "and a no-second-slerp explicit-route contract"
            ),
        }
        output = GATE / "report.json"
        output.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
        return 0
    except (OSError, KeyError, ValueError, json.JSONDecodeError, GateError) as error:
        print(f"G4C2_CONTRACT_PERFORMANCE_GATE_ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
