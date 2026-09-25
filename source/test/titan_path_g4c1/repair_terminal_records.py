#!/usr/bin/env python3
"""Repair missing STOP/COMPLETE metadata records without rerunning or changing MD."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BUILD_REPORT = HERE / "generated" / "matrix_build_report.json"
sys.path.insert(0, str(ROOT / "tools"))

from titan_path_controller import (  # noqa: E402
    TERMINAL_SCHEMA,
    _controller_contract,
    _input_hashes,
    atomic_write_json,
    load_and_validate_config,
    parse_key_value_file,
    read_xvg_boundary,
    sha256_file,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    build = json.loads(BUILD_REPORT.read_text(encoding="utf-8"))
    repaired: dict[str, object] = {}
    tool_hash = sha256(Path(__file__).resolve())
    for cell_id, build_cell in build["cells"].items():
        config_path = (ROOT / build_cell["config"]).resolve()
        config = load_and_validate_config(config_path)
        work_dir = Path(config["work_dir"])
        state_path = work_dir / "controller_state.json"
        if not state_path.is_file():
            continue
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("status") not in ("STOPPED", "COMPLETE"):
            continue
        if state.get("last_terminal_record"):
            continue
        if state.get("input_sha256") != _input_hashes(config):
            raise RuntimeError(f"input hashes drifted for {cell_id}")
        if state.get("controller_contract") != _controller_contract(config):
            raise RuntimeError(f"controller contract drifted for {cell_id}")
        segment = int(state["final_segment"])
        segment_dir = work_dir / "segments" / f"segment_{segment:04d}"
        current_tpr = Path(state["current_tpr"])
        checkpoint = work_dir / "run.cpt"
        gro = work_dir / "run.gro"
        pullx = work_dir / "pullx.xvg"
        status_path = segment_dir / "pull_direction_status.txt"
        for required in (current_tpr, checkpoint, gro, pullx, status_path):
            if not required.is_file() or required.stat().st_size == 0:
                raise RuntimeError(f"missing terminal evidence for {cell_id}: {required}")
        boundary = read_xvg_boundary(
            pullx,
            time_column=config["pullx_columns"]["time"],
            coordinate_column=config["pullx_columns"]["coordinate"],
            reference_column=config["pullx_columns"]["reference"],
            coordinate_component_columns=config["pullx_columns"][
                "coordinate_components"
            ],
            reference_com_columns=config["pullx_columns"]["reference_com"],
            ligand_com_columns=config["pullx_columns"]["ligand_com"],
        )
        status = parse_key_value_file(status_path)
        expected_outcome = "COMPLETE" if state["status"] == "COMPLETE" else "STOP"
        if status.get("STATUS") != expected_outcome:
            raise RuntimeError(f"terminal outcome mismatch for {cell_id}")
        if status.get("REASON") != state.get(
            "completion_reason" if expected_outcome == "COMPLETE" else "stop_reason"
        ):
            raise RuntimeError(f"terminal reason mismatch for {cell_id}")
        if abs(float(status["SNAPSHOT_TIME_PS"]) - boundary["time_ps"]) > 0.00200001:
            raise RuntimeError(f"terminal time mismatch for {cell_id}")
        md_hashes_before = {
            "checkpoint": sha256_file(checkpoint),
            "gro": sha256_file(gro),
            "pullx": sha256_file(pullx),
        }
        terminal_record = {
            "schema": TERMINAL_SCHEMA,
            "segment": segment,
            "outcome": expected_outcome,
            "reason": status.get("REASON", ""),
            "direction_kind": config["direction_kind"],
            "spring_handoff_mode": config["spring_handoff_mode"],
            "boundary": boundary,
            "titan_status": status,
            "metadata_repair": {
                "schema": "titan-path-terminal-metadata-repair-v1",
                "tool_sha256": tool_hash,
                "md_rerun": False,
            },
            "sha256": {
                "input_tpr": sha256_file(current_tpr),
                **md_hashes_before,
                "titan_status": sha256_file(status_path),
            },
        }
        terminal_path = segment_dir / "terminal_record.json"
        atomic_write_json(terminal_path, terminal_record)
        (segment_dir / "terminal.complete").write_text(
            sha256_file(terminal_path) + "\n",
            encoding="ascii",
            newline="\n",
        )
        state.update(
            checkpoint=str(checkpoint),
            terminal_gro=str(gro),
            last_terminal_record=str(terminal_path),
            terminal_boundary=boundary,
            metadata_repair={
                "schema": "titan-path-terminal-metadata-repair-v1",
                "tool_sha256": tool_hash,
                "md_rerun": False,
                "md_sha256_before": md_hashes_before,
            },
        )
        atomic_write_json(state_path, state)
        md_hashes_after = {
            "checkpoint": sha256_file(checkpoint),
            "gro": sha256_file(gro),
            "pullx": sha256_file(pullx),
        }
        if md_hashes_after != md_hashes_before:
            raise RuntimeError(f"metadata repair changed MD evidence for {cell_id}")
        repaired[cell_id] = {
            "terminal_record": str(terminal_path.relative_to(ROOT)).replace(
                "\\", "/"
            ),
            "terminal_record_sha256": sha256_file(terminal_path),
            "md_sha256_unchanged": True,
        }
    report = {
        "schema": "titan-path-g4c1-terminal-repair-report-v1",
        "status": "PASS",
        "md_rerun": False,
        "repaired_cell_count": len(repaired),
        "cells": repaired,
    }
    output = (
        ROOT
        / "docs"
        / "titan_path"
        / "gates"
        / "G4C1_4w52_2snapshot_2seed_smoke_20260825"
        / "terminal_repair_report.json"
    )
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
