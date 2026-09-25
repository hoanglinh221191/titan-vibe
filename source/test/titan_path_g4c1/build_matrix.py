#!/usr/bin/env python3
"""Materialize and preflight the 24-cell G4C1 calibration matrix."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GENERATED = HERE / "generated"
BASE_GENERATED = ROOT / "test" / "titan_path_g4c0" / "generated"
INITIAL_REPORT = (
    ROOT
    / "docs"
    / "titan_path"
    / "gates"
    / "G4C1_4w52_2snapshot_2seed_smoke_20260825"
    / "initial_states"
    / "initial_state_report.json"
)
sys.path.insert(0, str(ROOT / "tools"))

from titan_path_controller import (  # noqa: E402
    _controller_contract,
    _input_hashes,
    load_and_validate_config,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_from_generated(repo_relative: str) -> str:
    return os.path.relpath(
        (ROOT / repo_relative).resolve(), GENERATED.resolve()
    ).replace("\\", "/")


def main() -> int:
    matrix_path = HERE / "matrix_manifest.json"
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    if matrix.get("schema") != "titan-path-g4c1-matrix-manifest-v1":
        raise RuntimeError("matrix manifest schema mismatch")
    initial = json.loads(INITIAL_REPORT.read_text(encoding="utf-8"))
    if initial.get("status") != "PASS":
        raise RuntimeError("initial-state report is not PASS")
    states = initial["states"]
    state_ids = matrix["state_ids"]
    arm_ids = matrix["arm_ids"]
    if set(states) != set(state_ids):
        raise RuntimeError("matrix state IDs do not match the initial-state report")
    if len(state_ids) * len(arm_ids) != int(matrix["expected_cells"]):
        raise RuntimeError("matrix cardinality does not match expected_cells")
    GENERATED.mkdir(parents=True, exist_ok=True)
    cells: dict[str, object] = {}
    for state_id in state_ids:
        state_info = states[state_id]
        for arm_id in arm_ids:
            base_path = BASE_GENERATED / f"{arm_id}.json"
            config = json.loads(base_path.read_text(encoding="utf-8"))
            cell_id = f"{state_id}__{arm_id}"
            config["campaign_id"] = matrix["campaign_id"]
            config["campaign_cell"] = cell_id
            config["snapshot_id"] = state_info["snapshot_id"]
            config["seed_id"] = state_info["seed_id"]
            config["seed"] = state_info["seed"]
            config["work_dir"] = (
                "../../../docs/titan_path/gates/"
                + matrix["campaign_id"]
                + "/cells/"
                + cell_id
            )
            config["initial_tpr"] = relative_from_generated(
                state_info["initial_tpr"]
            )
            config["template_mdp"] = relative_from_generated(
                matrix["segment_template"]
            )
            config["titan_command"] = [
                item.replace(
                    "{config_dir}/../arg_",
                    "{config_dir}/../../titan_path_g4c0/arg_",
                )
                for item in config["titan_command"]
            ]
            config_path = GENERATED / f"{cell_id}.json"
            config_path.write_text(
                json.dumps(config, indent=2, sort_keys=True, allow_nan=False)
                + "\n",
                encoding="utf-8",
                newline="\n",
            )
            loaded = load_and_validate_config(config_path)
            cells[cell_id] = {
                "snapshot_id": state_info["snapshot_id"],
                "seed_id": state_info["seed_id"],
                "seed": state_info["seed"],
                "arm_id": arm_id,
                "config": str(config_path.relative_to(ROOT)).replace("\\", "/"),
                "config_sha256": sha256(config_path),
                "input_sha256": _input_hashes(loaded),
                "controller_contract": _controller_contract(loaded),
            }
    unique_initial_hashes = {
        cell["input_sha256"]["initial_tpr"] for cell in cells.values()
    }
    report = {
        "schema": "titan-path-g4c1-build-report-v1",
        "status": (
            "PREFLIGHT_OK"
            if len(cells) == int(matrix["expected_cells"])
            and len(unique_initial_hashes) == len(state_ids)
            else "FAIL"
        ),
        "campaign_id": matrix["campaign_id"],
        "scope": matrix["scope"],
        "matrix_manifest_sha256": sha256(matrix_path),
        "initial_state_report_sha256": sha256(INITIAL_REPORT),
        "cell_count": len(cells),
        "distinct_initial_tpr_hash_count": len(unique_initial_hashes),
        "cells": cells,
    }
    output = GENERATED / "matrix_build_report.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0 if report["status"] == "PREFLIGHT_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())
