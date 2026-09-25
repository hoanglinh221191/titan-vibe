#!/usr/bin/env python3
"""Materialize and preflight the frozen six-arm 4W52 implementation smoke."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GENERATED = HERE / "generated"
sys.path.insert(0, str(ROOT / "tools"))

from titan_path_controller import (  # noqa: E402
    CONFIG_SCHEMA,
    _controller_contract,
    _input_hashes,
    load_and_validate_config,
)


EXPECTED_ARMS = (
    "A_fixed",
    "B_original_legacy_reset",
    "C_original_continuous",
    "D_centerline_geometry",
    "E_centerline_energy",
    "F_corridor_energy",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_from_generated(repo_relative: str) -> str:
    absolute = (ROOT / repo_relative).resolve()
    return os.path.relpath(absolute, GENERATED.resolve()).replace("\\", "/")


def route_command(arg_file: str) -> list[str]:
    return [
        "python3",
        "{config_dir}/../../../tools/titan_path_run_snapshot.py",
        "--gmx",
        "gmx",
        "--titan",
        "{config_dir}/../../../titan_wsl",
        "--gro",
        "{gro}",
        "--snapshot",
        "{segment_dir}/snapshot.pdb",
        "--time-ps",
        "{time_ps}",
        "--old-direction",
        "{old_direction_x}",
        "{old_direction_y}",
        "{old_direction_z}",
        "--arg-file",
        "{config_dir}/" + relative_from_generated(arg_file),
        "--topology",
        "{config_dir}/../../../.codex_tmp/DTW_Path_Classification/T4L_benzene/opes_explore/template/lysozyme.top",
        "--index",
        "{config_dir}/../../../.codex_tmp/DTW_Path_Classification/T4L_benzene/opes_explore/template/index.ndx",
        "--route",
        "{config_dir}/../../titan_path_g4_4w52/hj_literature_route_a.tpath",
        "--route-id",
        "4w52-hj-literature-route-a",
        "--fflib",
        "{config_dir}/../../../example/manuscript_mm_pbsa_tabi/input/lib_ff",
    ]


def fixed_command(provider_id: str, direction: list[float], max_turn: float) -> list[str]:
    return [
        "python3",
        "{config_dir}/../../../tools/titan_direction_adapter.py",
        "--mode",
        "fixed",
        "--provider-id",
        provider_id,
        "--time-ps",
        "{time_ps}",
        "--max-turn-angle-deg",
        f"{max_turn:.17g}",
        "--direction",
        *(f"{value:.17g}" for value in direction),
    ]


def unconstrained_command(
    provider_id: str, arg_file: str, max_turn: float
) -> list[str]:
    return [
        "python3",
        "{config_dir}/../../../tools/titan_direction_adapter.py",
        "--mode",
        "unconstrained",
        "--provider-id",
        provider_id,
        "--time-ps",
        "{time_ps}",
        "--max-turn-angle-deg",
        f"{max_turn:.17g}",
        "--old-direction",
        "{old_direction_x}",
        "{old_direction_y}",
        "{old_direction_z}",
        "--gmx",
        "gmx",
        "--titan",
        "{config_dir}/../../../titan_wsl",
        "--gro",
        "{gro}",
        "--snapshot",
        "{segment_dir}/snapshot.pdb",
        "--arg-file",
        "{config_dir}/" + relative_from_generated(arg_file),
        "--topology",
        "{config_dir}/../../../.codex_tmp/DTW_Path_Classification/T4L_benzene/opes_explore/template/lysozyme.top",
        "--index",
        "{config_dir}/../../../.codex_tmp/DTW_Path_Classification/T4L_benzene/opes_explore/template/index.ndx",
        "--fflib",
        "{config_dir}/../../../example/manuscript_mm_pbsa_tabi/input/lib_ff",
    ]


def build_config(manifest: dict[str, object], arm_id: str, arm: dict[str, object]) -> dict[str, object]:
    direction_kind = str(arm["direction_kind"])
    max_turn = float(arm["max_turn_angle_deg"])
    initial_direction = [float(value) for value in manifest["initial_direction"]]  # type: ignore[arg-type]
    if direction_kind == "explicit_route":
        command = route_command(str(arm["arg_file"]))
    elif direction_kind == "unconstrained":
        command = unconstrained_command(
            str(arm["provider_id"]), str(arm["arg_file"]), max_turn
        )
    else:
        command = fixed_command(str(arm["provider_id"]), initial_direction, max_turn)
    config: dict[str, object] = {
        "schema": CONFIG_SCHEMA,
        "campaign_id": manifest["campaign_id"],
        "campaign_arm": arm_id,
        "work_dir": "../../../docs/titan_path/gates/"
        + str(manifest["campaign_id"])
        + "/"
        + str(arm.get("work_dir_suffix", arm_id)),
        "gmx_command": ["gmx"],
        "titan_command": command,
        "template_mdp": relative_from_generated(str(manifest["template_mdp"])),
        "topology": relative_from_generated(str(manifest["topology"])),
        "route_file": relative_from_generated(str(manifest["route_file"])),
        "index": relative_from_generated(str(manifest["index"])),
        "initial_tpr": relative_from_generated(str(manifest["initial_tpr"])),
        "initial_checkpoint": None,
        "direction_kind": direction_kind,
        "spring_handoff_mode": arm["spring_handoff_mode"],
        "initial_direction": initial_direction,
        "initial_step": 0,
        "initial_simulation_part": 1,
        "dt_ps": manifest["dt_ps"],
        "segment_steps": manifest["segment_steps"],
        "max_segments": manifest["max_segments"],
        "rate_nm_per_ps": manifest["rate_nm_per_ps"],
        "spring_k_kj_mol_nm2": manifest["spring_k_kj_mol_nm2"],
        "max_turn_angle_deg": max_turn,
        "pullx_columns": {
            "time": 0,
            "coordinate": 1,
            "reference": 2,
            "coordinate_components": [3, 4, 5],
            "reference_com": [6, 7, 8],
            "ligand_com": [9, 10, 11],
        },
        "mdrun_extra": [
            "-ntmpi",
            "1",
            "-ntomp",
            "1",
            "-reprod",
            "-pin",
            "off",
        ],
        "grompp_extra": [],
        "environment": {"GMX_MAXBACKUP": "-1"},
    }
    if direction_kind == "explicit_route":
        config["expected_route_id"] = manifest["route_id"]
        config["expected_reference_sha256"] = manifest["reference_sha256"]
    else:
        config["expected_provider_id"] = arm["provider_id"]
    return config


def main() -> int:
    manifest_path = HERE / "campaign_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "titan-path-g4c0-six-arm-manifest-v1":
        raise RuntimeError("campaign manifest schema mismatch")
    arms = manifest.get("arms")
    if not isinstance(arms, dict) or tuple(arms) != EXPECTED_ARMS:
        raise RuntimeError("campaign arms or ordering do not match the frozen A-F matrix")
    GENERATED.mkdir(parents=True, exist_ok=True)
    report_arms: dict[str, object] = {}
    for arm_id in EXPECTED_ARMS:
        arm = arms[arm_id]
        if not isinstance(arm, dict):
            raise RuntimeError(f"arm {arm_id} must be an object")
        config = build_config(manifest, arm_id, arm)
        config_path = GENERATED / f"{arm_id}.json"
        config_path.write_text(
            json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        loaded = load_and_validate_config(config_path)
        report_arms[arm_id] = {
            "description": arm["description"],
            "config": str(config_path.relative_to(ROOT)).replace("\\", "/"),
            "config_sha256": sha256(config_path),
            "input_sha256": _input_hashes(loaded),
            "controller_contract": _controller_contract(loaded),
        }
    report = {
        "schema": "titan-path-g4c0-build-report-v1",
        "status": "PREFLIGHT_OK",
        "campaign_id": manifest["campaign_id"],
        "scope": manifest["scope"],
        "manifest_sha256": sha256(manifest_path),
        "arm_count": len(report_arms),
        "arms": report_arms,
    }
    output = GENERATED / "campaign_build_report.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
