#!/usr/bin/env python3
"""Audit all 24 G4C1 calibration cells and their paired implementation invariants."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "test" / "titan_path_g4c0"))

from campaign_report import (  # noqa: E402
    gmx_check,
    read_xvg,
    row_at,
    row_before,
    vector_distance,
)
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


def boundary_delta(lhs: dict[str, Any], rhs: dict[str, Any]) -> float:
    maximum = 0.0
    for key, left in lhs.items():
        right = rhs[key]
        if isinstance(left, (int, float)):
            maximum = max(maximum, abs(float(left) - float(right)))
        elif isinstance(left, list):
            maximum = max(maximum, vector_distance(left, right))
    return maximum


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", default="gmx")
    args = parser.parse_args(argv)
    try:
        build_path = HERE / "generated" / "matrix_build_report.json"
        run_path = HERE / "generated" / "matrix_run_report.json"
        manifest_path = HERE / "matrix_manifest.json"
        build = json.loads(build_path.read_text(encoding="utf-8"))
        run = json.loads(run_path.read_text(encoding="utf-8"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if build.get("status") != "PREFLIGHT_OK" or run.get("status") != "RUN_COMPLETE":
            raise RuntimeError("matrix build/run reports are not complete")
        cells: dict[str, Any] = {}
        paired: dict[str, dict[str, Any]] = {}
        all_cells_pass = True
        max_continuous_extension_residual = 0.0
        max_legacy_extension = 0.0
        min_route_forward_archive: int | None = None
        min_route_clearance = math.inf
        max_registration_rmsd = 0.0
        gmx_pass_count = 0

        for cell_id, build_cell in build["cells"].items():
            config_path = (ROOT / build_cell["config"]).resolve()
            config = load_and_validate_config(config_path)
            cell_dir = Path(config["work_dir"])
            state = json.loads(
                (cell_dir / "controller_state.json").read_text(encoding="utf-8")
            )
            method_stopped = state.get("status") == "STOPPED"
            expected_record_count = 0 if method_stopped else 2
            records = [
                json.loads(
                    (
                        cell_dir
                        / "segments"
                        / f"segment_{index:04d}"
                        / "segment_record.json"
                    ).read_text(encoding="utf-8")
                )
                for index in range(expected_record_count)
            ]
            pullx = read_xvg(cell_dir / "pullx.xvg")
            pullf = read_xvg(cell_dir / "pullf.xvg")
            times = [row[0] for row in pullx]
            force_times = [row[0] for row in pullf]
            expected_end_ps = 2.0 if method_stopped else 4.0
            expected_rows = 1001 if method_stopped else 2001
            grid_ok = (
                len(times) == expected_rows
                and times == force_times
                and abs(times[0]) <= 1.0e-12
                and abs(times[-1] - expected_end_ps) <= 1.0e-9
                and all(right > left for left, right in zip(times, times[1:]))
            )
            common_state_ok = (
                state.get("schema") == "titan-path-controller-state-v2"
                and state.get("input_sha256") == _input_hashes(config)
                and state.get("input_sha256") == build_cell["input_sha256"]
                and state.get("controller_contract") == _controller_contract(config)
                and state.get("controller_contract") == build_cell["controller_contract"]
            )
            if method_stopped:
                terminal_path = Path(state["last_terminal_record"])
                terminal = json.loads(terminal_path.read_text(encoding="utf-8"))
                state_ok = common_state_ok and (
                    state.get("next_segment") == 0
                    and state.get("stop_reason") == "NO_VALID_FORWARD_DIRECTION"
                    and state.get("checkpoint") == str(cell_dir / "run.cpt")
                    and terminal.get("schema")
                    == "titan-path-controller-terminal-v1"
                    and terminal.get("outcome") == "STOP"
                    and terminal.get("reason") == state.get("stop_reason")
                    and terminal.get("reason_detail")
                    == state.get("stop_reason_detail")
                    and (terminal_path.parent / "terminal.complete").is_file()
                )
            else:
                terminal = None
                state_ok = common_state_ok and (
                    state.get("status") == "INCOMPLETE_MAX_SEGMENTS"
                    and state.get("next_segment") == 2
                )
            records_ok = True
            route_records = 0
            for index, record in enumerate(records):
                segment_dir = cell_dir / "segments" / f"segment_{index:04d}"
                mdrun_log = (segment_dir / "mdrun.log").read_text(
                    encoding="utf-8", errors="replace"
                )
                grompp_log = (segment_dir / "grompp.log").read_text(
                    encoding="utf-8", errors="replace"
                )
                next_mdp = (segment_dir / "next.mdp").read_text(encoding="utf-8")
                record_ok = (
                    record.get("schema") == "titan-path-controller-segment-v2"
                    and record.get("direction_kind") == config["direction_kind"]
                    and record.get("spring_handoff_mode")
                    == config["spring_handoff_mode"]
                    and (segment_dir / "segment.complete").is_file()
                    and "Writing final coordinates." in mdrun_log
                    and "Fatal error" not in mdrun_log
                    and "Fatal error" not in grompp_log
                    and "WARNING" not in grompp_log
                    and "pcoupl = C-rescale" in next_mdp
                    and abs(float(record["component_old_coordinate_residual_nm"]))
                    <= 5.0e-5
                )
                if config["spring_handoff_mode"] == "force_continuous":
                    record_ok = record_ok and (
                        abs(float(record["handoff"]["extension_residual_nm"]))
                        <= 1.0e-12
                        and abs(float(record["handoff"]["energy_residual_kj_mol"]))
                        <= 1.0e-12
                        and abs(
                            float(
                                record["handoff"][
                                    "force_magnitude_residual_kj_mol_nm"
                                ]
                            )
                        )
                        <= 1.0e-10
                    )
                else:
                    record_ok = record_ok and (
                        record["handoff"]["schema_version"]
                        == "titan-path-legacy-spring-reset-v1"
                        and record["force_continuous_counterfactual"] is not None
                        and record["mdp_updates"]["pull-coord1-start"] == "yes"
                        and abs(
                            float(record["mdp_updates"]["pull-coord1-init"])
                            + float(config["rate_nm_per_ps"])
                            * float(record["boundary"]["time_ps"])
                        )
                        <= 1.0e-15
                    )
                if config["direction_kind"] == "explicit_route":
                    route_records += 1
                    selection = record["validated_route_selection"]
                    record_ok = record_ok and (
                        selection is not None
                        and selection["mode"] == "Advance"
                        and selection["reported_forward_archive_size"] > 0
                        and selection["target_progress_A"]
                        > selection["path_progress_A"] + 1.0e-6
                    )
                    forward = int(selection["reported_forward_archive_size"])
                    min_route_forward_archive = (
                        forward
                        if min_route_forward_archive is None
                        else min(min_route_forward_archive, forward)
                    )
                    min_route_clearance = min(
                        min_route_clearance,
                        float(record["titan_status"]["TARGET_MIN_CLEARANCE_A"]),
                    )
                    max_registration_rmsd = max(
                        max_registration_rmsd,
                        float(record["titan_status"]["REGISTRATION_RMSD_A"]),
                    )
                records_ok = records_ok and record_ok

            if method_stopped:
                assert terminal is not None
                first_boundary = terminal["boundary"]
                first_requested_direction = None
                first_source_vector_sha256 = None
                observed_extension = None
                observed_handoff_residual = None
                observed_handoff_ok = True
                force_jump = None
                if config["direction_kind"] == "explicit_route":
                    max_registration_rmsd = max(
                        max_registration_rmsd,
                        float(terminal["titan_status"]["REGISTRATION_RMSD_A"]),
                    )
            else:
                first = records[0]
                first_boundary = first["boundary"]
                first_requested_direction = first["requested_direction_unit"]
                first_source_vector_sha256 = first["titan_status"].get(
                    "SOURCE_VECTOR_SHA256"
                )
                handoff_time = float(first["boundary"]["time_ps"])
                post_x = row_at(pullx, handoff_time)
                pre_f = row_before(pullf, handoff_time)
                post_f = row_at(pullf, handoff_time)
                observed_extension = post_x[1] - post_x[2]
                signed_extension = float(first["handoff"]["signed_extension_nm"])
                if config["spring_handoff_mode"] == "legacy_reset":
                    observed_handoff_residual = observed_extension
                    observed_handoff_ok = abs(observed_extension) <= 5.0e-6
                    max_legacy_extension = max(
                        max_legacy_extension, abs(observed_extension)
                    )
                else:
                    observed_handoff_residual = observed_extension - signed_extension
                    observed_handoff_ok = abs(observed_handoff_residual) <= 5.0e-6
                    max_continuous_extension_residual = max(
                        max_continuous_extension_residual,
                        abs(observed_handoff_residual),
                    )
                force_jump = abs(post_f[1] - pre_f[1])
            trajectory_check = gmx_check(args.gmx, "-f", cell_dir / "run.xtc")
            energy_check = gmx_check(args.gmx, "-e", cell_dir / "run.edr")
            gmx_pass_count += int(trajectory_check["status"] == "PASS")
            gmx_pass_count += int(energy_check["status"] == "PASS")
            required = ("run.cpt", "run.gro", "run.xtc", "run.edr", "run.log")
            files_ok = all(
                (cell_dir / name).is_file() and (cell_dir / name).stat().st_size > 0
                for name in required
            )
            cell_pass = (
                state_ok
                and grid_ok
                and records_ok
                and observed_handoff_ok
                and files_ok
                and trajectory_check["status"] == "PASS"
                and energy_check["status"] == "PASS"
            )
            all_cells_pass = all_cells_pass and cell_pass
            summary = {
                "status": "PASS" if cell_pass else "FAIL",
                "snapshot_id": build_cell["snapshot_id"],
                "seed_id": build_cell["seed_id"],
                "arm_id": build_cell["arm_id"],
                "controller_outcome": state["status"],
                "stop_reason": state.get("stop_reason"),
                "stop_reason_detail": state.get("stop_reason_detail"),
                "state_contract_ok": state_ok,
                "time_grid_ok": grid_ok,
                "record_count": len(records),
                "route_record_count": route_records,
                "observed_first_handoff_extension_nm": observed_extension,
                "observed_first_handoff_residual_nm": observed_handoff_residual,
                "observed_first_handoff_force_jump_kj_mol_nm": force_jump,
                "first_boundary": first_boundary,
                "first_requested_direction_unit": first_requested_direction,
                "first_source_vector_sha256": first_source_vector_sha256,
                "final_direction_unit": (
                    None if method_stopped else state["direction_unit"]
                ),
                "gmx_trajectory": trajectory_check["status"],
                "gmx_energy": energy_check["status"],
                "sha256": {
                    "controller_state": sha256(cell_dir / "controller_state.json"),
                    "pullx": sha256(cell_dir / "pullx.xvg"),
                    "pullf": sha256(cell_dir / "pullf.xvg"),
                },
            }
            cells[cell_id] = summary
            paired.setdefault(
                build_cell["snapshot_id"] + "_" + build_cell["seed_id"], {}
            )[build_cell["arm_id"]] = summary

        pair_reports: dict[str, Any] = {}
        all_pairs_pass = True
        legacy_to_continuous_jump_ratios: list[float] = []
        for pair_id, arms in paired.items():
            reference = arms[manifest["arm_ids"][0]]["first_boundary"]
            max_delta = max(
                boundary_delta(reference, arms[arm]["first_boundary"])
                for arm in manifest["arm_ids"][1:]
            )
            bc_direction_delta = vector_distance(
                arms["B_original_legacy_reset"][
                    "first_requested_direction_unit"
                ],
                arms["C_original_continuous"][
                    "first_requested_direction_unit"
                ],
            )
            bc_source_match = (
                arms["B_original_legacy_reset"]["first_source_vector_sha256"]
                == arms["C_original_continuous"]["first_source_vector_sha256"]
                and arms["B_original_legacy_reset"][
                    "first_source_vector_sha256"
                ]
                is not None
            )
            legacy_jump = float(
                arms["B_original_legacy_reset"][
                    "observed_first_handoff_force_jump_kj_mol_nm"
                ]
            )
            continuous_jump = float(
                arms["C_original_continuous"][
                    "observed_first_handoff_force_jump_kj_mol_nm"
                ]
            )
            ratio = legacy_jump / continuous_jump if continuous_jump > 0.0 else math.inf
            legacy_to_continuous_jump_ratios.append(ratio)
            route_vectors = {
                arm: arms[arm]["final_direction_unit"]
                for arm in (
                    "D_centerline_geometry",
                    "E_centerline_energy",
                    "F_corridor_energy",
                )
            }
            if all(vector is not None for vector in route_vectors.values()):
                route_distances = {
                    "D_vs_E": vector_distance(
                        route_vectors["D_centerline_geometry"],
                        route_vectors["E_centerline_energy"],
                    ),
                    "E_vs_F": vector_distance(
                        route_vectors["E_centerline_energy"],
                        route_vectors["F_corridor_energy"],
                    ),
                    "D_vs_F": vector_distance(
                        route_vectors["D_centerline_geometry"],
                        route_vectors["F_corridor_energy"],
                    ),
                }
            else:
                route_distances = None
            pair_pass = (
                max_delta <= 1.0e-12
                and bc_direction_delta <= 1.0e-12
                and bc_source_match
                and all(arms[arm]["status"] == "PASS" for arm in manifest["arm_ids"])
            )
            all_pairs_pass = all_pairs_pass and pair_pass
            pair_reports[pair_id] = {
                "status": "PASS" if pair_pass else "FAIL",
                "paired_initial_boundary_max_delta": max_delta,
                "B_C_direction_delta": bc_direction_delta,
                "B_C_source_vector_hash_match": bc_source_match,
                "legacy_force_jump_kj_mol_nm": legacy_jump,
                "continuous_force_jump_kj_mol_nm": continuous_jump,
                "legacy_to_continuous_jump_ratio": ratio,
                "route_final_direction_distances": route_distances,
                "method_stops": {
                    arm: arms[arm]["stop_reason"]
                    for arm in manifest["arm_ids"]
                    if arms[arm]["controller_outcome"] == "STOPPED"
                },
                "method_stop_details": {
                    arm: arms[arm]["stop_reason_detail"]
                    for arm in manifest["arm_ids"]
                    if arms[arm]["controller_outcome"] == "STOPPED"
                    and arms[arm]["stop_reason_detail"] is not None
                },
            }

        passed = (
            all_cells_pass
            and all_pairs_pass
            and len(cells) == int(manifest["expected_cells"])
            and gmx_pass_count == 2 * int(manifest["expected_cells"])
        )
        outcome_counts: dict[str, int] = {}
        stop_reason_counts: dict[str, int] = {}
        stop_detail_counts: dict[str, int] = {}
        for cell in cells.values():
            outcome = str(cell["controller_outcome"])
            outcome_counts[outcome] = outcome_counts.get(outcome, 0) + 1
            if cell["stop_reason"] is not None:
                reason = str(cell["stop_reason"])
                stop_reason_counts[reason] = stop_reason_counts.get(reason, 0) + 1
            if cell["stop_reason_detail"] is not None:
                detail = str(cell["stop_reason_detail"])
                stop_detail_counts[detail] = stop_detail_counts.get(detail, 0) + 1
        report = {
            "schema": "titan-path-g4c1-matrix-report-v1",
            "status": "PASS" if passed else "FAIL",
            "campaign_id": manifest["campaign_id"],
            "scope": manifest["scope"],
            "claim_boundary": "calibration implementation evidence only; do not estimate method effects, completion probabilities, work distributions, or publishable uncertainty from these cells",
            "matrix_manifest_sha256": sha256(manifest_path),
            "matrix_build_report_sha256": sha256(build_path),
            "matrix_run_report_sha256": sha256(run_path),
            "cell_count": len(cells),
            "pair_count": len(pair_reports),
            "gmx_check_pass_count": gmx_pass_count,
            "controller_outcome_counts": outcome_counts,
            "method_stop_reason_counts": stop_reason_counts,
            "method_stop_detail_reason_counts": stop_detail_counts,
            "max_force_continuous_extension_residual_nm": max_continuous_extension_residual,
            "max_legacy_post_reset_extension_nm": max_legacy_extension,
            "min_route_forward_archive_size": min_route_forward_archive,
            "min_route_target_clearance_A": min_route_clearance,
            "max_route_registration_rmsd_A": max_registration_rmsd,
            "legacy_to_continuous_jump_ratio_range": [
                min(legacy_to_continuous_jump_ratios),
                max(legacy_to_continuous_jump_ratios),
            ],
            "pairs": pair_reports,
            "cells": cells,
        }
        output = (
            ROOT
            / "docs"
            / "titan_path"
            / "gates"
            / manifest["campaign_id"]
            / "report.json"
        )
        output.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
        return 0 if passed else 1
    except (OSError, RuntimeError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"TITAN_PATH_G4C1_MATRIX_REPORT_ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
