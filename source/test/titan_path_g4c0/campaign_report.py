#!/usr/bin/env python3
"""Audit the paired six-arm 4W52 implementation smoke without scientific inference."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "tools"))

from titan_path_controller import (  # noqa: E402
    _controller_contract,
    _input_hashes,
    load_and_validate_config,
)


ARMS = (
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


def read_xvg(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith(("#", "@")):
            continue
        row = [float(value) for value in line.split()]
        if not all(math.isfinite(value) for value in row):
            raise RuntimeError(f"non-finite XVG row in {path}")
        rows.append(row)
    if not rows:
        raise RuntimeError(f"no numeric XVG rows in {path}")
    return rows


def row_at(rows: list[list[float]], time_ps: float) -> list[float]:
    matching = [row for row in rows if abs(row[0] - time_ps) <= 1.0e-9]
    if len(matching) != 1:
        raise RuntimeError(f"expected one row at {time_ps} ps; observed {len(matching)}")
    return matching[0]


def row_before(rows: list[list[float]], time_ps: float) -> list[float]:
    prior = [row for row in rows if row[0] < time_ps - 1.0e-12]
    if not prior:
        raise RuntimeError(f"no row before {time_ps} ps")
    return prior[-1]


def gmx_check(gmx: str, flag: str, path: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [gmx, "check", flag, str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    passed = completed.returncode == 0 and "Fatal error" not in completed.stdout
    return {
        "status": "PASS" if passed else "FAIL",
        "returncode": completed.returncode,
        "output_tail": completed.stdout.splitlines()[-12:],
    }


def vector_distance(lhs: list[float], rhs: list[float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(lhs, rhs)))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument(
        "--campaign-root",
        type=Path,
        default=ROOT
        / "docs"
        / "titan_path"
        / "gates"
        / "G4C0_4w52_six_arm_smoke_20260825",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)

    try:
        campaign_root = args.campaign_root.resolve()
        build_report_path = HERE / "generated" / "campaign_build_report.json"
        build_report = json.loads(build_report_path.read_text(encoding="utf-8"))
        manifest = json.loads(
            (HERE / "campaign_manifest.json").read_text(encoding="utf-8")
        )
        arm_reports: dict[str, Any] = {}
        first_boundaries: dict[str, dict[str, Any]] = {}
        first_requested: dict[str, list[float]] = {}
        first_source_hash: dict[str, str | None] = {}
        all_pass = True
        expected_rows = int(manifest["max_segments"]) * int(
            manifest["segment_steps"]
        ) + 1
        expected_end_ps = (
            int(manifest["max_segments"])
            * int(manifest["segment_steps"])
            * float(manifest["dt_ps"])
        )

        for arm_id in ARMS:
            config_path = HERE / "generated" / f"{arm_id}.json"
            config = load_and_validate_config(config_path)
            arm_dir = Path(config["work_dir"])
            state = json.loads(
                (arm_dir / "controller_state.json").read_text(encoding="utf-8")
            )
            records = [
                json.loads(
                    (
                        arm_dir
                        / "segments"
                        / f"segment_{index:04d}"
                        / "segment_record.json"
                    ).read_text(encoding="utf-8")
                )
                for index in range(int(manifest["max_segments"]))
            ]
            pullx = read_xvg(arm_dir / "pullx.xvg")
            pullf = read_xvg(arm_dir / "pullf.xvg")
            times = [row[0] for row in pullx]
            force_times = [row[0] for row in pullf]
            grid_ok = (
                len(times) == expected_rows
                and times == force_times
                and abs(times[0]) <= 1.0e-12
                and abs(times[-1] - expected_end_ps) <= 1.0e-9
                and all(right > left for left, right in zip(times, times[1:]))
            )
            state_ok = (
                state.get("schema") == "titan-path-controller-state-v2"
                and state.get("status") == "INCOMPLETE_MAX_SEGMENTS"
                and state.get("next_segment") == int(manifest["max_segments"])
                and state.get("input_sha256") == _input_hashes(config)
                and state.get("controller_contract") == _controller_contract(config)
            )
            record_files_ok = True
            record_summaries: list[dict[str, Any]] = []
            for index, record in enumerate(records):
                segment = arm_dir / "segments" / f"segment_{index:04d}"
                mdrun_log = (segment / "mdrun.log").read_text(
                    encoding="utf-8", errors="replace"
                )
                grompp_log = (segment / "grompp.log").read_text(
                    encoding="utf-8", errors="replace"
                )
                record_ok = (
                    record.get("schema") == "titan-path-controller-segment-v2"
                    and record.get("direction_kind") == config["direction_kind"]
                    and record.get("spring_handoff_mode")
                    == config["spring_handoff_mode"]
                    and (segment / "segment.complete").is_file()
                    and "Writing final coordinates." in mdrun_log
                    and "Fatal error" not in mdrun_log
                    and "Fatal error" not in grompp_log
                    and "WARNING" not in grompp_log
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
                    selection = record["validated_route_selection"]
                    record_ok = record_ok and (
                        selection is not None
                        and selection["mode"] == "Advance"
                        and selection["reported_forward_archive_size"] > 0
                        and selection["target_progress_A"]
                        > selection["path_progress_A"] + 1.0e-6
                    )
                record_files_ok = record_files_ok and record_ok
                record_summaries.append(
                    {
                        "segment": index,
                        "status": "PASS" if record_ok else "FAIL",
                        "boundary_time_ps": record["boundary"]["time_ps"],
                        "direction_change_deg": record["handoff"][
                            "direction_change_deg"
                        ],
                        "path_progress_A": record["titan_status"].get(
                            "PATH_PROGRESS_A"
                        ),
                        "target_progress_A": record["titan_status"].get(
                            "TARGET_PROGRESS_A"
                        ),
                        "forward_archive_size": record["titan_status"].get(
                            "FORWARD_ARCHIVE_SIZE"
                        ),
                    }
                )

            first = records[0]
            first_boundaries[arm_id] = first["boundary"]
            first_requested[arm_id] = first["requested_direction_unit"]
            first_source_hash[arm_id] = first["titan_status"].get(
                "SOURCE_VECTOR_SHA256"
            )
            handoff_time = float(first["boundary"]["time_ps"])
            post_x = row_at(pullx, handoff_time)
            pre_f = row_before(pullf, handoff_time)
            post_f = row_at(pullf, handoff_time)
            observed_extension = post_x[1] - post_x[2]
            signed_extension = float(first["handoff"]["signed_extension_nm"])
            if config["spring_handoff_mode"] == "legacy_reset":
                observed_residual = observed_extension
                handoff_observed_ok = abs(observed_extension) <= 5.0e-6
            else:
                observed_residual = observed_extension - signed_extension
                handoff_observed_ok = abs(observed_residual) <= 5.0e-6
            force_jump = abs(post_f[1] - pre_f[1])
            trajectory_check = gmx_check(args.gmx, "-f", arm_dir / "run.xtc")
            energy_check = gmx_check(args.gmx, "-e", arm_dir / "run.edr")
            required = ("run.cpt", "run.gro", "run.xtc", "run.edr", "run.log")
            required_ok = all(
                (arm_dir / name).is_file() and (arm_dir / name).stat().st_size > 0
                for name in required
            )
            arm_passed = (
                state_ok
                and grid_ok
                and record_files_ok
                and handoff_observed_ok
                and required_ok
                and trajectory_check["status"] == "PASS"
                and energy_check["status"] == "PASS"
            )
            all_pass = all_pass and arm_passed
            arm_reports[arm_id] = {
                "status": "PASS" if arm_passed else "FAIL",
                "description": manifest["arms"][arm_id]["description"],
                "direction_kind": config["direction_kind"],
                "spring_handoff_mode": config["spring_handoff_mode"],
                "state_contract_ok": state_ok,
                "time_grid_ok": grid_ok,
                "pull_rows": len(times),
                "time_end_ps": times[-1],
                "record_checks": record_summaries,
                "observed_first_handoff_extension_nm": observed_extension,
                "observed_first_handoff_residual_nm": observed_residual,
                "observed_first_handoff_force_jump_kj_mol_nm": force_jump,
                "gmx_check_trajectory": trajectory_check,
                "gmx_check_energy": energy_check,
                "sha256": {
                    name: sha256(arm_dir / name)
                    for name in (*required, "pullx.xvg", "pullf.xvg")
                },
            }

        reference_boundary = first_boundaries[ARMS[0]]
        paired_boundary_max_delta = 0.0
        for arm_id in ARMS[1:]:
            for key, reference_value in reference_boundary.items():
                value = first_boundaries[arm_id][key]
                if isinstance(reference_value, (int, float)):
                    paired_boundary_max_delta = max(
                        paired_boundary_max_delta,
                        abs(float(value) - float(reference_value)),
                    )
                elif isinstance(reference_value, list):
                    paired_boundary_max_delta = max(
                        paired_boundary_max_delta,
                        vector_distance(value, reference_value),
                    )
        paired_initial_boundary_ok = paired_boundary_max_delta <= 1.0e-12
        original_direction_delta = vector_distance(
            first_requested["B_original_legacy_reset"],
            first_requested["C_original_continuous"],
        )
        original_provider_paired_ok = (
            original_direction_delta <= 1.0e-12
            and first_source_hash["B_original_legacy_reset"]
            == first_source_hash["C_original_continuous"]
            and first_source_hash["B_original_legacy_reset"] is not None
        )
        legacy_jump = arm_reports["B_original_legacy_reset"][
            "observed_first_handoff_force_jump_kj_mol_nm"
        ]
        continuous_jump = arm_reports["C_original_continuous"][
            "observed_first_handoff_force_jump_kj_mol_nm"
        ]
        jump_reduction = (
            float(legacy_jump) / float(continuous_jump)
            if float(continuous_jump) > 0.0
            else None
        )
        legacy_artifact_reproduced = (
            float(legacy_jump) > 10.0
            and float(continuous_jump) < 0.5 * float(legacy_jump)
        )
        route_final_directions = {
            arm_id: json.loads(
                (campaign_root / arm_id / "controller_state.json").read_text(
                    encoding="utf-8"
                )
            )["direction_unit"]
            for arm_id in (
                "D_centerline_geometry",
                "E_centerline_energy",
                "F_corridor_energy",
            )
        }
        route_pair_distances = {
            "D_vs_E": vector_distance(
                route_final_directions["D_centerline_geometry"],
                route_final_directions["E_centerline_energy"],
            ),
            "E_vs_F": vector_distance(
                route_final_directions["E_centerline_energy"],
                route_final_directions["F_corridor_energy"],
            ),
            "D_vs_F": vector_distance(
                route_final_directions["D_centerline_geometry"],
                route_final_directions["F_corridor_energy"],
            ),
        }
        route_arms_distinct = all(value > 1.0e-6 for value in route_pair_distances.values())
        all_pass = (
            all_pass
            and paired_initial_boundary_ok
            and original_provider_paired_ok
            and legacy_artifact_reproduced
            and route_arms_distinct
        )
        report = {
            "schema": "titan-path-g4c0-six-arm-report-v1",
            "status": "PASS" if all_pass else "FAIL",
            "campaign_id": manifest["campaign_id"],
            "scope": manifest["scope"],
            "claim_boundary": "implementation smoke only; not a replica and not scientific evidence of method superiority",
            "manifest_sha256": sha256(HERE / "campaign_manifest.json"),
            "build_report_sha256": sha256(build_report_path),
            "build_preflight_status": build_report["status"],
            "paired_initial_boundary_ok": paired_initial_boundary_ok,
            "paired_initial_boundary_max_delta": paired_boundary_max_delta,
            "original_B_C_provider_paired_ok": original_provider_paired_ok,
            "original_B_C_direction_delta": original_direction_delta,
            "legacy_reset_artifact_reproduced": legacy_artifact_reproduced,
            "legacy_force_jump_kj_mol_nm": legacy_jump,
            "continuous_force_jump_kj_mol_nm": continuous_jump,
            "observed_force_jump_reduction_factor": jump_reduction,
            "route_arm_final_direction_distances": route_pair_distances,
            "route_arms_distinct": route_arms_distinct,
            "arms": arm_reports,
        }
        output = (
            args.output.resolve()
            if args.output is not None
            else campaign_root / "report.json"
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
        return 0 if all_pass else 1
    except (OSError, RuntimeError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"TITAN_PATH_G4C0_REPORT_ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
