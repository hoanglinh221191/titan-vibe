#!/usr/bin/env python3
"""Audit a segmented 4W52 TITAN-Path mechanical/dynamic pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


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
        row = [float(token) for token in line.split()]
        if not all(math.isfinite(value) for value in row):
            raise RuntimeError(f"non-finite value in {path}")
        rows.append(row)
    if not rows:
        raise RuntimeError(f"no numeric rows in {path}")
    return rows


def read_key_value(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="strict").splitlines():
        line = raw.strip()
        if not line:
            continue
        if "=" not in line:
            raise RuntimeError(f"malformed key-value line in {path}: {line!r}")
        key, value = line.split("=", 1)
        if key in values:
            raise RuntimeError(f"duplicate key {key!r} in {path}")
        values[key] = value
    return values


def gmx_check(gmx: str, flag: str, path: Path) -> dict[str, object]:
    completed = subprocess.run(
        [gmx, "check", flag, str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    output = completed.stdout
    passed = completed.returncode == 0 and "Fatal error" not in output
    return {
        "status": "PASS" if passed else "FAIL",
        "returncode": completed.returncode,
        "output_tail": output.splitlines()[-20:],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=root
        / "docs"
        / "titan_path"
        / "gates"
        / "G4_4w52_pilot_20260825"
        / "run_003",
    )
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--expected-segments", type=int)
    args = parser.parse_args(argv)

    try:
        work = args.work_dir.resolve()
        state = json.loads((work / "controller_state.json").read_text(encoding="utf-8"))
        record_paths = sorted((work / "segments").glob("segment_*/segment_record.json"))
        segment_indices = [int(path.parent.name.removeprefix("segment_")) for path in record_paths]
        if not record_paths or segment_indices != list(range(len(record_paths))):
            raise RuntimeError("segment records must be present and contiguous from index zero")
        if args.expected_segments is not None and len(record_paths) != args.expected_segments:
            raise RuntimeError(
                f"found {len(record_paths)} segment records; expected {args.expected_segments}"
            )
        records = [
            json.loads(path.read_text(encoding="utf-8")) for path in record_paths
        ]
        pullx = read_xvg(work / "pullx.xvg")
        pullf = read_xvg(work / "pullf.xvg")
        times = [row[0] for row in pullx]
        pullf_times = [row[0] for row in pullf]
        strict_time = all(b > a for a, b in zip(times, times[1:]))
        force_time_match = pullf_times == times
        inferred_dt_ps = times[1] - times[0] if len(times) > 1 else math.nan
        expected_end_ps = float(records[-1]["boundary"]["time_ps"])
        expected_rows = int(round(expected_end_ps / inferred_dt_ps)) + 1
        expected_grid = (
            math.isfinite(inferred_dt_ps)
            and inferred_dt_ps > 0.0
            and len(times) == expected_rows
            and abs(times[0]) <= 1.0e-12
            and abs(times[-1] - expected_end_ps) <= 1.0e-8
        )

        handoff_checks: list[dict[str, object]] = []
        records_ok = True
        for index, record in enumerate(records):
            handoff = record["handoff"]
            status = record["titan_status"]
            segment_dir = record_paths[index].parent
            mdrun_log = (segment_dir / "mdrun.log").read_text(
                encoding="utf-8", errors="replace"
            )
            grompp_log = (segment_dir / "grompp.log").read_text(
                encoding="utf-8", errors="replace"
            )
            item = {
                "segment": index,
                "boundary_time_ps": record["boundary"]["time_ps"],
                "direction_change_deg": handoff["direction_change_deg"],
                "extension_residual_nm": handoff["extension_residual_nm"],
                "energy_residual_kj_mol": handoff["energy_residual_kj_mol"],
                "force_magnitude_residual_kj_mol_nm": handoff[
                    "force_magnitude_residual_kj_mol_nm"
                ],
                "registration_rmsd_A": float(status["REGISTRATION_RMSD_A"]),
                "target_clearance_A": float(status["TARGET_MIN_CLEARANCE_A"]),
                "path_progress_A": float(status["PATH_PROGRESS_A"]),
                "target_progress_A": float(status["TARGET_PROGRESS_A"]),
                "mdrun_finished": "Writing final coordinates." in mdrun_log
                and "Fatal error" not in mdrun_log,
                "grompp_clean": "Fatal error" not in grompp_log
                and "WARNING" not in grompp_log,
                "segment_complete_marker": (segment_dir / "segment.complete").is_file(),
            }
            observed_extension_residual_nm: float | None = None
            if index + 1 < len(records):
                boundary_time = float(record["boundary"]["time_ps"])
                matching = [row for row in pullx if abs(row[0] - boundary_time) <= 1.0e-9]
                if len(matching) != 1:
                    raise RuntimeError(
                        f"expected one post-handoff pullx row at {boundary_time} ps"
                    )
                observed_extension = matching[0][1] - matching[0][2]
                observed_extension_residual_nm = (
                    observed_extension - float(handoff["signed_extension_nm"])
                )
                item["observed_post_handoff_extension_nm"] = observed_extension
                item["observed_extension_residual_nm"] = observed_extension_residual_nm
            else:
                item["observed_post_handoff_extension_nm"] = None
                item["observed_extension_residual_nm"] = None
            passed = (
                abs(float(handoff["extension_residual_nm"])) <= 1.0e-12
                and abs(float(handoff["energy_residual_kj_mol"])) <= 1.0e-12
                and abs(float(handoff["force_magnitude_residual_kj_mol_nm"]))
                <= 1.0e-10
                and float(handoff["direction_change_deg"]) <= 20.0 + 1.0e-10
                and float(status["REGISTRATION_RMSD_A"]) <= 2.0
                and float(status["TARGET_MIN_CLEARANCE_A"]) >= 1.5
                and abs(
                    float(status["SNAPSHOT_TIME_PS"])
                    - float(record["boundary"]["time_ps"])
                )
                <= 0.00200001
                and item["mdrun_finished"]
                and item["grompp_clean"]
                and item["segment_complete_marker"]
                and (
                    observed_extension_residual_nm is None
                    or abs(observed_extension_residual_nm) <= 5.0e-6
                )
            )
            item["status"] = "PASS" if passed else "FAIL"
            handoff_checks.append(item)
            records_ok = records_ok and passed

        correction_path = work / "forward_selection_correction.json"
        correction_check: dict[str, object] = {"status": "NOT_APPLICABLE"}
        if correction_path.is_file():
            correction = json.loads(correction_path.read_text(encoding="utf-8"))
            replacement_tpr = (work / str(correction["replacement_tpr"])).resolve()
            corrected_status_path = (
                work / str(correction["corrected_status_file"])
            ).resolve()
            handoff_record_path = (
                work / str(correction["handoff_record_file"])
            ).resolve()
            corrected_status = read_key_value(corrected_status_path)
            corrected_handoff = json.loads(
                handoff_record_path.read_text(encoding="utf-8")
            )["handoff"]
            boundary_time = float(correction["boundary_time_ps"])
            before_records = [
                record
                for record in records
                if abs(float(record["boundary"]["time_ps"]) - boundary_time) <= 1.0e-9
            ]
            if len(before_records) != 1:
                raise RuntimeError(
                    f"expected one recorded pre-correction boundary at {boundary_time} ps"
                )
            before_record = before_records[0]
            replacement_hash = sha256(replacement_tpr)
            after_records = [
                record
                for record in records
                if record["sha256"]["input_tpr"] == replacement_hash
            ]
            if len(after_records) != 1:
                raise RuntimeError("replacement TPR must bind exactly one executed segment")
            after_record = after_records[0]
            corrected_path_progress = float(corrected_status["PATH_PROGRESS_A"])
            corrected_target_progress = float(corrected_status["TARGET_PROGRESS_A"])
            post_path_progress = float(after_record["titan_status"]["PATH_PROGRESS_A"])
            post_target_progress = float(after_record["titan_status"]["TARGET_PROGRESS_A"])
            matching_boundary = [
                row for row in pullx if abs(row[0] - boundary_time) <= 1.0e-9
            ]
            if len(matching_boundary) != 1:
                raise RuntimeError(
                    f"expected one corrected pullx row at {boundary_time} ps"
                )
            observed_extension = matching_boundary[0][1] - matching_boundary[0][2]
            observed_residual = observed_extension - float(
                corrected_handoff["signed_extension_nm"]
            )
            hashes_ok = (
                replacement_hash == correction["sha256"]["replacement_tpr"]
                and sha256(corrected_status_path)
                == correction["sha256"]["corrected_status"]
                and sha256(handoff_record_path)
                == correction["sha256"]["handoff_record"]
            )
            selection_before_ok = (
                abs(
                    float(before_record["titan_status"]["TARGET_PROGRESS_A"])
                    - float(correction["selection_before"]["target_progress_a"])
                )
                <= 1.0e-8
            )
            static_forward_ok = (
                int(corrected_status["FORWARD_ARCHIVE_SIZE"]) > 0
                and corrected_target_progress > corrected_path_progress + 1.0e-6
            )
            dynamic_forward_ok = post_path_progress > corrected_path_progress + 1.0e-6
            post_selection_ok = (
                int(after_record["titan_status"]["FORWARD_ARCHIVE_SIZE"]) > 0
                and post_target_progress > post_path_progress + 1.0e-6
            )
            handoff_ok = (
                abs(float(corrected_handoff["extension_residual_nm"])) <= 1.0e-12
                and abs(float(corrected_handoff["energy_residual_kj_mol"])) <= 1.0e-12
                and abs(
                    float(corrected_handoff["force_magnitude_residual_kj_mol_nm"])
                )
                <= 1.0e-10
                and abs(observed_residual) <= 5.0e-6
            )
            correction_passed = (
                correction.get("schema") == "titan-path-controller-correction-v1"
                and hashes_ok
                and selection_before_ok
                and static_forward_ok
                and dynamic_forward_ok
                and post_selection_ok
                and handoff_ok
            )
            correction_check = {
                "status": "PASS" if correction_passed else "FAIL",
                "boundary_time_ps": boundary_time,
                "replacement_tpr_sha256": replacement_hash,
                "hashes_ok": hashes_ok,
                "pre_fix_target_progress_A": float(
                    before_record["titan_status"]["TARGET_PROGRESS_A"]
                ),
                "corrected_target_progress_A": corrected_target_progress,
                "corrected_forward_archive_size": int(
                    corrected_status["FORWARD_ARCHIVE_SIZE"]
                ),
                "post_segment_path_progress_A": post_path_progress,
                "post_segment_target_progress_A": post_target_progress,
                "post_segment_forward_archive_size": int(
                    after_record["titan_status"]["FORWARD_ARCHIVE_SIZE"]
                ),
                "observed_extension_residual_nm": observed_residual,
                "selection_before_ok": selection_before_ok,
                "static_forward_selection_ok": static_forward_ok,
                "dynamic_forward_progress_ok": dynamic_forward_ok,
                "post_segment_forward_selection_ok": post_selection_ok,
                "force_continuous_handoff_ok": handoff_ok,
            }

        trajectory_check = gmx_check(args.gmx, "-f", work / "run.xtc")
        energy_check = gmx_check(args.gmx, "-e", work / "run.edr")
        required = ("run.cpt", "run.gro", "run.xtc", "run.edr", "run.log")
        files_ok = all((work / name).is_file() and (work / name).stat().st_size > 0 for name in required)
        state_ok = (
            state.get("status") == "INCOMPLETE_MAX_SEGMENTS"
            and state.get("next_segment") == len(records)
        )
        passed = (
            records_ok
            and strict_time
            and force_time_match
            and expected_grid
            and files_ok
            and state_ok
            and trajectory_check["status"] == "PASS"
            and energy_check["status"] == "PASS"
            and correction_check["status"] in ("PASS", "NOT_APPLICABLE")
        )
        report = {
            "schema": "titan-path-g4-4w52-mechanical-pilot-report-v3",
            "status": "PASS" if passed else "FAIL",
            "scope": (
                f"{len(records)} segments ending at {expected_end_ps:g} ps; mechanical/dynamic "
                "pilot only, not publishable sampling"
            ),
            "expected_controller_stop": "INCOMPLETE_MAX_SEGMENTS",
            "controller_state_ok": state_ok,
            "segment_count": len(records),
            "time_start_ps": times[0],
            "time_end_ps": times[-1],
            "inferred_dt_ps": inferred_dt_ps,
            "pull_rows": len(times),
            "strictly_increasing_time": strict_time,
            "pull_force_time_grid_equal": force_time_match,
            "expected_time_grid": expected_grid,
            "all_required_outputs_present": files_ok,
            "handoff_checks": handoff_checks,
            "forward_selection_correction": correction_check,
            "gmx_check_trajectory": trajectory_check,
            "gmx_check_energy": energy_check,
            "sha256": {
                name: sha256(work / name)
                for name in (*required, "pullx.xvg", "pullf.xvg", "controller_state.json")
            },
        }
        output = args.report.resolve() if args.report else work / "report.json"
        output.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
        return 0 if passed else 1
    except (OSError, RuntimeError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"TITAN_PATH_G4_PILOT_REPORT_ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
