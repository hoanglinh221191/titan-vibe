#!/usr/bin/env python3
"""Run the preflighted G4C1 cells with bounded local concurrency."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BUILD_REPORT = HERE / "generated" / "matrix_build_report.json"
ALLOWED_METHOD_STOPS = {
    "NO_VALID_FORWARD_DIRECTION",
    "NO_ENERGETICALLY_SELECTABLE_FORWARD_TARGET",
    "NO_REACHABLE_FORWARD_TARGET",
    "TERMINAL_REACHED_WITHOUT_PHYSICAL_EXIT",
}


def classify_state(config_path: Path, returncode: int | None) -> dict[str, object]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    work_dir = (config_path.parent / config["work_dir"]).resolve()
    state_path = work_dir / "controller_state.json"
    if not state_path.is_file():
        return {
            "accepted_terminal_state": False,
            "classification": "MISSING_STATE",
            "controller_status": None,
            "next_segment": None,
            "stop_reason": None,
        }
    state = json.loads(state_path.read_text(encoding="utf-8"))
    status = state.get("status")
    next_segment = state.get("next_segment")
    stop_reason = state.get("stop_reason")
    expected_max = int(config["max_segments"])
    if status == "INCOMPLETE_MAX_SEGMENTS" and next_segment == expected_max:
        accepted = returncode in (None, 2)
        classification = "EXPECTED_MAX_SEGMENTS"
    elif status == "STOPPED" and stop_reason in ALLOWED_METHOD_STOPS:
        accepted = returncode in (None, 2)
        classification = "PRESPECIFIED_METHOD_STOP"
    elif status == "COMPLETE":
        accepted = returncode in (None, 0)
        classification = "ROUTE_COMPLETE"
    else:
        accepted = False
        classification = "UNEXPECTED_CONTROLLER_STATE"
    return {
        "accepted_terminal_state": accepted,
        "classification": classification,
        "controller_status": status,
        "next_segment": next_segment,
        "stop_reason": stop_reason,
    }


def run_cell(cell_id: str, config_path: Path) -> dict[str, object]:
    command = [
        sys.executable,
        str(ROOT / "tools" / "titan_path_controller.py"),
        "--config",
        str(config_path),
        "--run",
    ]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    log = HERE / "generated" / f"{cell_id}.runner.log"
    log.write_text(
        "COMMAND " + " ".join(command) + "\n" + completed.stdout,
        encoding="utf-8",
        newline="\n",
    )
    classification = classify_state(config_path, completed.returncode)
    return {
        "cell_id": cell_id,
        "returncode": completed.returncode,
        "runner_log": str(log.relative_to(ROOT)).replace("\\", "/"),
        **classification,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--audit-existing", action="store_true")
    args = parser.parse_args(argv)
    if args.max_workers < 1 or args.max_workers > 8:
        raise SystemExit("--max-workers must be in [1, 8]")
    build = json.loads(BUILD_REPORT.read_text(encoding="utf-8"))
    if build.get("status") != "PREFLIGHT_OK":
        raise SystemExit("matrix build report is not PREFLIGHT_OK")
    results: dict[str, object] = {}
    if args.audit_existing:
        for cell_id, cell in build["cells"].items():
            config_path = (ROOT / cell["config"]).resolve()
            results[cell_id] = {
                "cell_id": cell_id,
                "returncode": None,
                "runner_log": str(
                    (HERE / "generated" / f"{cell_id}.runner.log").relative_to(ROOT)
                ).replace("\\", "/"),
                **classify_state(config_path, None),
            }
    else:
        futures = {}
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            for cell_id, cell in build["cells"].items():
                config_path = (ROOT / cell["config"]).resolve()
                work_dir = (config_path.parent / json.loads(
                    config_path.read_text(encoding="utf-8")
                )["work_dir"]).resolve()
                if (work_dir / "controller_state.json").exists():
                    raise SystemExit(
                        f"cell already has controller state; refusing implicit overwrite: {cell_id}"
                    )
                futures[executor.submit(run_cell, cell_id, config_path)] = cell_id
            for future in as_completed(futures):
                result = future.result()
                results[result["cell_id"]] = result
                print(
                    f"{result['cell_id']} returncode={result['returncode']} "
                    f"classification={result['classification']} "
                    f"accepted={result['accepted_terminal_state']}",
                    flush=True,
                )
    ordered = {cell_id: results[cell_id] for cell_id in build["cells"]}
    passed = all(result["accepted_terminal_state"] for result in ordered.values())
    report = {
        "schema": "titan-path-g4c1-run-report-v1",
        "status": "RUN_COMPLETE" if passed else "FAIL",
        "max_workers": args.max_workers,
        "audit_existing": args.audit_existing,
        "cell_count": len(ordered),
        "cells": ordered,
    }
    output = HERE / "generated" / "matrix_run_report.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
