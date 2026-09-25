#!/usr/bin/env python3
"""Build four hash-bound G4C1 calibration TPRs from two snapshots and two seeds."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "tools"))

from titan_path_handoff import patch_mdp_text  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(
    command: list[str], *, cwd: Path, log: Path, stdin: str | None = None
) -> None:
    completed = subprocess.run(
        command,
        cwd=cwd,
        input=stdin,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    log.write_text(
        "COMMAND " + " ".join(command) + "\n" + completed.stdout,
        encoding="utf-8",
        newline="\n",
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}; see {log}"
        )


def gro_atom_count(path: Path) -> int:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"truncated GRO file: {path}")
    return int(lines[1].strip())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT
        / "docs"
        / "titan_path"
        / "gates"
        / "G4C1_4w52_2snapshot_2seed_smoke_20260825"
        / "initial_states",
    )
    args = parser.parse_args(argv)
    try:
        manifest_path = HERE / "initial_state_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("schema") != "titan-path-g4c1-initial-state-manifest-v1":
            raise RuntimeError("initial-state manifest schema mismatch")
        output_root = args.output_root.resolve()
        output_root.mkdir(parents=True, exist_ok=True)
        tpr_for_whole = (ROOT / manifest["whole_molecule_tpr"]).resolve()
        topology = (ROOT / manifest["topology"]).resolve()
        index = (ROOT / manifest["index"]).resolve()
        template_path = (ROOT / manifest["template_mdp"]).resolve()
        for required in (tpr_for_whole, topology, index, template_path):
            if not required.is_file():
                raise RuntimeError(f"required input is absent: {required}")
        template = template_path.read_text(encoding="utf-8")
        report_states: dict[str, object] = {}
        expected_atoms: int | None = None
        for snapshot_id, snapshot_info in manifest["snapshots"].items():
            source = (ROOT / snapshot_info["source"]).resolve()
            if not source.is_file():
                raise RuntimeError(f"snapshot source is absent: {source}")
            snapshot_dir = output_root / snapshot_id
            snapshot_dir.mkdir(parents=True, exist_ok=True)
            whole = snapshot_dir / "whole.gro"
            run(
                [
                    args.gmx,
                    "trjconv",
                    "-s",
                    str(tpr_for_whole),
                    "-f",
                    str(source),
                    "-o",
                    str(whole),
                    "-pbc",
                    "mol",
                ],
                cwd=snapshot_dir,
                log=snapshot_dir / "trjconv.log",
                stdin="0\n",
            )
            atom_count = gro_atom_count(whole)
            if expected_atoms is None:
                expected_atoms = atom_count
            if atom_count != expected_atoms:
                raise RuntimeError("whole-molecule snapshot atom counts differ")
            for seed_id, seed_value in manifest["seeds"].items():
                state_id = f"{snapshot_id}_{seed_id}"
                state_dir = output_root / state_id
                state_dir.mkdir(parents=True, exist_ok=True)
                mdp = state_dir / "initial.mdp"
                mdp.write_text(
                    patch_mdp_text(
                        template,
                        {
                            "gen-vel": "yes",
                            "continuation": "no",
                            "gen-seed": str(seed_value),
                            "ld-seed": str(seed_value),
                            "nsteps": "1000",
                        },
                    ),
                    encoding="utf-8",
                    newline="\n",
                )
                initial_tpr = state_dir / "initial.tpr"
                run(
                    [
                        args.gmx,
                        "grompp",
                        "-f",
                        str(mdp),
                        "-c",
                        str(whole),
                        "-p",
                        str(topology),
                        "-n",
                        str(index),
                        "-o",
                        str(initial_tpr),
                        "-po",
                        str(state_dir / "mdout.mdp"),
                        "-maxwarn",
                        "0",
                    ],
                    cwd=state_dir,
                    log=state_dir / "grompp.log",
                )
                run(
                    [args.gmx, "dump", "-s", str(initial_tpr)],
                    cwd=state_dir,
                    log=state_dir / "tpr.dump",
                )
                report_states[state_id] = {
                    "snapshot_id": snapshot_id,
                    "seed_id": seed_id,
                    "seed": seed_value,
                    "source_time_ps": snapshot_info["source_time_ps"],
                    "source_sha256": sha256(source),
                    "whole_gro_sha256": sha256(whole),
                    "mdp_sha256": sha256(mdp),
                    "tpr_sha256": sha256(initial_tpr),
                    "atom_count": atom_count,
                    "initial_tpr": str(initial_tpr.relative_to(ROOT)).replace(
                        "\\", "/"
                    ),
                }
        tpr_hashes = [state["tpr_sha256"] for state in report_states.values()]
        report = {
            "schema": "titan-path-g4c1-initial-state-report-v1",
            "status": "PASS" if len(set(tpr_hashes)) == len(tpr_hashes) else "FAIL",
            "scope": manifest["scope"],
            "manifest_sha256": sha256(manifest_path),
            "whole_molecule_tpr_sha256": sha256(tpr_for_whole),
            "state_count": len(report_states),
            "distinct_tpr_hash_count": len(set(tpr_hashes)),
            "states": report_states,
        }
        output = output_root / "initial_state_report.json"
        output.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
        return 0 if report["status"] == "PASS" else 1
    except (OSError, RuntimeError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"TITAN_PATH_G4C1_INITIAL_STATE_ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
