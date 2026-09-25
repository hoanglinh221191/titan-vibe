#!/usr/bin/env python3
"""G2: quantify legacy spring reset versus force-continuous direction changes."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import shutil
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "test"))

from titan_path_gromacs_g1_gate import prepare_output, run, sha256, xvg_rows  # noqa: E402
from titan_path_handoff import (  # noqa: E402
    compute_force_continuous_handoff,
    mdp_updates_for_handoff,
    patch_mdp_text,
)


def copy_checkpoint_family(source: Path, destination: Path) -> None:
    destination.mkdir()
    for name in (
        "run.cpt",
        "run.edr",
        "run.gro",
        "run.log",
        "run.trr",
        "run.xtc",
        "pullx.xvg",
        "pullf.xvg",
    ):
        shutil.copy2(source / name, destination / name)


def continue_one_step(
    *,
    gmx: str,
    fixture: Path,
    branch: Path,
    mdp: Path,
) -> None:
    run(
        [
            gmx,
            "grompp",
            "-f",
            str(mdp),
            "-c",
            "run.gro",
            "-t",
            "run.cpt",
            "-p",
            str(fixture / "topol.top"),
            "-n",
            str(fixture / "index.ndx"),
            "-o",
            "next.tpr",
            "-po",
            "mdout.mdp",
            "-maxwarn",
            "0",
        ],
        cwd=branch,
        log=branch / "grompp.txt",
    )
    run(
        [
            gmx,
            "mdrun",
            "-s",
            "next.tpr",
            "-deffnm",
            "run",
            "-px",
            "pullx.xvg",
            "-pf",
            "pullf.xvg",
            "-cpo",
            "run.cpt",
            "-cpi",
            "run.cpt",
            "-reprod",
            "-ntmpi",
            "1",
            "-ntomp",
            "1",
        ],
        cwd=branch,
        log=branch / "mdrun.txt",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    fixture = ROOT / "test" / "titan_path_gromacs_g1_fixture"
    output = args.output.resolve()
    prepare_output(output)
    baseline = output / "baseline_100_steps"
    baseline.mkdir()
    version = run(
        [args.gmx, "--version"], cwd=output, log=output / "gmx_version.txt"
    )
    run(
        [
            args.gmx,
            "grompp",
            "-f",
            str(fixture / "segment1.mdp"),
            "-c",
            str(fixture / "system.gro"),
            "-p",
            str(fixture / "topol.top"),
            "-n",
            str(fixture / "index.ndx"),
            "-o",
            "run.tpr",
            "-po",
            "mdout.mdp",
            "-maxwarn",
            "0",
        ],
        cwd=baseline,
        log=baseline / "grompp.txt",
    )
    run(
        [
            args.gmx,
            "mdrun",
            "-s",
            "run.tpr",
            "-deffnm",
            "run",
            "-px",
            "pullx.xvg",
            "-pf",
            "pullf.xvg",
            "-cpo",
            "run.cpt",
            "-reprod",
            "-ntmpi",
            "1",
            "-ntomp",
            "1",
        ],
        cwd=baseline,
        log=baseline / "mdrun.txt",
    )

    base_pullx = xvg_rows(baseline / "pullx.xvg")
    base_pullf = xvg_rows(baseline / "pullf.xvg")
    fields = base_pullx[-1].split()
    if len(fields) != 12:
        raise RuntimeError("unexpected baseline pullx column count")
    time_ps = float(fields[0])
    old_coordinate = float(fields[1])
    old_reference = float(fields[2])
    reference_com = tuple(float(value) for value in fields[6:9])
    ligand_com = tuple(float(value) for value in fields[9:12])
    old_extension = old_coordinate - old_reference
    old_force = float(base_pullf[-1].split()[1])
    template = (fixture / "segment1.mdp").read_text(encoding="utf-8")

    results: list[dict[str, object]] = []
    all_pass = True
    for angle_deg in (5.0, 10.0, 20.0, 40.0):
        radians = math.radians(angle_deg)
        direction = (math.cos(radians), math.sin(radians), 0.0)
        continuous = output / f"angle_{int(angle_deg):02d}_continuous"
        legacy = output / f"angle_{int(angle_deg):02d}_legacy_reset"
        copy_checkpoint_family(baseline, continuous)
        copy_checkpoint_family(baseline, legacy)

        handoff = compute_force_continuous_handoff(
            time_ps=time_ps,
            old_coordinate_nm=old_coordinate,
            old_reference_nm=old_reference,
            new_direction=direction,
            old_direction=(1.0, 0.0, 0.0),
            reference_com_nm=reference_com,
            ligand_com_nm=ligand_com,
            rate_nm_per_ps=0.001,
            spring_k_kj_mol_nm2=100.0,
        )
        continuous_updates = mdp_updates_for_handoff(
            handoff, init_step=100, nsteps=1, simulation_part=2
        )
        continuous_mdp = continuous / "next.mdp"
        continuous_mdp.write_text(
            patch_mdp_text(template, continuous_updates),
            encoding="utf-8",
            newline="\n",
        )
        legacy_updates = dict(continuous_updates)
        legacy_updates["pull-coord1-start"] = "yes"
        legacy_updates["pull-coord1-init"] = f"{-0.001 * time_ps:.17g}"
        legacy_mdp = legacy / "next.mdp"
        legacy_mdp.write_text(
            patch_mdp_text(template, legacy_updates),
            encoding="utf-8",
            newline="\n",
        )
        continue_one_step(
            gmx=args.gmx, fixture=fixture, branch=continuous, mdp=continuous_mdp
        )
        continue_one_step(
            gmx=args.gmx, fixture=fixture, branch=legacy, mdp=legacy_mdp
        )

        continuous_x = xvg_rows(continuous / "pullx.xvg")[-1].split()
        legacy_x = xvg_rows(legacy / "pullx.xvg")[-1].split()
        continuous_force = float(xvg_rows(continuous / "pullf.xvg")[-1].split()[1])
        legacy_force = float(xvg_rows(legacy / "pullf.xvg")[-1].split()[1])
        continuous_extension = float(continuous_x[1]) - float(continuous_x[2])
        legacy_extension = float(legacy_x[1]) - float(legacy_x[2])
        continuous_jump = abs(continuous_force - old_force)
        legacy_jump = abs(legacy_force - old_force)
        passed = (
            handoff.extension_residual_nm == 0.0
            and handoff.energy_residual_kj_mol == 0.0
            and handoff.force_magnitude_residual_kj_mol_nm == 0.0
            and abs(continuous_extension - old_extension) < 0.002
            and abs(legacy_extension) < 0.002
            and continuous_jump < 0.1 * legacy_jump
            and legacy_jump > 10.0
        )
        all_pass = all_pass and passed
        results.append(
            {
                "angle_deg": angle_deg,
                "status": "PASS" if passed else "FAIL",
                "old_extension_nm": old_extension,
                "old_force_kj_mol_nm": old_force,
                "handoff": asdict(handoff),
                "continuous_first_new_row": {
                    "time_ps": float(continuous_x[0]),
                    "coordinate_nm": float(continuous_x[1]),
                    "reference_nm": float(continuous_x[2]),
                    "extension_nm": continuous_extension,
                    "force_kj_mol_nm": continuous_force,
                    "force_jump_magnitude": continuous_jump,
                },
                "legacy_first_new_row": {
                    "time_ps": float(legacy_x[0]),
                    "coordinate_nm": float(legacy_x[1]),
                    "reference_nm": float(legacy_x[2]),
                    "extension_nm": legacy_extension,
                    "force_kj_mol_nm": legacy_force,
                    "force_jump_magnitude": legacy_jump,
                },
                "observed_jump_reduction_factor": (
                    legacy_jump / continuous_jump
                    if continuous_jump > 0.0
                    else None
                ),
            }
        )

    report = {
        "schema": "titan-path-gromacs-g2-report-v1",
        "status": "PASS" if all_pass else "FAIL",
        "gromacs_version_banner": version.splitlines()[:12],
        "baseline": {
            "time_ps": time_ps,
            "coordinate_nm": old_coordinate,
            "reference_nm": old_reference,
            "extension_nm": old_extension,
            "force_kj_mol_nm": old_force,
            "checkpoint_sha256": sha256(baseline / "run.cpt"),
        },
        "angles": results,
    }
    (output / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
