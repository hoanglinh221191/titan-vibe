#!/usr/bin/env python3
"""Native two-particle fork gate. Only this explicit test runs tiny fixture MD."""
import argparse
import json
import math
import re
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "test"))
import titan_path_controller as c
import titan_path_isolated_segment as isolated
from titan_path_native_contract import parse_checkpoint_dump
from titan_path_handoff import compute_force_continuous_handoff, mdp_updates_for_handoff, patch_mdp_text
from titan_path_controller_native_smoke import NativeControllerSmoke


def compare_native_states(reference: str, actual: str) -> dict:
    def values(text):
        return {(kind, int(index)): float(value) for kind, index, value in re.findall(r"^\s*(x|v)\[\s*(\d+)\]\s*=\s*([-+0-9.eE]+)\s*$", text, re.MULTILINE)}
    wanted, observed = values(reference), values(actual)
    if len(wanted) != 12 or wanted != observed:
        raise RuntimeError("native fork differs from independent same-prefix reference x/v")
    return {"status": "PASS", "position_velocity_components": len(wanted), "max_absolute_difference": 0, "precision_scope": "native gmx dump printed x/v precision"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--audit-existing", action="store_true", help="Read retained completed native outputs only; never run MD")
    args = parser.parse_args()
    root = args.output.resolve()
    if args.audit_existing:
        result = isolated.read_verified_json(root / "isolated/result.json")
        for role, name in result["output_map"].items():
            if c.sha256_file(Path(name)) != result["output_hashes"][role]:
                raise RuntimeError("retained native output changed")
        ref = root / "independent_same_prefix_reference/final_checkpoint.dump"
        actual = root / "isolated_final_checkpoint.dump"
        parity = compare_native_states(ref.read_text(), actual.read_text())
        isolated.source_hashes(isolated.load_config(root / "prepared/config.json"))
        audit = {"status": "PASS_NATIVE_TWO_PARTICLE_FORK", "audit_only_no_new_md": True, "original_gate_parser_failure": "Expected vector rows, but native CPT prints flattened scalar x/v rows", "independent_reference_parity": parity, "reference_dump_sha256": c.sha256_file(ref), "actual_dump_sha256": c.sha256_file(actual), "result": result}
        isolated.write_json(root / "native_gate_post_audit.json", audit)
        print(json.dumps({key: value for key, value in audit.items() if key != "result"}, indent=2))
        return
    if root.exists():
        raise RuntimeError("native gate output already exists; retain it and use a new identity")
    root.mkdir(parents=True)
    source, prepared = root / "source", root / "prepared"
    source.mkdir()
    helper = NativeControllerSmoke()
    helper.setUp()
    try:
        shutil.copytree(helper.root, prepared)
        config, status = helper.config, helper.status
    finally:
        helper.doCleanups()
    def run(command, directory, name):
        completed = subprocess.run(command, cwd=directory, capture_output=True, text=True)
        (directory / name).write_text(completed.stdout + completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"native command failed; retained {directory / name}")
        return completed.stdout
    fixture = ROOT / "test/titan_path_gromacs_g1_fixture"
    run([args.gmx, "grompp", "-f", str(fixture / "segment1.mdp"), "-c", str(fixture / "system.gro"), "-p", str(fixture / "topol.top"), "-n", str(fixture / "index.ndx"), "-o", str(source / "run.tpr"), "-maxwarn", "0"], source, "grompp.log")
    run([args.gmx, "mdrun", "-s", str(source / "run.tpr"), "-deffnm", str(source / "run"), "-px", str(source / "pullx.xvg"), "-pf", str(source / "pullf.xvg"), "-cpo", str(source / "run.cpt"), "-ntmpi", "1", "-ntomp", "1", "-reprod"], source, "mdrun_command.log")
    seed = parse_checkpoint_dump(run([args.gmx, "dump", "-cp", str(source / "run.cpt")], source, "source_checkpoint.dump"))
    boundary = c.read_xvg_boundary(source / "pullx.xvg", time_column=0, coordinate_column=1, reference_column=2, coordinate_component_columns=[3,4,5], reference_com_columns=[6,7,8], ligand_com_columns=[9,10,11])
    for old, new in (("run.cpt", "checkpoint.cpt"), ("run.gro", "seed.gro")):
        shutil.copy2(source / old, prepared / new)
    shutil.copy2(fixture / "topol.top", prepared / "topol.top")
    shutil.copy2(fixture / "index.ndx", prepared / "index.ndx")
    run([args.gmx, "editconf", "-f", str(source / "run.gro"), "-o", str(prepared / "snapshot.pdb")], prepared, "editconf.log")
    direction = [math.cos(math.radians(20)), math.sin(math.radians(20)), 0]
    handoff = compute_force_continuous_handoff(time_ps=boundary["time_ps"], old_coordinate_nm=boundary["coordinate_nm"], old_reference_nm=boundary["reference_nm"], old_direction=[1,0,0], new_direction=direction, reference_com_nm=[0,0,0], ligand_com_nm=boundary["coordinate_components_nm"], rate_nm_per_ps=1, spring_k_kj_mol_nm2=100)
    updates = mdp_updates_for_handoff(handoff, init_step=100, nsteps=50, simulation_part=2)
    updates["dt"] = ".001"
    (prepared / "template.mdp").write_text(patch_mdp_text((fixture / "segment1.mdp").read_text(), updates))
    run([args.gmx, "grompp", "-f", str(prepared / "template.mdp"), "-c", str(prepared / "seed.gro"), "-t", str(prepared / "checkpoint.cpt"), "-p", str(prepared / "topol.top"), "-n", str(prepared / "index.ndx"), "-o", str(prepared / "initial.tpr"), "-maxwarn", "0"], prepared, "grompp.log")
    status.update(SNAPSHOT_TIME_PS=str(boundary["time_ps"]), NEXT_VECTOR_X_A=str(.5 * direction[0]), NEXT_VECTOR_Y_A=str(.5 * direction[1]), NEXT_VECTOR_ANGLE_DEG="20")
    (prepared / "status.txt").write_text("\n".join(f"{k}={v}" for k,v in status.items()) + "\n")
    (prepared / "vector.dat").write_text(f"{.5 * direction[0]:.4f} {.5 * direction[1]:.4f} 0.0000 Angle:20.0000 -- Modul:0.5000 Angstrom\n# SYNTHETIC NAV input; native mechanics only\n")
    (prepared / "boundary.json").write_text(json.dumps({"boundary": boundary, "sha256": {"checkpoint": c.sha256_file(prepared / "checkpoint.cpt")}}))
    config.update(schema=isolated.SCHEMA, cell_id="two_particle_fixture", work_dir=str(root / "isolated"), source_work_dir=str(source), gmx_command=[args.gmx], initial_direction=direction, initial_step=100, initial_simulation_part=2, dt_ps=.001, rate_nm_per_ps=1, spring_k_kj_mol_nm2=100, segment_steps=100, max_turn_angle_deg=20, mdrun_extra=["-ntmpi","1","-ntomp","1","-reprod"], environment={})
    binding = {"schema": "titan-path-initial-nav-binding-v1", "source_replay_bundle_sha256": "1" * 64, "source_replay_manifest_sha256": "2" * 64, "old_direction": [1,0,0], "source_work_dir": str(source), "source_checkpoint_output_files": seed["output_files"], "source_files_sha256": {str(path): c.sha256_file(path) for path in source.iterdir() if path.is_file()}, "sha256": {key: c.sha256_file(prepared / config["initial_nav"][key]) for key in ("status_file","vector_file","snapshot_file","boundary_record")}}
    for key in ("initial_tpr","index","topology","template_mdp","route_file"):
        binding["sha256"][key] = c.sha256_file(prepared / config[key])
    binding["sha256"]["checkpoint"] = c.sha256_file(prepared / "checkpoint.cpt")
    (prepared / "binding.json").write_text(json.dumps(binding, indent=2))
    config_path = prepared / "config.json"
    config_path.write_text(json.dumps(config, indent=2))
    result = isolated.run_segment(isolated.load_config(config_path))
    reference = root / "independent_same_prefix_reference"
    reference.mkdir()
    common = ["-deffnm", str(reference / "run"), "-px", str(reference / "pullx.xvg"), "-pf", str(reference / "pullf.xvg"), "-cpo", str(reference / "run.cpt"), "-ntmpi", "1", "-ntomp", "1", "-reprod"]
    run([args.gmx, "mdrun", "-s", str(source / "run.tpr"), *common], reference, "baseline_command.log")
    run([args.gmx, "mdrun", "-s", str(prepared / "initial.tpr"), "-cpi", str(reference / "run.cpt"), *common], reference, "continuation_command.log")
    ref_dump = run([args.gmx, "dump", "-cp", str(reference / "run.cpt")], reference, "final_checkpoint.dump")
    actual_dump = run([args.gmx, "dump", "-cp", result["output_map"]["checkpoint"]], root, "isolated_final_checkpoint.dump")
    parity = compare_native_states(ref_dump, actual_dump)
    isolated.source_hashes(isolated.load_config(config_path))
    parity["reference_checkpoint_sha256"] = c.sha256_file(reference / "run.cpt")
    report = {"schema": "titan-path-isolated-native-gate-v1", "status": "PASS_NATIVE_TWO_PARTICLE_FORK", "synthetic_nav_input": True, "b1_md_executed": False, "fixture_md_executed": True, "result": result, "independent_reference_parity": parity}
    isolated.write_json(root / "native_gate_report.json", report)
    print(json.dumps({"status": report["status"], "steps": result["steps_executed"], "source_hashes_unchanged": result["source_hashes_unchanged"], "output_map": result["output_map"]}, indent=2))


if __name__ == "__main__":
    main()
