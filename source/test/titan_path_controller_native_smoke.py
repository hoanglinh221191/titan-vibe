#!/usr/bin/env python3
"""Synthetic no-MD tests of the native binary and initial NAV handoff boundary."""
from __future__ import annotations

import json
import io
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import titan_path_controller as c
from titan_path_native_contract import NativeContractError, parse_tpr_dump, validate_tpr_contract, parse_checkpoint_dump, validate_checkpoint_contract, expected_pull_membership


def native_tpr_text(expected: dict) -> str:
    """Spelling/structure follows retained gmx dump 2026.3 output."""
    group_blocks = "   pull-group 0:\n     atom: not available\n"
    for i, atoms in enumerate(expected["pull_group_atoms"], 1):
        pbcatom = expected["pull_group_pbcatoms"][i - 1]
        pbcatom = atoms[0] if pbcatom is None else pbcatom
        group_blocks += f"   pull-group {i}:\n     atom ({len(atoms)}):\n        atom[" + ",".join(str(j) for j in range(len(atoms))) + "] = {" + ",".join(map(str, atoms)) + f"}}\n     weight: not available\n     pbcatom = {pbcatom}\n"
    return f"""inputrec:
   integrator = md
   tinit = {expected['tinit_ps']}
   dt = {expected['dt_ps']}
   nsteps = {expected['nsteps']}
   init-step = {expected['init_step']}
   simulation-part = {expected['simulation_part']}
   pull = true
   pull-print-COM = true
   pull-print-ref-value = true
   pull-print-components = true
   pull-nstxout = 1
   pull-nstfout = 1
   pull-xout-average = false
   pull-fout-average = false
   pull-pbc-ref-prev-step-com = {str(expected['pull_pbc_ref_prev_step_com']).lower()}
   pull-ngroups = 3
{group_blocks}\
   pull-ncoords = 1
   pull-coord 0:
   type = umbrella
   geometry = direction
   group[0] = 1
   group[1] = 2
   dim[0]=1
   dim[1]=1
   dim[2]=1
   origin[0]=0
   origin[1]=0
   origin[2]=0
   vec[0]={expected['direction_unit'][0]}
   vec[1]={expected['direction_unit'][1]}
   vec[2]={expected['direction_unit'][2]}
   start = false
   init = {expected['reference_init_nm']}
   rate = {expected['rate_nm_per_ps']}
   k = {expected['spring_k_kj_mol_nm2']}
   kB = {expected['spring_k_kj_mol_nm2']}
   free-energy = no
topology:
   natoms = {expected['natoms']}
"""


def checkpoint_text(work: Path, step=10000) -> str:
    outputs = [str(work / name) for name in ("run.log", "pullx.xvg", "pullf.xvg", "run.edr")]
    return f"checkpoint file version = 23\n#atoms = 3\nsimulation part # = 2\nstep = {step}\nt = {step * .002:.6f}\nnumber of output files = 4\n" + "\n".join("output filename = " + name for name in outputs) + "\n"


class NativeControllerSmoke(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.work = self.root / "run"
        self.work.mkdir()
        for name in ("run.log", "pullx.xvg", "pullf.xvg", "run.edr"):
            (self.work / name).write_text("SYNTHETIC retained append output\n")
        for name in ("initial.tpr", "checkpoint.cpt", "topol.top", "route.tpath", "index.ndx", "snapshot.pdb"):
            (self.root / name).write_text("SYNTHETIC " + name)
        (self.root / "template.mdp").write_text("integrator = md\npull-coord1-k = 123\npull-ngroups = 2\npull-ncoords = 1\npull-group1-name = Protein\npull-group2-name = BENZ\n")
        (self.root / "index.ndx").write_text("[ Protein ]\n1 2\n[ BENZ ]\n3\n")
        (self.root / "snapshot.pdb").write_text("ATOM      1\nATOM      2\nHETATM    3\n")
        self.status = {
            "SCHEMA": "TITAN_PATH_STATUS_V1", "STATUS": "CONTINUE", "REASON": "OK", "REASON_DETAIL": "NONE",
            "ROUTE_KIND": "EXPLICIT", "ROUTE_ID": "test", "REFERENCE_SHA256": "0" * 64,
            "ROUTE_FILE_SHA256": c.sha256_file(self.root / "route.tpath"), "MODE": "Advance", "MAX_TURN_ANGLE_DEG": "20",
            "USED_OLD_VECTOR": "1", "TARGET_CLEARANCE_VECTOR": "NEXT_VECTOR", "NEXT_VECTOR_ANGLE_DEG": "0",
            "NEXT_VECTOR_X_A": "0.5", "NEXT_VECTOR_Y_A": "0", "NEXT_VECTOR_Z_A": "0", "SNAPSHOT_TIME_PS": "20",
            "PATH_PROGRESS_A": "1", "TARGET_PROGRESS_A": "1.5", "APPLIED_PROGRESS_A": "1.5", "APPLIED_DELTA_PROGRESS_A": "0.5", "APPLIED_PATH_DISTANCE_A": "0.5",
            "VALID_ARCHIVE_SIZE": "1", "FORWARD_ARCHIVE_SIZE": "1", "EPSILON_ELIGIBLE_SIZE": "1",
            "NAVIGATION_ENABLED": "1", "NAVIGATION_COST_KIND": "ALGORITHMIC_DIMENSIONLESS", "NAVIGATION_COST_EPSILON": ".02", "TARGET_NAVIGATION_TUBE_COST": ".01", "TARGET_NAVIGATION_BACKTRACK_COST": "0", "TARGET_NAVIGATION_COST": ".01", "NAVIGATION_ELIGIBLE_SIZE": "1", "ENERGY_ELIGIBLE_MIN_NAVIGATION_COST": "0",
            "FIXED_PROBE_ENABLED": "1", "CONFIGURED_PROBE_LENGTH_A": "1", "CONFIGURED_PROBE_LENGTHS_A": "1,0.5,0.25", "NEXT_REFERENCE_DISPLACEMENT_A": "0.5",
            "CLEARANCE_MODE": "swept-linear-center", "CLEARANCE_SCOPE": "CONTINUOUS_LINEAR_TRANSLATION",
            "CLEARANCE_COORDINATE_FRAME": "UPSTREAM_PBC_RECONSTRUCTED_CARTESIAN", "CLEARANCE_PERIODIC_DISTANCE_POLICY": "DIRECT_IN_RECONSTRUCTED_FRAME", "CLEARANCE_PBC_PREPROCESSOR": "pbcn",
            "SWEPT_CLEARANCE_PROBE_EVALUATIONS": "1", "SWEPT_CLEARANCE_PAIR_CHECKS": "4", "SWEPT_CLEARANCE_REJECTIONS": "0", "SWEPT_INTERIOR_ONLY_REJECTIONS": "0",
            "TARGET_SWEPT_MIN_CLEARANCE_A": "1.6", "TARGET_ENDPOINT_MIN_CLEARANCE_A": "1.6", "TARGET_MIN_CLEARANCE_A": "1.6", "MIN_CLEARANCE_THRESHOLD_A": "1.5", "TARGET_SWEPT_MIN_FRACTION": "1", "TARGET_SWEPT_MIN_COMPLETE": "1",
        }
        (self.root / "vector.dat").write_text("0.5 0 0\n")
        self.write_status()
        self.boundary_record = {"boundary": {"time_ps": 20, "coordinate_nm": .1, "reference_nm": .2, "coordinate_components_nm": [.1, 0, 0]}, "sha256": {"checkpoint": c.sha256_file(self.root / "checkpoint.cpt")}}
        (self.root / "boundary.json").write_text(json.dumps(self.boundary_record))
        self.config = {"schema": c.CONFIG_SCHEMA, "work_dir": "run", "gmx_command": ["SYNTHETIC_GMX"], "titan_command": ["SYNTHETIC_ADAPTER", "--profile", "nav-p2", "--previous-status", "{previous_status}", "--previous-status-sha256", "{previous_status_sha256}", "--progress-window-A", "{route_continuity_window_A}"],
                       "template_mdp": "template.mdp", "topology": "topol.top", "route_file": "route.tpath", "index": "index.ndx", "initial_tpr": "initial.tpr", "initial_checkpoint": "checkpoint.cpt",
                       "initial_direction": [1, 0, 0], "initial_step": 10000, "initial_simulation_part": 2, "dt_ps": .002, "segment_steps": 10000, "max_segments": 1, "rate_nm_per_ps": .005, "spring_k_kj_mol_nm2": 600,
                       "expected_route_id": "test", "expected_reference_sha256": "0" * 64, "expected_clearance_mode": "swept-linear-center", "reference_displacement_contract": "match_titan_probe", "allowed_probe_lengths_A": [1, .5, .25], "recovery_max_turn_angle_deg": 60, "route_continuity_window_A": 1, "checkpoint_output_policy": "same_prefix",
                       "pullx_columns": {"time": 0, "coordinate": 1, "reference": 2, "coordinate_components": [3, 4, 5], "reference_com": [6, 7, 8], "ligand_com": [9, 10, 11]},
                       "initial_nav": {"status_file": "status.txt", "vector_file": "vector.dat", "snapshot_file": "snapshot.pdb", "boundary_record": "boundary.json", "binding_file": "binding.json", "old_direction": [1, 0, 0]}}
        self.rebind()

    def write_status(self):
        (self.root / "status.txt").write_text("\n".join(f"{k}={v}" for k, v in self.status.items()) + "\n")

    def rebind(self):
        binding = {"schema": "titan-path-initial-nav-binding-v1", "source_replay_bundle_sha256": "1" * 64, "source_replay_manifest_sha256": "2" * 64, "old_direction": self.config["initial_nav"]["old_direction"], "sha256": {key: c.sha256_file(self.root / self.config["initial_nav"][key]) for key in ("status_file", "vector_file", "snapshot_file", "boundary_record")}}
        binding["sha256"]["checkpoint"] = c.sha256_file(self.root / "checkpoint.cpt")
        for key in ("initial_tpr", "index", "topology", "template_mdp", "route_file"):
            binding["sha256"][key] = c.sha256_file(self.root / self.config[key])
        (self.root / "binding.json").write_text(json.dumps(binding))

    def load(self):
        path = self.root / "config.json"
        path.write_text(json.dumps(self.config))
        return c.load_and_validate_config(path)

    def test_initial_half_probe_and_native_semantics(self):
        validated = c.validate_initial_nav(self.load())
        expected = validated["expected_tpr"]
        self.assertEqual(expected["nsteps"], 5000)
        self.assertAlmostEqual(expected["reference_init_nm"], .1)
        self.assertAlmostEqual(validated["handoff"]["signed_extension_nm"], -.1)
        validate_tpr_contract(parse_tpr_dump(native_tpr_text(expected)), expected)
        for key, wrong in (("dt_ps", .004), ("nsteps", 10000), ("init_step", 0), ("simulation_part", 1), ("direction_unit", [0, 1, 0]), ("reference_init_nm", 0), ("spring_k_kj_mol_nm2", 123), ("rate_nm_per_ps", .01), ("tinit_ps", 20)):
            with self.subTest(field=key), self.assertRaises(NativeContractError):
                validate_tpr_contract(parse_tpr_dump(native_tpr_text({**expected, key: wrong})), expected)
        for before, after in (("start = false", "start = true"), ("geometry = direction", "geometry = distance"), ("pull-nstxout = 1", "pull-nstxout = 5"), ("pull-ncoords = 1", "pull-ncoords = 2"), ("kB = 600", "kB = 12")):
            with self.subTest(change=after), self.assertRaises(NativeContractError):
                validate_tpr_contract(parse_tpr_dump(native_tpr_text(expected).replace(before, after)), expected)

    def test_ledger_and_selected_vector_drift_rejected(self):
        loaded = self.load()
        (self.root / "snapshot.pdb").write_text("ATOM  different1\nATOM  different2\nHETATM different3\n")
        with self.assertRaisesRegex(c.ControllerError, "snapshot_file"):
            c.validate_initial_nav(loaded)
        self.rebind()
        self.config["initial_direction"] = [0, 1, 0]
        with self.assertRaisesRegex(c.ControllerError, "initial_direction"):
            c.validate_initial_nav(self.load())
        self.config["initial_direction"] = [1, 0, 0]
        self.config["initial_step"] = 0
        with self.assertRaisesRegex(c.ControllerError, "initial_step"):
            c.validate_initial_nav(self.load())

    def test_initial_recovery_uses_60_degree_limit(self):
        direction = [math.cos(math.pi / 3), math.sin(math.pi / 3), 0]
        self.config["initial_direction"] = direction
        self.status.update(MODE="Recovery", MAX_TURN_ANGLE_DEG="60", NEXT_VECTOR_ANGLE_DEG="60", NEXT_VECTOR_X_A=str(.5 * direction[0]), NEXT_VECTOR_Y_A=str(.5 * direction[1]), PATH_DISTANCE_A="4", APPLIED_PATH_DISTANCE_A="3.5")
        (self.root / "vector.dat").write_text(f"{.5 * direction[0]} {.5 * direction[1]} 0\n")
        self.write_status(); self.rebind()
        validated = c.validate_initial_nav(self.load())
        self.assertAlmostEqual(validated["handoff"]["direction_change_deg"], 60)
        self.status["MODE"] = "Advance"
        self.write_status(); self.rebind()
        with self.assertRaisesRegex(c.ControllerError, "turn"):
            c.validate_initial_nav(self.load())

    def test_checkpoint_prefix_and_time_are_required(self):
        observed = parse_checkpoint_dump(checkpoint_text(self.work))
        validate_checkpoint_contract(observed, work_dir=self.work, step=10000, time_ps=20)
        for changed in ({"step": 9999}, {"time_ps": 19.99}, {"output_files": [str(self.root / "other" / Path(name).name) for name in observed["output_files"]]}, {"output_files": [Path(name).name for name in observed["output_files"]]}):
            with self.assertRaises(NativeContractError):
                validate_checkpoint_contract({**observed, **changed}, work_dir=self.work, step=10000, time_ps=20)
        with self.assertRaisesRegex(NativeContractError, "atom count"):
            validate_checkpoint_contract(observed, work_dir=self.work, step=10000, time_ps=20, natoms=4)

    def test_actual_native_atom_membership_and_count_must_match(self):
        expected = c.validate_initial_nav(self.load())["expected_tpr"]
        self.assertEqual(expected["pull_group_atoms"], [[0, 1], [2]])
        for wrong in ({"pull_group_atoms": [[0], [1, 2]]}, {"natoms": 4}):
            with self.assertRaises(NativeContractError):
                validate_tpr_contract(parse_tpr_dump(native_tpr_text({**expected, **wrong})), expected)
        compressed = native_tpr_text(expected).replace("atom[0,1] = {0,1}", "atom[0,...,1] = {0,...,1}")
        validate_tpr_contract(parse_tpr_dump(compressed), expected)
        for before, after in (("{0,...,1}", "{0,0}"), ("pull-group 1:", "pull-group 2:"), ("{0,...,1}", "{0,...,999999999}"), ("atom (2):", "atom (3):"), ("weight: not available", "weight (2):\n weight[0] = 1")):
            with self.subTest(change=after), self.assertRaises(NativeContractError):
                parse_tpr_dump(compressed.replace(before, after, 1))

    def test_ndx_duplicate_and_selected_only_snapshot_are_rejected(self):
        for content in ("[ Protein ]\n1 1\n[ BENZ ]\n3\n", "[ Protein ]\n1 2\n[ Protein ]\n3\n", "[ Protein ]\n1 2\n[ BENZ ]\n2\n"):
            (self.root / "index.ndx").write_text(content)
            with self.assertRaises(NativeContractError):
                expected_pull_membership(self.root / "template.mdp", self.root / "index.ndx", self.root / "snapshot.pdb")
        (self.root / "index.ndx").write_text("[ Protein ]\n1 2\n[ BENZ ]\n3\n")
        (self.root / "snapshot.pdb").write_text("ATOM      1\nATOM      2\n")
        with self.assertRaisesRegex(NativeContractError, "full-system"):
            expected_pull_membership(self.root / "template.mdp", self.root / "index.ndx", self.root / "snapshot.pdb")

    def test_producer_binding_freezes_initial_tpr_and_engine_inputs(self):
        loaded = self.load()
        for key in ("initial_tpr", "index", "topology", "template_mdp", "route_file"):
            path = loaded[key]
            original = path.read_bytes()
            path.write_bytes(original + b"\nSYNTHETIC alteration")
            with self.subTest(field=key), self.assertRaisesRegex(c.ControllerError, key):
                c.validate_initial_nav(loaded)
            path.write_bytes(original)

    def test_cli_preflight_audits_native_inputs_without_state_or_md(self):
        loaded = self.load()
        expected = c.validate_initial_nav(loaded)["expected_tpr"]
        calls = []
        def fake(command, *, log, **kwargs):
            calls.append(command)
            self.assertIn("dump", command)
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(native_tpr_text(expected) if "-s" in command else checkpoint_text(self.work))
            return 0
        output = io.StringIO()
        with patch.object(c, "_run_command", side_effect=fake), patch("sys.stdout", output):
            self.assertEqual(c.main(["--config", str(self.root / "config.json"), "--preflight"]), 0)
        report = json.loads(output.getvalue())
        self.assertEqual(report["status"], "PREFLIGHT_OK")
        self.assertEqual(report["initial_nav_validation"]["native"]["observed_tpr"]["nsteps"], 5000)
        self.assertEqual(len(calls), 2)
        self.assertFalse((self.work / "controller_state.json").exists())

    def test_explicit_pbc_atom_and_weight_policy_are_bound(self):
        template = self.root / "template.mdp"
        template.write_text(template.read_text() + "pull-group1-pbcatom = 2\npull-pbc-ref-prev-step-com = yes\n")
        self.rebind()
        expected = c.validate_initial_nav(self.load())["expected_tpr"]
        self.assertEqual(expected["pull_group_pbcatoms"], [1, None])
        validate_tpr_contract(parse_tpr_dump(native_tpr_text(expected)), expected)
        for before, after in (("pbcatom = 1", "pbcatom = 0"), ("pull-pbc-ref-prev-step-com = true", "pull-pbc-ref-prev-step-com = false")):
            with self.assertRaises(NativeContractError):
                validate_tpr_contract(parse_tpr_dump(native_tpr_text(expected).replace(before, after, 1)), expected)
        template.write_text(template.read_text() + "pull-group1-weights = 1 1\n")
        with self.assertRaisesRegex(c.ControllerError, "custom MDP pull weights"):
            self.load()

    def test_invalid_native_tpr_never_reaches_mdrun(self):
        loaded = self.load()
        expected = c.validate_initial_nav(loaded)["expected_tpr"]
        calls = []
        def fake(command, *, log, **kwargs):
            calls.append(command)
            self.assertNotIn("mdrun", command)
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(native_tpr_text({**expected, "nsteps": 10000}))
            return 0
        with patch.object(c, "_run_command", side_effect=fake), self.assertRaisesRegex(c.ControllerError, "nsteps"):
            c.run_controller(loaded, resume=False)
        self.assertEqual(len(calls), 1)
        self.assertFalse((self.work / "controller_state.json").exists())

    def test_native_prefix_failure_never_reaches_mdrun(self):
        loaded = self.load()
        expected = c.validate_initial_nav(loaded)["expected_tpr"]
        calls = []
        def fake(command, *, log, **kwargs):
            calls.append(command)
            self.assertNotIn("mdrun", command)
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(native_tpr_text(expected) if "-s" in command else checkpoint_text(self.root / "other"))
            return 0
        with patch.object(c, "_run_command", side_effect=fake), self.assertRaisesRegex(c.ControllerError, "checkpoint output identity"):
            c.run_controller(loaded, resume=False)
        self.assertEqual(len(calls), 2)

    def test_strict_config_cannot_disable_native_requirements(self):
        for key, wrong in (("initial_nav", None), ("route_continuity_window_A", 0), ("checkpoint_output_policy", "noappend"), ("mdrun_extra", ["-noappend"]), ("titan_command", ["legacy_adapter"])):
            original = self.config.get(key)
            self.config[key] = wrong
            with self.subTest(field=key), self.assertRaises(c.ControllerError):
                self.load()
            if original is None:
                del self.config[key]
            else:
                self.config[key] = original

    def test_mocked_strict_two_segment_cycle_keeps_probe_and_previous_status(self):
        self.config["max_segments"] = 3
        loaded = self.load()
        initial_expected = c.validate_initial_nav(loaded)["expected_tpr"]
        md_steps = []
        previous_paths = []
        def value(command, flag):
            return Path(command[command.index(flag) + 1])
        def expected_for(tpr):
            return initial_expected if tpr == loaded["initial_tpr"] else json.loads(tpr.read_text())
        def fake(command, *, cwd, log, **kwargs):
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text("SYNTHETIC " + " ".join(command) + "\n")
            if "dump" in command:
                if "-s" in command:
                    log.write_text(native_tpr_text(expected_for(value(command, "-s"))))
                else:
                    checkpoint = value(command, "-cp")
                    step = 10000 if checkpoint == loaded["initial_checkpoint"] else int(checkpoint.read_text())
                    log.write_text(checkpoint_text(self.work, step))
            elif "mdrun" in command:
                expected = expected_for(value(command, "-s"))
                md_steps.append(expected["nsteps"])
                end = expected["init_step"] + expected["nsteps"]
                time = end * .002
                ref = expected["reference_init_nm"] + .005 * time
                coord = ref - .1
                value(command, "-cpo").write_text(str(end))
                (self.work / "run.gro").write_text("SYNTHETIC gro")
                value(command, "-px").write_text(f"{time} {coord} {ref} {coord} 0 0 0 0 0 {coord} 0 0\n")
            elif "grompp" in command:
                updates = dict(line.split("=", 1) for line in value(command, "-f").read_text().splitlines() if "=" in line)
                updates = {key.strip().lower(): val.strip() for key, val in updates.items()}
                self.assertEqual(float(updates["pull-coord1-k"]), 600)
                expected = {**loaded["native_pull_membership"], "dt_ps": float(updates["dt"]), "tinit_ps": float(updates["tinit"]), "init_step": int(updates["init-step"]), "nsteps": int(updates["nsteps"]), "simulation_part": int(updates["simulation-part"]), "direction_unit": [float(x) for x in updates["pull-coord1-vec"].split()], "reference_init_nm": float(updates["pull-coord1-init"]), "spring_k_kj_mol_nm2": float(updates["pull-coord1-k"]), "rate_nm_per_ps": float(updates["pull-coord1-rate"])}
                value(command, "-o").write_text(json.dumps(expected))
            else:
                prior = value(command, "--previous-status")
                previous_paths.append(prior)
                self.assertEqual(command[command.index("--previous-status-sha256") + 1], c.sha256_file(prior))
                self.assertEqual(command[command.index("--progress-window-A") + 1], "1.0")
                status = {**self.status, "SNAPSHOT_TIME_PS": "30" if len(md_steps) == 1 else "35"}
                if len(md_steps) == 1:
                    status.update(NEXT_VECTOR_X_A=".25", NEXT_REFERENCE_DISPLACEMENT_A=".25", PATH_PROGRESS_A="1.5", TARGET_PROGRESS_A="1.75", APPLIED_PROGRESS_A="1.75", APPLIED_DELTA_PROGRESS_A=".25", ROUTE_TARGET_VECTOR_X_A="1", ROUTE_TARGET_VECTOR_Y_A="0", ROUTE_TARGET_VECTOR_Z_A="0", TURN_CAPPED="0")
                else:
                    status.update(STATUS="STOP", REASON="NO_FORWARD_PROGRESS_CANDIDATE", REASON_DETAIL="FORWARD_PROGRESS_EXHAUSTED", VALID_ARCHIVE_SIZE="236", FORWARD_ARCHIVE_SIZE="0", FORWARD_PROGRESS_REJECTIONS="236", EPSILON_ELIGIBLE_SIZE="0", NAVIGATION_ELIGIBLE_SIZE="0")
                (cwd / "pull_direction_status.txt").write_text("\n".join(f"{k}={v}" for k, v in status.items()) + "\n")
            return 0.0
        with patch.object(c, "_run_command", side_effect=fake):
            state = c.run_controller(loaded, resume=False)
        self.assertEqual(md_steps, [5000, 2500])
        self.assertEqual(state["status"], "STOPPED")
        self.assertEqual(state["stop_reason_detail"], "FORWARD_PROGRESS_EXHAUSTED")
        self.assertEqual(previous_paths, [self.root / "status.txt", self.work / "segments/segment_0000/pull_direction_status.txt"])
        with patch.object(c, "_run_command", side_effect=AssertionError("terminal must not execute")):
            self.assertEqual(c.run_controller(loaded, resume=True)["status"], "STOPPED")

    def test_interrupted_checkpoint_range_and_prepared_tpr_hash(self):
        loaded = self.load()
        expected = c.validate_initial_nav(loaded)["expected_tpr"]
        step = 12500
        def fake_dump(command, *, log, **kwargs):
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(native_tpr_text(expected) if "-s" in command else checkpoint_text(self.work, step))
            return 0.0
        with patch.object(c, "_run_command", side_effect=fake_dump):
            c.audit_native_inputs(loaded, tpr=loaded["initial_tpr"], expected=expected, checkpoint=loaded["initial_checkpoint"], log_dir=self.work / "audit", allow_interrupted_checkpoint=True)
            step = 15001
            with self.assertRaisesRegex(c.ControllerError, "outside"):
                c.audit_native_inputs(loaded, tpr=loaded["initial_tpr"], expected=expected, checkpoint=loaded["initial_checkpoint"], log_dir=self.work / "audit", allow_interrupted_checkpoint=True)
        prepared = self.work / "changed.tpr"
        prepared.write_text("SYNTHETIC original prepared TPR")
        state = c._initial_state(loaded)
        state.update(current_tpr=str(prepared), current_tpr_expected=expected, current_tpr_sha256=c.sha256_file(prepared))
        c.atomic_write_json(self.work / "controller_state.json", state)
        prepared.write_text("SYNTHETIC changed prepared TPR")
        with patch.object(c, "_run_command", side_effect=AssertionError("changed TPR must not execute")), self.assertRaisesRegex(c.ControllerError, "prepared TPR changed"):
            c.run_controller(loaded, resume=True)

    def test_interrupted_first_segment_adopts_newer_bound_checkpoint(self):
        loaded = self.load()
        expected = c.validate_initial_nav(loaded)["expected_tpr"]
        state = c._initial_state(loaded)
        state.update(current_tpr_expected=expected, current_tpr_sha256=c.sha256_file(loaded["initial_tpr"]))
        c.atomic_write_json(self.work / "controller_state.json", state)
        (self.work / "run.cpt").write_text("12500")
        class MockReachedMD(RuntimeError):
            pass
        def fake(command, *, log, **kwargs):
            if "mdrun" in command:
                self.assertEqual(command[command.index("-cpi") + 1], str(self.work / "run.cpt"))
                raise MockReachedMD("mock stopped before execution")
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(native_tpr_text(expected) if "-s" in command else checkpoint_text(self.work, 12500))
            return 0
        with patch.object(c, "_run_command", side_effect=fake), self.assertRaises(MockReachedMD):
            c.run_controller(loaded, resume=True)
        adopted = json.loads((self.work / "controller_state.json").read_text())
        self.assertEqual(adopted["checkpoint"], str(self.work / "run.cpt"))
        self.assertEqual(adopted["checkpoint_adoptions"][0]["reason"], "INTERRUPTED_SEGMENT0_BOUND_CHECKPOINT_RECOVERY")

    def test_graceful_short_md_does_not_advance_or_call_titan(self):
        loaded = self.load()
        expected = c.validate_initial_nav(loaded)["expected_tpr"]
        def fake(command, *, log, **kwargs):
            log.parent.mkdir(parents=True, exist_ok=True)
            if "mdrun" in command:
                (self.work / "run.cpt").write_text("12500")
                (self.work / "run.gro").write_text("SYNTHETIC shortened output")
                (self.work / "pullx.xvg").write_text("25 .125 .225 .125 0 0 0 0 0 .125 0 0\n")
                log.write_text("SYNTHETIC normal exit before requested endpoint")
            elif "dump" in command:
                log.write_text(native_tpr_text(expected) if "-s" in command else checkpoint_text(self.work))
            else:
                self.fail("shortened MD must not run TITAN or grompp")
            return 0
        with patch.object(c, "_run_command", side_effect=fake), self.assertRaisesRegex(c.ControllerError, "exact NAV segment boundary"):
            c.run_controller(loaded, resume=False)
        state = json.loads((self.work / "controller_state.json").read_text())
        self.assertEqual(state["next_segment"], 0)
        self.assertEqual(state["status"], "READY")
        self.assertEqual((self.work / "run.cpt").read_text(), "12500")


if __name__ == "__main__":
    unittest.main(verbosity=2)
