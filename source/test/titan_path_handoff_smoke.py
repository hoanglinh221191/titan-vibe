#!/usr/bin/env python3

import importlib.util
import math
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "tools" / "titan_path_handoff.py"
SPEC = importlib.util.spec_from_file_location("titan_path_handoff", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot load {MODULE_PATH}")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class TitanPathHandoffSmoke(unittest.TestCase):
    def test_signed_extension_energy_and_force_magnitude_are_preserved(self):
        result = MODULE.compute_force_continuous_handoff(
            time_ps=20.0,
            old_coordinate_nm=1.2,
            old_reference_nm=1.0,
            old_direction=(1.0, 0.0, 0.0),
            new_direction=(0.0, 3.0, 0.0),
            reference_com_nm=(0.0, 0.0, 0.0),
            ligand_com_nm=(0.0, 2.0, 0.0),
            rate_nm_per_ps=0.005,
            spring_k_kj_mol_nm2=600.0,
        )
        self.assertAlmostEqual(result.signed_extension_nm, 0.2, places=14)
        self.assertAlmostEqual(result.new_coordinate_nm, 2.0, places=14)
        self.assertAlmostEqual(result.new_reference_at_handoff_nm, 1.8, places=14)
        self.assertAlmostEqual(result.new_reference_init_nm, 1.7, places=14)
        self.assertAlmostEqual(result.spring_energy_before_kj_mol, 12.0, places=12)
        self.assertAlmostEqual(result.spring_energy_after_kj_mol, 12.0, places=12)
        self.assertAlmostEqual(result.scalar_force_before_kj_mol_nm, -120.0, places=12)
        self.assertAlmostEqual(result.scalar_force_after_kj_mol_nm, -120.0, places=12)
        self.assertAlmostEqual(result.direction_change_deg, 90.0, places=12)
        self.assertAlmostEqual(result.extension_residual_nm, 0.0, places=14)

    def test_negative_extension_is_not_replaced_by_magnitude(self):
        result = MODULE.compute_force_continuous_handoff(
            time_ps=5.0,
            old_coordinate_nm=0.7,
            old_reference_nm=0.9,
            new_direction=(2.0, 0.0, 0.0),
            reference_com_nm=(1.0, 0.0, 0.0),
            ligand_com_nm=(1.5, 0.0, 0.0),
            rate_nm_per_ps=0.0,
            spring_k_kj_mol_nm2=100.0,
        )
        self.assertAlmostEqual(result.signed_extension_nm, -0.2, places=14)
        self.assertAlmostEqual(result.new_reference_at_handoff_nm, 0.7, places=14)
        self.assertAlmostEqual(result.scalar_force_after_kj_mol_nm, 20.0, places=12)

    def test_legacy_reset_records_the_intended_spring_discontinuity(self):
        result = MODULE.compute_legacy_reset_handoff(
            time_ps=20.0,
            old_coordinate_nm=0.8,
            old_reference_nm=1.0,
            new_direction=(0.0, 1.0, 0.0),
            old_direction=(1.0, 0.0, 0.0),
            reference_com_nm=(0.0, 0.0, 0.0),
            ligand_com_nm=(1.0, 2.0, 0.0),
            rate_nm_per_ps=0.005,
            spring_k_kj_mol_nm2=600.0,
        )
        self.assertEqual(result.schema_version, "titan-path-legacy-spring-reset-v1")
        self.assertAlmostEqual(result.signed_extension_nm, -0.2)
        self.assertAlmostEqual(result.new_coordinate_nm, 2.0)
        self.assertAlmostEqual(result.new_reference_at_handoff_nm, 2.0)
        self.assertAlmostEqual(result.spring_energy_before_kj_mol, 12.0)
        self.assertAlmostEqual(result.spring_energy_after_kj_mol, 0.0)
        self.assertAlmostEqual(result.extension_residual_nm, 0.2)
        self.assertAlmostEqual(result.force_magnitude_residual_kj_mol_nm, -120.0)

    def test_slerp_respects_increment_bound(self):
        angle = math.radians(40.0)
        target = (math.cos(angle), math.sin(angle), 0.0)
        directions = MODULE.slerp_bounded((1.0, 0.0, 0.0), target, 15.0)
        self.assertEqual(len(directions), 3)
        previous = (1.0, 0.0, 0.0)
        for direction in directions:
            self.assertLessEqual(MODULE.angle_degrees(previous, direction), 15.0 + 1.0e-12)
            previous = direction
        self.assertLess(MODULE.angle_degrees(directions[-1], target), 1.0e-10)

    def test_invalid_direction_fails_closed(self):
        with self.assertRaises(MODULE.HandoffError):
            MODULE.compute_force_continuous_handoff(
                time_ps=0.0,
                old_coordinate_nm=0.0,
                old_reference_nm=0.0,
                new_direction=(0.0, 0.0, 0.0),
                reference_com_nm=(0.0, 0.0, 0.0),
                ligand_com_nm=(1.0, 0.0, 0.0),
                rate_nm_per_ps=0.0,
                spring_k_kj_mol_nm2=1.0,
            )

    def test_mdp_patch_sets_continuation_contract(self):
        result = MODULE.compute_force_continuous_handoff(
            time_ps=20.0,
            old_coordinate_nm=1.2,
            old_reference_nm=1.0,
            new_direction=(0.0, 1.0, 0.0),
            reference_com_nm=(0.0, 0.0, 0.0),
            ligand_com_nm=(0.0, 2.0, 0.0),
            rate_nm_per_ps=0.005,
            spring_k_kj_mol_nm2=600.0,
        )
        updates = MODULE.mdp_updates_for_handoff(
            result, init_step=10000, nsteps=10000, simulation_part=2
        )
        template = (
            "integrator = md\n"
            "continuation = no ; legacy value\n"
            "gen_vel = yes\n"
            "pull-coord1-start = yes\n"
            "pull-coord1-vec = 1 0 0\n"
            "pull-coord1-k = 123\npull-coord1-kB = 456\n"
        )
        patched = MODULE.patch_mdp_text(template, updates)
        self.assertIn("continuation = yes ; legacy value", patched)
        self.assertIn("gen_vel = no", patched)
        self.assertIn("pull-coord1-start = no", patched)
        self.assertIn("pull-coord1-dim = Y Y Y", patched)
        self.assertIn("pull-print-ref-value = yes", patched)
        self.assertIn("init-step = 10000", patched)
        self.assertIn("simulation-part = 2", patched)
        self.assertIn("pull-coord1-k = 600", patched)
        self.assertIn("pull-coord1-kB = 600", patched)
        self.assertIn("pull-coord1-rate = 0.005", patched)
        self.assertIn("pull-nstxout = 1", patched)
        self.assertIn("pull-nstfout = 1", patched)

    def test_duplicate_control_key_is_rejected(self):
        with self.assertRaises(MODULE.HandoffError):
            MODULE.patch_mdp_text(
                "continuation = no\ncontinuation = yes\n", {"continuation": "yes"}
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
