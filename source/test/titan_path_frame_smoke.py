"""Synthetic numerical frame tests; no physical qualification is asserted."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"))
from titan_path_frame import (choose_core, distance_statistics, frame_metrics,
    medoid_index, training_indices, violations)


class FrameSmoke(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(619)
        self.ref = self.rng.normal(size=(24, 3))*10
        self.points = self.rng.normal(size=(5, 3))*5
        self.xyz = self.ref[None, :, :]+self.rng.normal(size=(25, 24, 3))*.02
        self.ids = np.arange(24)

    def test_streaming_moments_match_explicit(self):
        mean, std = distance_statistics(self.xyz)
        d = np.linalg.norm(self.xyz[:, :, None, :]-self.xyz[:, None, :, :], axis=3)
        np.testing.assert_allclose(mean, d.mean(axis=0), atol=1e-13)
        np.testing.assert_allclose(std, d.std(axis=0), atol=1e-13)

    def test_training_validation_disjoint_and_excludes_transient(self):
        t, v = training_indices(np.arange(0, 1001, 2))
        self.assertEqual((t[0], t[-1], v[0], v[-1]), (100, 300, 301, 500))
        self.assertFalse(set(t)&set(v))

    def test_fraction_count_and_sorted_ids(self):
        c = choose_core([self.xyz], self.ids)
        self.assertEqual(len(c["selected_indices"]), 12)
        self.assertEqual(c["selected_source_ids"], sorted(c["selected_source_ids"]))

    def test_branch_order_invariance(self):
        other = self.xyz+self.rng.normal(size=self.xyz.shape)*.01
        self.assertEqual(choose_core([self.xyz, other], self.ids),
                         {**choose_core([other, self.xyz], self.ids),
                          "branch_scores_A": choose_core([self.xyz, other], self.ids)["branch_scores_A"]})

    def test_rigid_motion_selection_invariance(self):
        q, _ = np.linalg.qr(self.rng.normal(size=(3, 3)))
        moved = self.xyz@q+[22, -8, 13]
        a, b = choose_core([self.xyz], self.ids), choose_core([moved], self.ids)
        self.assertEqual(a["selected_indices"], b["selected_indices"])
        np.testing.assert_allclose(a["worst_branch_scores_A"], b["worst_branch_scores_A"], atol=1e-13)

    def test_exact_score_tie_uses_atom_id(self):
        c = choose_core([np.repeat(self.ref[None, :, :], 5, axis=0)], self.ids[::-1])
        self.assertEqual(c["selected_source_ids"], list(range(12)))

    def test_moving_atoms_penalized(self):
        xyz = np.repeat(self.ref[None, :, :], 25, axis=0)
        xyz[:, 12:, :] += self.rng.normal(size=(25, 12, 3))*5
        c = choose_core([xyz], self.ids)
        self.assertGreater(np.mean(c["worst_branch_scores_A"][12:]), np.mean(c["worst_branch_scores_A"][:12]))

    def test_rigid_transform_zero_metric(self):
        q, _ = np.linalg.qr(self.rng.normal(size=(3, 3)))
        if np.linalg.det(q) < 0: q[:, 0] *= -1
        m = frame_metrics(self.ref, self.ref@q+[5, 2, -7], self.points,
                          np.arange(12), self.ref[:12])
        for k in ("core_rmsd_reference_A", "core_deformation_medoid_A", "split_max_route_disagreement_A", "loo_max_route_shift_A"):
            self.assertLess(m[k], 1e-10)

    def test_metrics_covary_under_common_rigid_motion(self):
        current = self.xyz[0]
        a = frame_metrics(self.ref, current, self.points, np.arange(12), self.ref[:12])
        q, _ = np.linalg.qr(self.rng.normal(size=(3, 3)))
        b = frame_metrics(self.ref@q+3, current@q+3, self.points@q+3, np.arange(12), self.ref[:12]@q+3)
        for k in a: np.testing.assert_allclose(a[k], b[k], atol=1e-11)

    def test_collinear_is_not_qualified(self):
        xyz = np.c_[np.arange(24), np.zeros((24, 2))]
        with self.assertRaises(ValueError):
            frame_metrics(xyz, xyz, self.points, np.arange(12), xyz[:12])

    def test_invalid_duplicate_ids(self):
        with self.assertRaises(ValueError): choose_core([self.xyz], np.zeros(24, dtype=int))

    def test_nonfinite_coordinate_rejected(self):
        self.xyz[0, 0, 0] = np.nan
        with self.assertRaises(ValueError): distance_statistics(self.xyz)

    def test_medoid_is_an_actual_frame(self):
        i = medoid_index(self.xyz)
        self.assertTrue(0 <= i < len(self.xyz))
        self.assertEqual(i, medoid_index(self.xyz))

    def test_bounds_do_not_turn_equality_into_failure(self):
        m = {"core_rmsd_reference_A": 1., "core_deformation_medoid_A": 1.,
             "split_max_route_disagreement_A": .25, "loo_max_route_shift_A": .25}
        p = {"core_fit_rmsd_max_A": 1., "core_fit_deformation_from_training_medoid_max_A": 1.,
             "split_core_max_route_disagreement_A": .25, "leave_one_out_max_route_shift_A": .25}
        self.assertEqual(violations(m, p), [])
        m["split_max_route_disagreement_A"] += 1e-6
        self.assertEqual(violations(m, p), ["SPLIT_CORE_FRAME_SENSITIVE"])


if __name__ == "__main__": unittest.main()
