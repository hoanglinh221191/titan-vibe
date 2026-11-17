import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_constrained_string import StringPolicy, optimize_string, redistribute, checked_nodes


def geometry(p):
    if np.max(abs(p[:, 1:])) > .2: return None
    return [{'delta_progress': float(x[0]), 'path_distance': 0., 'minimum_clearance': 2.,
             'angle_deg': 0., 'pull_distance': float(np.linalg.norm(x))} for x in p[1:]]


def score(points):
    return [float(np.dot(p, p)) for p in points]


def rank(records):
    return min(range(len(records)), key=lambda i: records[i]['energy_key'])


class StringTests(unittest.TestCase):
    def setUp(self): self.p = np.array([[0., 0., 0.], [.2, .1, 0.], [.4, 0., 0.]])
    def test_redistribute_endpoints(self):
        q = redistribute(self.p)
        np.testing.assert_array_equal(q[[0, -1]], self.p[[0, -1]])
    def test_bad_shape(self):
        with self.assertRaises(ValueError): checked_nodes([[0, 0, 0]])
    def test_nonfinite(self):
        self.p[1, 0] = np.nan
        with self.assertRaises(ValueError): checked_nodes(self.p)
    def test_duplicate_edge(self):
        self.p[1] = self.p[0]
        with self.assertRaises(ValueError): checked_nodes(self.p)
    def test_root(self):
        self.p[0, 0] = .1
        with self.assertRaises(ValueError): checked_nodes(self.p)
    def test_policy(self):
        for kwargs in ({'iterations': 0}, {'energy_budget': 12001}, {'radii_A': (float('nan'),)}, {'iterations': True}):
            with self.assertRaises(ValueError): StringPolicy(**kwargs)
    def test_infeasible_no_energy(self):
        def fail(_): raise AssertionError('energy must not run')
        r = optimize_string(self.p, lambda _: None, fail, rank)
        self.assertEqual(r['status'], 'NO_FEASIBLE_INITIAL_STRING')
    def test_budget(self):
        r = optimize_string(self.p, geometry, score, rank, StringPolicy(energy_budget=2))
        self.assertEqual(r['status'], 'ENERGY_BUDGET_EXHAUSTED')
        self.assertEqual(r['energy_evaluations'], 2)
        np.testing.assert_array_equal(r['nodes_A'], self.p)
    def test_improves_and_no_mutation(self):
        before = self.p.copy()
        r = optimize_string(self.p, geometry, score, rank)
        self.assertLess(sum(r['final_energy_keys']), sum(r['initial_energy_keys']))
        self.assertFalse(r['md_executed']); self.assertFalse(r['md_authorized'])
        np.testing.assert_array_equal(self.p, before)
        np.testing.assert_array_equal(r['nodes_A'][0], [0, 0, 0])
    def test_bad_energy(self):
        with self.assertRaises(ValueError): optimize_string(self.p, geometry, lambda p: [np.nan]*len(p), rank)
    def test_bad_rank(self):
        with self.assertRaises(ValueError): optimize_string(self.p, geometry, score, lambda _: -1)
    def test_record_count(self):
        with self.assertRaises(ValueError): optimize_string(self.p, lambda _: [], score, rank)
    def test_reparam_never_bypasses_geometry(self):
        def fixed(p):
            return geometry(p) if np.array_equal(p, self.p) else None
        r = optimize_string(self.p, fixed, score, rank)
        np.testing.assert_array_equal(r['nodes_A'], self.p)
        self.assertTrue(all(not t['reparameterized'] for t in r['trace']))
    def test_finite_score_count(self):
        with self.assertRaises(ValueError): optimize_string(self.p, geometry, lambda _: [], rank)
    def test_deterministic(self):
        a = optimize_string(self.p, geometry, score, rank, StringPolicy(iterations=3))
        b = optimize_string(self.p, geometry, score, rank, StringPolicy(iterations=3))
        self.assertEqual(a, b)
    def test_no_unevaluated_final_nodes(self):
        r = optimize_string(self.p, geometry, score, rank, StringPolicy(energy_budget=30))
        self.assertTrue(all(x is not None and np.isfinite(x) for x in r['final_energy_keys']))
    def test_reparam_checked_after_coordinate_moves(self):
        calls = []
        def check(p):
            calls.append(p.copy())
            return geometry(p)
        r = optimize_string(self.p, check, score, rank, StringPolicy(iterations=2))
        self.assertTrue(any(np.array_equal(p, np.array(r['nodes_A'])) for p in calls))
        self.assertGreater(len(calls), 10)
    def test_short_initial_budget_no_false_success(self):
        r = optimize_string(self.p, geometry, score, rank, StringPolicy(energy_budget=1))
        self.assertEqual(r['status'], 'ENERGY_BUDGET_EXHAUSTED')
        self.assertEqual(r['energy_evaluations'], 0)
        self.assertFalse(r['md_authorized'])
    def test_coupled_constraint_rechecked(self):
        p = np.array([[0., 0., 0.], [.2, 0., 0.], [.4, 0., 0.]])
        def coupled(q):
            if np.sum(q[:, 1]) > .040000001: return None
            return geometry(q)
        r = optimize_string(p, coupled, lambda q: [-x[1] for x in q], rank, StringPolicy(iterations=2))
        self.assertLessEqual(np.sum(np.array(r['nodes_A'])[:, 1]), .040000001)
    def test_score_cache_accounts_unique_poses(self):
        received = []
        def counted(points):
            received.extend(tuple(p) for p in points)
            return score(points)
        r = optimize_string(self.p, geometry, counted, rank, StringPolicy(iterations=2))
        self.assertEqual(len(received), len(set(received)))
        self.assertEqual(r['energy_evaluations'], len(received))
    def test_callback_cannot_mutate_input(self):
        before = self.p.copy()
        def mutate(q):
            q[:] = 42.
            return None
        optimize_string(self.p, mutate, score, rank)
        np.testing.assert_array_equal(self.p, before)


if __name__ == '__main__': unittest.main()
