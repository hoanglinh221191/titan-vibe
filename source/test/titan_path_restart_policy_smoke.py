"""No-MD policy, geometry, provenance and mutation gates."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import titan_path_restart_policy as policy
import titan_path_benchmark_runtime as runtime

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/'.codex_tmp/b1_r2_native_restart_20260914'


class PolicyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); (self.root/'inputs').mkdir()
        self.config = policy.b1_policy()
        self.cell = {'source_branch': 'branch_a', 'rng_lock': {'source_checkpoint_sha256':
                     self.config['branches']['branch_a']['equilibrium_checkpoint_sha256']}}

    def enable(self):
        path = self.root/policy.POLICY_PATH; path.write_text(json.dumps(self.config))
        (self.root/'PROTOCOL.json').write_text(json.dumps({'restart_policy_sha256': policy.sha(path)}))
        return policy.for_cell(self.root, self.cell)

    def test_historical_mode_is_unchanged(self):
        self.assertIsNone(policy.for_cell(self.root, self.cell))
        self.assertEqual(policy.mdp_updates(None), {})

    def test_explicit_branch_anchors(self):
        p = self.enable()
        self.assertEqual(policy.mdp_updates(p)['pull-group1-pbcatom'], '1605')
        self.assertEqual(policy.mdp_updates(p)['pull-group2-pbcatom'], '2609')
        self.assertEqual(policy.mdp_updates(p)['pull-pbc-ref-prev-step-com'], 'no')

    def test_missing_declared_policy_rejected(self):
        self.enable(); (self.root/policy.POLICY_PATH).unlink()
        with self.assertRaisesRegex(ValueError, 'BINDING'):
            policy.for_cell(self.root, self.cell)

    def test_unbound_file_rejected(self):
        (self.root/policy.POLICY_PATH).write_text(json.dumps(self.config))
        with self.assertRaisesRegex(ValueError, 'BINDING'):
            policy.for_cell(self.root, self.cell)

    def test_changed_policy_hash_rejected(self):
        self.enable(); (self.root/policy.POLICY_PATH).write_text('{}')
        with self.assertRaises(ValueError):
            policy.for_cell(self.root, self.cell)

    def test_relaxed_margin_even_rehashed_rejected(self):
        self.config['pbc_margin_nm'] = 0
        with self.assertRaisesRegex(ValueError, 'SCHEMA_OR_VALUES'):
            self.enable()

    def test_wrong_branch_checkpoint_rejected(self):
        self.enable(); self.cell['rng_lock']['source_checkpoint_sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'BRANCH_OR_SEED'):
            policy.for_cell(self.root, self.cell)

    def test_outside_group_anchor_rejected(self):
        with self.assertRaisesRegex(ValueError, 'MEMBERSHIP'):
            policy.group_geometry(np.zeros((3, 3)), [0, 1], 2, np.eye(3)*8, .25, .0001)

    def test_safe_central_reference(self):
        result = policy.group_geometry(np.array([[0., 0, 0], [2., 0, 0]]), [0, 1], 0, np.eye(3)*8, .25, .0001)
        self.assertEqual(result['remaining_margin_nm'], 2.)

    def test_marginal_reference_holds(self):
        with self.assertRaisesRegex(ValueError, 'MARGIN_HOLD'):
            policy.group_geometry(np.array([[0., 0, 0], [3.8, 0, 0]]), [0, 1], 0, np.eye(3)*8, .25, .0001)

    def test_individual_periodic_shift_preserves_guard(self):
        box = np.array([[8., 0, 0], [0, 8., 0], [4., 4., 5.656854]])
        x = np.array([[1., 1, 1], [2., 1, 1]])
        a = policy.group_geometry(x, [0, 1], 0, box, .25, .0001)
        x[1] += 3*box[2]-2*box[0]
        b = policy.group_geometry(x, [0, 1], 0, box, .25, .0001)
        self.assertAlmostEqual(a['max_radius_nm'], b['max_radius_nm'], places=12)

    def test_bad_box_rejected(self):
        for box in (np.zeros((3, 3)), np.diag([1., 1., -1.]), np.full((3, 3), np.nan)):
            with self.assertRaises(ValueError):
                policy.box_geometry(box)

    def test_native_old_default_reference_not_accepted_as_fixed(self):
        p = self.enable(); native = (RAW/'outputs/trajectories/orange__block_00__adaptive/segment_002/prepared/tpr.dump').read_text()
        with self.assertRaisesRegex(ValueError, 'ANCHOR_CHANGED'):
            policy.geometry_guard(native, p)

    def test_rehashing_anchor_change_rejected(self):
        self.config['branches']['branch_a']['pbcatoms_0based'][0] = 1301
        with self.assertRaises(ValueError):
            self.enable()

    def test_previous_com_mode_switch_rejected(self):
        self.config['previous_com_mode'] = True
        with self.assertRaises(ValueError):
            self.enable()


class NativeVerletTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        source = RAW/'outputs/trajectories/blue__block_03__tangent/segment_042/prepared'
        for name in ('tpr.dump', 'segment.tpr', 'grompp.log'):
            (self.path/name).write_bytes((source/name).read_bytes())

    def mutate(self, old, new):
        p = self.path/'tpr.dump'; text = p.read_text(); self.assertIn(old, text); p.write_text(text.replace(old, new, 1))

    def test_native_computed_radius_evidence(self):
        self.assertEqual(policy.neighbor_proof(self.path)['rlist_nm'], 1.222)

    def test_radius_without_matching_producer_rejected(self):
        self.mutate('rlist                          = 1.222', 'rlist                          = 1.223')
        with self.assertRaisesRegex(ValueError, 'PRODUCER_MISMATCH'):
            policy.neighbor_proof(self.path)

    def test_missing_producer_rejected(self):
        (self.path/'grompp.log').unlink()
        with self.assertRaisesRegex(ValueError, 'LOG_ABSENT'):
            policy.neighbor_proof(self.path)

    def test_ambiguous_producer_rejected(self):
        (self.path/'tpr.grompp.log').write_bytes((self.path/'grompp.log').read_bytes())
        with self.assertRaisesRegex(ValueError, 'AMBIGUOUS'):
            policy.neighbor_proof(self.path)

    def test_changed_physical_cutoff_rejected(self):
        self.mutate('rcoulomb                       = 1.2', 'rcoulomb                       = 1.1')
        with self.assertRaisesRegex(ValueError, 'FIXED_FIELD'):
            policy.neighbor_proof(self.path)

    def test_looser_buffer_tolerance_rejected(self):
        self.mutate('verlet-buffer-tolerance        = 0.005', 'verlet-buffer-tolerance        = 0.01')
        with self.assertRaisesRegex(ValueError, 'FIXED_FIELD'):
            policy.neighbor_proof(self.path)

    def test_nstlist_change_rejected(self):
        self.mutate('nstlist                        = 20', 'nstlist                        = 80')
        with self.assertRaisesRegex(ValueError, 'NSTLIST'):
            policy.neighbor_proof(self.path)

    def test_producer_tolerance_change_rejected(self):
        p = self.path/'grompp.log'; p.write_text(p.read_text().replace('0.005 kJ/mol/ps', '0.05 kJ/mol/ps'))
        with self.assertRaisesRegex(ValueError, 'PRODUCER_TOLERANCE'):
            policy.neighbor_proof(self.path)

    def test_duplicate_producer_line_rejected(self):
        p = self.path/'grompp.log'; p.write_text(p.read_text()*2)
        with self.assertRaises(ValueError):
            policy.neighbor_proof(self.path)


if __name__ == '__main__':
    unittest.main()
