"""No-MD design/accounting and read-only export regression tests."""
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import titan_path_benchmark_design as design
import titan_path_export_b1_equilibrium as exporter
from titan_path_route_approach import approach_text


class DesignTests(unittest.TestCase):
    def setUp(self):
        self.routes = {'straight': [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                       'curved': [[0, 0, 0], [1, 1, 0], [2, 0, 0]]}
        self.selection = design.select_routes(self.routes)
        self.rows = design.planned_rows(self.selection, {'branch_a': 'a' * 64, 'branch_b': 'b' * 64}, 'test')

    def test_approved_budget(self):
        b = design.budget_record()
        self.assertEqual((b['trajectory_count'], b['maximum_trajectory_ps'], b['maximum_total_new_md_ps']), (30, 500, 15000))
        self.assertEqual((b['md_partition'], b['md_walltime']), ('batch', '72:00:00'))
        self.assertEqual((b['auxiliary_partition'], b['auxiliary_walltime']), ('interact', '08:00:00'))
        self.assertFalse(b['automatic_retry']); self.assertFalse(b['automatic_budget_extension'])

    def test_invalid_budget(self):
        for value in (True, 0, -1, 1.5):
            with self.assertRaises(ValueError): design.Budget(max_segments=value)
        for value in (True, 0, float('nan'), float('inf')):
            with self.assertRaises(ValueError): design.Budget(segment_ps=value)

    def test_route_geometry_analytic(self):
        self.assertEqual(self.selection['relatively_direct'], 'straight')
        self.assertEqual(self.selection['relatively_curved'], 'curved')
        self.assertAlmostEqual(self.selection['metrics']['curved']['tortuosity'], 2 ** .5)
        self.assertAlmostEqual(self.selection['metrics']['curved']['maximum_turn_deg'], 90)
        self.assertFalse(self.selection['selection_uses_previous_md_outcomes'])

    def test_rotation_translation_scale_do_not_change_selection(self):
        transformed = {name: [[-3 * p[1] + 100, 3 * p[0] - 20, 3 * p[2] + 5] for p in pts]
                       for name, pts in self.routes.items()}
        other = design.select_routes(transformed)
        self.assertEqual(other['relatively_direct'], self.selection['relatively_direct'])
        self.assertEqual(other['relatively_curved'], self.selection['relatively_curved'])

    def test_degenerate_route_fails(self):
        for pts in ([[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [1, 0, 0], [0, 0, 0]]):
            with self.assertRaises(ValueError): design.route_geometry(pts)

    def test_tangent_outgoing_vertex_rotation_and_fixed_probe(self):
        import numpy as np
        points = [[0, 0, 0], [1, 0, 0], [1, 1, 0]]
        self.assertEqual(design.tangent_displacement(points, np.eye(3), 0), [.25, 0, 0])
        self.assertEqual(design.tangent_displacement(points, np.eye(3), 1), [0, .25, 0])
        self.assertEqual(design.tangent_displacement(points, np.eye(3), 2), [0, .25, 0])
        rotation = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        self.assertEqual(design.tangent_displacement(points, rotation, 0), [0, .25, 0])

    def test_invalid_tangent_inputs_fail(self):
        import numpy as np
        for progress in (-1, 3, float('nan')):
            with self.assertRaises(ValueError): design.tangent_displacement(self.routes['straight'], np.eye(3), progress)
        with self.assertRaises(ValueError): design.tangent_displacement(self.routes['straight'], np.eye(3)*2, 0)
        with self.assertRaises(ValueError): design.tangent_displacement(self.routes['straight'], np.eye(3), 0, probe_A=.5)

    def test_seed_blocks_matched_across_routes_and_methods(self):
        self.assertEqual(len(self.rows), 30)
        seeds = []
        for block in range(5):
            rows = [r for r in self.rows if r['block'] == block]
            self.assertEqual(len(rows), 6)
            self.assertEqual(len({r['rng_lock']['ld_seed'] for r in rows}), 1)
            self.assertEqual(len({r['rng_lock']['source_checkpoint_sha256'] for r in rows}), 1)
            seeds.append(rows[0]['rng_lock']['ld_seed'])
        self.assertEqual(len(set(seeds)), 5)

    def test_new_run_identity_changes_seeds(self):
        other = design.planned_rows(self.selection, {'branch_a': 'a' * 64, 'branch_b': 'b' * 64}, 'test2')
        self.assertNotEqual(other[0]['rng_lock']['ld_seed'], self.rows[0]['rng_lock']['ld_seed'])

    def test_missing_results_kept_in_denominator(self):
        result = design.account_results(self.rows, [])
        self.assertFalse(result['complete_accounting'])
        for method in design.METHODS:
            self.assertEqual(result['methods'][method]['target_count_over_planned'], [0, 10])
            self.assertEqual(result['methods'][method]['missing'], 10)

    def test_one_method_stop_does_not_label_other_arms(self):
        result = design.account_results(self.rows, [{'trajectory_id': self.rows[0]['trajectory_id'], 'status': 'METHOD_STOP'}])
        self.assertEqual(result['methods']['adaptive']['outcomes']['METHOD_STOP'], 1)
        self.assertEqual(result['methods']['fixed']['outcomes']['METHOD_STOP'], 0)
        self.assertEqual(result['methods']['tangent']['missing'], 10)

    def test_technical_failure_not_method_failure_or_success(self):
        records = [{'trajectory_id': r['trajectory_id'], 'status': 'TECHNICAL_FAILURE'} for r in self.rows]
        result = design.account_results(self.rows, records)
        self.assertTrue(result['complete_accounting'])
        self.assertFalse(result['physical_success_evidence_validated_here'])
        self.assertEqual(result['methods']['adaptive']['technical_or_missing'], 10)
        self.assertEqual(result['methods']['adaptive']['outcomes']['METHOD_STOP'], 0)

    def test_duplicate_unknown_and_nonterminal_rejected(self):
        valid = {'trajectory_id': self.rows[0]['trajectory_id'], 'status': 'TARGET_REACHED'}
        for records in ([valid, valid], [{'trajectory_id': 'other', 'status': 'TARGET_REACHED'}],
                        [{'trajectory_id': valid['trajectory_id'], 'status': 'RUNNING'}]):
            with self.assertRaises(ValueError): design.account_results(self.rows, records)


class ExportTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix='B1-EXPORT-SYNTHETIC-')
        self.addCleanup(tmp.cleanup); self.root = Path(tmp.name)
        self.source = self.root / 'source'; self.source.mkdir()
        self.data = b'SYNTHETIC CHECKPOINT, NOT MD'
        (self.source / 'a.cpt').write_bytes(self.data)
        self.expected = {'a.cpt': (hashlib.sha256(self.data).hexdigest(), len(self.data))}
        self.output = self.root / 'export.tar.gz'

    def test_round_trip_export_does_not_touch_source(self):
        before = (self.source / 'a.cpt').stat().st_mtime_ns
        result = exporter.export(self.source, self.output, self.expected)
        self.assertEqual(result['status'], 'EXPORT_READY_NO_JOB')
        with tarfile.open(self.output) as archive:
            self.assertEqual(archive.extractfile('a.cpt').read(), self.data)
            self.assertFalse(json.load(archive.extractfile('EXPORT.json'))['md_executed'])
        self.assertEqual((self.source / 'a.cpt').stat().st_mtime_ns, before)

    def test_missing_and_changed_source_rejected_before_output(self):
        for data in (b'x', b'Z' * len(self.data)):
            (self.source / 'a.cpt').write_bytes(data)
            with self.assertRaises(ValueError): exporter.export(self.source, self.output, self.expected)
            self.assertFalse(self.output.exists())
        (self.source / 'a.cpt').unlink()
        with self.assertRaises(FileNotFoundError): exporter.export(self.source, self.output, self.expected)

    def test_no_overwrite(self):
        exporter.export(self.source, self.output, self.expected)
        first = self.output.read_bytes()
        with self.assertRaises(ValueError): exporter.export(self.source, self.output, self.expected)
        self.assertEqual(self.output.read_bytes(), first)

    def test_escaping_input_rejected(self):
        for name in ('../a.cpt', '/a.cpt', 'C:/a.cpt', 'bad\\a.cpt'):
            with self.assertRaises(ValueError): exporter.read_sources(self.source, {name: self.expected['a.cpt']})


class ApproachTests(unittest.TestCase):
    def text(self):
        return ('TITAN_PATH 1\nroute_id = synthetic\nbranch_id = test\nsource = synthetic reference\n'
                'units = angstrom\nregistration = rigid\nstart_index = 0\nterminal_index = 2\n'
                'waypoint = 0 0 0 0 2\nwaypoint = 1 1 0 0 2\nwaypoint = 2 1 1 0 2\n'
                'anchor = 1 - 1 CA 0 0 0\n')

    def test_prefix_original_suffix_radius_and_identity(self):
        text, points, identity = approach_text(self.text())
        self.assertEqual(points, [[-2.,0,0],[0.,0,0],[1.,0,0],[1.,1,0]])
        self.assertEqual(identity, 'synthetic-approach2a-v1')
        self.assertIn('terminal_index = 3',text)
        self.assertIn('anchor = 1 - 1 CA 0 0 0',text)
        self.assertIn('waypoint = 1 0 0 0 2',text)

    def test_double_prefix_or_changed_length_rejected(self):
        text, _, _ = approach_text(self.text())
        with self.assertRaises(ValueError): approach_text(text)
        with self.assertRaises(ValueError): approach_text(self.text(), 3.)

    def test_malformed_indices_units_duplicate_or_degenerate_fail(self):
        for before, after in (('terminal_index = 2','terminal_index = 8'), ('units = angstrom','units = nm'),
                              ('start_index = 0','start_index = 1'), ('waypoint = 1 1 0 0 2','waypoint = 1 0 0 0 2'),
                              ('route_id = synthetic','route_id = synthetic\nroute_id = duplicate')):
            with self.assertRaises(ValueError): approach_text(self.text().replace(before,after))

    def test_three_arm_runtime_rng_does_not_use_method_in_identity(self):
        from unittest.mock import patch
        from titan_path_paired_loop import rng_for_cell
        from titan_path_tpr_pair_contract import seed_lock
        lock=seed_lock('benchmark/block_00','a'*64)
        seeds=[]
        with patch('titan_path_paired_loop.sha',return_value='a'*64):
            for method in design.METHODS:
                seeds.append(rng_for_cell(Path('synthetic'),{'cell':'blue__'+method,'rng_pair_id':'benchmark/block_00','rng_lock':lock})['ld_seed'])
        self.assertEqual(len(set(seeds)),1)
        bad={**lock,'ld_seed':lock['ld_seed']+1}
        with patch('titan_path_paired_loop.sha',return_value='a'*64),self.assertRaises(ValueError):
            rng_for_cell(Path('synthetic'),{'cell':'blue__adaptive','rng_pair_id':'benchmark/block_00','rng_lock':bad})


if __name__ == '__main__':
    unittest.main(verbosity=2)
