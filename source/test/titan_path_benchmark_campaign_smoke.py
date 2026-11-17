"""Synthetic R3 gate/collection tests; all process launches mocked, no MD."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import titan_path_benchmark_campaign as c
import titan_path_benchmark_quota as q
import titan_path_benchmark_runtime as engine
import titan_path_restart_policy as restart
import titan_path_benchmark_quota_smoke as helpers


class CampaignTests(unittest.TestCase):
    fake = helpers.QuotaTests.fake
    submissions = helpers.QuotaTests.submissions
    env = helpers.QuotaTests.env

    def setUp(self):
        helpers.QuotaTests.setUp(self)
        self.cells = []
        for block in range(5):
            branch = 'branch_a' if block % 2 == 0 else 'branch_b'
            for route in ('orange', 'blue'):
                for method in ('adaptive', 'tangent', 'fixed'):
                    self.cells.append({'cell': f'{route}__block_{block:02d}__{method}',
                        'block': block, 'route_label': route, 'method': method, 'source_branch': branch,
                        'rng_lock': {'source_checkpoint_sha256': restart.b1_policy()['branches'][branch]['equilibrium_checkpoint_sha256']}})
        q.publish(self.root / 'PLAN.json', {'cells': self.cells})
        q.publish(self.root / 'CAMPAIGN_POLICY.json', c.POLICY)
        q.publish(self.root / restart.POLICY_PATH, restart.b1_policy())
        q.publish(self.root / 'PROTOCOL.json', {
            'campaign_policy_sha256': q.digest(self.root / 'CAMPAIGN_POLICY.json'),
            'restart_policy_sha256': q.digest(self.root / restart.POLICY_PATH)})
        self.assertEqual(c.verify(self.root), self.cells)

    def complete(self, key, status='METHOD_STOP', rc=0):
        place = q.directory(self.root, key)
        job = q.read(place / 'ACCEPTED.json')['job_id']
        if key == 'build':
            q.publish(self.root / 'generated/PREPARE_READY.json', {'status': 'PASS_PREPARE_NO_MD', 'md_executed': False})
        else:
            cell = self.cells[int(key[3:])]
            count = 100 if status == 'BUDGET_EXHAUSTED' else 0
            q.publish(c.result_path(self.root, int(key[3:])), {'trajectory_id': cell['cell'],
                'method': cell['method'], 'route': cell['route_label'], 'block': cell['block'],
                'status': status, 'error': None, 'scientific_pass': False,
                'completed_segments': count, 'segments': [{}] * count})
        q.publish(place / 'DONE.json', {'stage': key, 'job_id': job, 'returncode': rc, 'error': None,
            'scientific_pass': False, 'completion_evidence': c.completion_evidence(self.root, key)})

    def prefix(self, count):
        q.advance(self.root); self.complete('build')
        for i in range(count):
            q.advance(self.root)
            self.complete(f'md_{i:02d}', status=sorted(c.NORMAL)[i % 4])

    def test_first_only_build(self):
        self.assertIn('SUBMITTED build', q.advance(self.root))
        self.assertEqual(len(self.submissions()), 1)

    def test_pilot_waits_for_all_twelve(self):
        self.prefix(11)
        self.assertIn('SUBMITTED md_11', q.advance(self.root))
        self.assertFalse((self.root / 'generated/PILOT_GATE.json').exists())
        self.active = {q.read(q.directory(self.root, 'md_11') / 'ACCEPTED.json')['job_id']}
        self.assertIn('WAIT_ACTIVE', q.advance(self.root))
        self.assertFalse(any('--array=12-12' in x for x in self.submissions()))

    def test_twelve_normal_terminal_outcomes_unlock_remaining(self):
        self.prefix(12)
        self.assertIn('SUBMITTED md_12', q.advance(self.root))
        pilot = q.read(self.root / 'generated/PILOT_GATE.json')
        self.assertEqual(len(pilot['evidence']), 13)
        self.assertFalse(pilot['physical_success_required'])
        self.assertFalse(pilot['scientific_pass'])

    def test_full_thirty_once_then_collection(self):
        self.prefix(30)
        self.assertIn('SUBMITTED collect', q.advance(self.root))
        self.assertEqual(len(self.submissions()), 32)
        for i in range(30):
            self.assertEqual(sum(f'--array={i}-{i}' in x for x in self.submissions()), 1)

    def test_pilot_failure_stops_md_and_latches(self):
        self.prefix(2); q.advance(self.root); self.complete('md_02', 'TECHNICAL_FAILURE', 1)
        self.assertIn('SUBMITTED collect', q.advance(self.root))
        self.assertTrue((self.root / 'generated/CAMPAIGN_HALT.json').exists())
        self.assertNotIn('--array=3-3', sum(self.submissions(), []))
        self.assertEqual(q.advance(self.root), 'HOLD_COLLECTION_INTERRUPTED_NO_RETRY')
        self.assertNotIn('--array=3-3', sum(self.submissions(), []))

    def test_production_failure_also_stops(self):
        self.prefix(15); q.advance(self.root); self.complete('md_15', 'TECHNICAL_FAILURE', 1)
        self.assertIn('SUBMITTED collect', q.advance(self.root))
        self.assertNotIn('--array=16-16', sum(self.submissions(), []))

    def test_worker_zero_rc_but_failure_status_rejected(self):
        self.prefix(0); q.advance(self.root); self.complete('md_00', 'TECHNICAL_FAILURE')
        self.assertIn('SUBMITTED collect', q.advance(self.root))

    def test_result_mutation_detected(self):
        self.prefix(1)
        path = c.result_path(self.root, 0)
        q.publish(path, {**q.read(path), 'status': 'TARGET_REACHED'}, replace=True)
        self.assertIn('SUBMITTED collect', q.advance(self.root))

    def test_missing_result_detected(self):
        self.prefix(1); c.result_path(self.root, 0).unlink()
        self.assertIn('SUBMITTED collect', q.advance(self.root))

    def test_corrupt_done_detected(self):
        self.prefix(1)
        (q.directory(self.root, 'md_00') / 'DONE.json').write_text('{')
        self.assertIn('SUBMITTED collect', q.advance(self.root))
        self.assertTrue((self.root / 'generated/CAMPAIGN_HALT.json').exists())
        self.assertEqual(self.submissions()[-1][-1], 'collect')

    def test_interrupted_md_collects_no_retry(self):
        self.prefix(0); q.advance(self.root)
        self.assertIn('SUBMITTED collect', q.advance(self.root))
        self.assertFalse(any('--array=1-1' in x for x in self.submissions()))

    def test_interrupted_build_collects_no_md(self):
        q.advance(self.root)
        self.assertIn('SUBMITTED collect', q.advance(self.root))
        self.assertFalse(any(x.startswith('--array=') for x in sum(self.submissions(), [])))

    def test_unknown_acceptance_holds_no_collector(self):
        self.prefix(0); self.behavior = 'timeout'; q.advance(self.root)
        self.behavior = 'ok'
        self.assertEqual(q.advance(self.root), 'HOLD_ACCEPTANCE_UNKNOWN')
        self.assertFalse((self.root / 'generated/CAMPAIGN_HALT.json').exists())

    def test_pilot_gate_mutation_halts(self):
        self.prefix(12); self.behavior = 'quota'; q.advance(self.root)
        q.publish(self.root / 'generated/PILOT_GATE.json', {'fake': True}, replace=True)
        self.behavior = 'ok'
        self.assertIn('SUBMITTED collect', q.advance(self.root))

    def test_worker_cannot_skip_pilot(self):
        self.prefix(0)
        with self.assertRaisesRegex(ValueError, 'OUT_OF_ORDER'):
            c.worker_guard(self.root, 'md_12')

    def test_policy_mutation_rejected_before_submit(self):
        q.publish(self.root / 'CAMPAIGN_POLICY.json', {**c.POLICY, 'pilot_count': 1}, replace=True)
        with self.assertRaisesRegex(ValueError, 'POLICY_CHANGED'): q.advance(self.root)
        self.assertFalse(self.submissions())

    def test_restart_policy_required(self):
        (self.root / restart.POLICY_PATH).unlink()
        with self.assertRaises(ValueError): q.advance(self.root)
        self.assertFalse(self.submissions())

    def test_plan_order_required(self):
        q.publish(self.root / 'PLAN.json', {'cells': self.cells[::-1]}, replace=True)
        with self.assertRaisesRegex(ValueError, 'ORDER'): q.advance(self.root)

    def test_real_wrapper_failure_halts(self):
        self.prefix(0); q.advance(self.root)
        with patch.dict(os.environ, self.env('md_00', '70002'), clear=True), patch.object(engine, 'verify'), patch.object(engine, 'run_trajectory', return_value=1):
            self.assertEqual(q.execute_stage(self.root), 1)
        self.assertEqual(self.submissions()[-1][-1], 'collect')

    def test_actual_collector_halted_denominator(self):
        self.prefix(2); q.advance(self.root); self.complete('md_02', 'TECHNICAL_FAILURE', 1)
        q.advance(self.root)
        report=self.collect_fixture()
        self.assertEqual(report['accounting']['planned_count'], 30)
        self.assertEqual(report['accounting']['technical_terminal_count'], 2)
        self.assertEqual(report['accounting']['technical_failure_count'], 1)
        self.assertEqual(report['accounting']['not_run_count'], 27)
        self.assertFalse(report['scientific_pass'])

    def collect_fixture(self):
        plan = {'cells': self.cells}
        with patch.object(engine, 'verify', return_value=plan), patch.object(engine, 'read', side_effect=lambda p: {'collect_job': '99'} if p.name == 'SUBMISSION.json' else q.read(p)), patch.dict(os.environ, {'SLURM_JOB_ID': '99'}):
            for name in ('RUNTIME_ROUTES.json', 'PACKAGE_LOCK.json'):
                q.publish(self.root / name, {'synthetic': True})
            engine.collect(self.root)
            engine.audit(self.root)
        return q.read(self.root / 'generated/COLLECTION.json')

    def test_corrupt_pilot_still_downloads_diagnostics(self):
        self.prefix(12); self.behavior='quota'; q.advance(self.root)
        (self.root/'generated/PILOT_GATE.json').write_text('{')
        self.behavior='ok'; q.advance(self.root)
        report=self.collect_fixture()
        self.assertEqual(report['status'],'FAILED_WITH_DIAGNOSTICS')
        self.assertIn('unreadable_receipt',report['pilot_gate'])
        self.assertEqual(report['accounting']['not_run_count'],18)

    def test_halt_cannot_be_reported_as_pass_with_thirty_valid_results(self):
        self.prefix(30)
        c.halt(self.root,'md_29','synthetic campaign-integrity fault')
        q.advance(self.root)
        report=self.collect_fixture()
        self.assertEqual(report['status'],'FAILED_WITH_DIAGNOSTICS')
        self.assertEqual(report['accounting']['technical_terminal_count'],30)
        self.assertTrue(report['campaign_halted'])

    def test_full_collection_no_scientific_pass(self):
        self.prefix(30); q.advance(self.root)
        reports = [q.read(c.result_path(self.root, i)) for i in range(30)]
        reports = c.classify_collection(self.root, self.cells, reports)
        accounting = c.accounting(self.cells, reports)
        self.assertEqual(accounting['technical_terminal_count'], 30)
        self.assertEqual(accounting['not_run_count'], 0)
        self.assertFalse(accounting['physical_success_evidence_validated_here'])


if __name__ == '__main__':
    unittest.main()
