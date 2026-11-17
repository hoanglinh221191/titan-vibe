"""No Slurm/MD: adversarial dispatch and per-job authorization checks."""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import titan_path_benchmark_quota as q


class QuotaTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        (self.root / 'FILES.sha256').write_text('immutable fixture')
        self.calls = []
        self.counter = 70000
        self.active = set()
        self.behavior = 'ok'
        self.runner = patch.object(q.subprocess, 'run', side_effect=self.fake)
        self.runner.start(); self.addCleanup(self.runner.stop)

    def fake(self, command, **kwargs):
        self.calls.append(command)
        if command[0] == 'squeue':
            return subprocess.CompletedProcess(command, 0, '\n'.join(sorted(self.active)), '')
        if self.behavior == 'quota':
            return subprocess.CompletedProcess(command, 1, '', 'sbatch: error: AssocMaxSubmitJobLimit')
        if self.behavior == 'timeout':
            raise subprocess.TimeoutExpired(command, 60)
        if self.behavior == 'ambiguous':
            return subprocess.CompletedProcess(command, 0, 'not a job ID', '')
        self.counter += 1
        return subprocess.CompletedProcess(command, 0, str(self.counter) + '\n', '')

    def submissions(self):
        return [c for c in self.calls if c[0] == 'sbatch']

    def done(self, key, ready=False):
        place = q.directory(self.root, key)
        job = q.read(place / 'ACCEPTED.json')['job_id']
        q.publish(place / 'DONE.json', {'job_id': job, 'returncode': 0})
        if ready:
            q.publish(self.root / 'generated/PREPARE_READY.json', {'synthetic': True})
        return job

    def env(self, key, job='70001'):
        req = q.request(self.root, key)
        intent = q.directory(self.root, key) / 'attempt_0000_INTENT.json'
        if not intent.exists():
            q.publish(intent, {'synthetic': True})
        value = {'B1_STAGE_KEY': key, 'SLURM_JOB_ID': job, 'SLURM_JOB_NAME': req['job_name'],
                 'SLURM_JOB_PARTITION': 'batch' if key.startswith('md_') else 'interact'}
        if key.startswith('md_'):
            value.update(SLURM_ARRAY_JOB_ID=job, SLURM_ARRAY_TASK_ID=str(req['index']))
        return value

    def test_first_submit_only_build(self):
        self.assertIn('SUBMITTED build', q.advance(self.root))
        self.assertEqual(len(self.submissions()), 1)
        self.assertFalse(any('--array=' in x for x in self.submissions()[0]))
        self.assertTrue((self.root / 'generated/SUBMISSION.json').is_file())

    def test_active_never_resubmitted(self):
        q.advance(self.root); self.active = {'70001'}
        self.assertIn('WAIT_ACTIVE', q.advance(self.root))
        self.assertEqual(len(self.submissions()), 1)

    def test_quota_rejection_explicit_retry(self):
        self.behavior = 'quota'
        self.assertIn('WAIT_QUOTA', q.advance(self.root))
        self.behavior = 'ok'; q.advance(self.root)
        place = q.directory(self.root, 'build')
        self.assertEqual(len(list(place.glob('*_INTENT.json'))), 2)
        self.assertEqual(q.read(place / 'attempt_0000_RESULT.json')['state'], 'REJECTED')

    def test_single_index_and_dependency(self):
        q.advance(self.root); job = self.done('build', ready=True)
        q.advance(self.root, predecessor=job)
        self.assertIn('--array=0-0', self.submissions()[-1])
        self.assertIn('--dependency=afterany:' + job, self.submissions()[-1])

    def test_all_thirty_untouched_indices_once(self):
        q.advance(self.root); self.done('build', ready=True)
        for i in range(30):
            q.advance(self.root)
            self.assertIn(f'--array={i}-{i}', self.submissions()[-1])
            self.done(f'md_{i:02d}')
        q.advance(self.root)
        self.assertEqual(len(self.submissions()), 32)
        self.assertTrue(self.submissions()[-1][-3].endswith('collect.slurm'))

    def test_failed_build_skips_md_and_collects(self):
        q.advance(self.root); self.done('build')
        self.assertIn('SUBMITTED collect', q.advance(self.root))
        self.assertFalse(any('--array=' in part for c in self.submissions() for part in c))

    def test_lost_build_not_restarted(self):
        q.advance(self.root)
        self.assertIn('SUBMITTED collect', q.advance(self.root))
        self.assertTrue((q.directory(self.root, 'build') / 'INTERRUPTED.json').exists())

    def test_lost_md_moves_only_to_next(self):
        q.advance(self.root); self.done('build', ready=True); q.advance(self.root)
        self.assertIn('SUBMITTED md_01', q.advance(self.root))
        self.assertTrue((q.directory(self.root, 'md_00') / 'INTERRUPTED.json').exists())

    def test_unknown_timeout_no_blind_retry(self):
        self.behavior = 'timeout'; q.advance(self.root)
        self.behavior = 'ok'
        self.assertEqual(q.advance(self.root), 'HOLD_ACCEPTANCE_UNKNOWN')
        self.assertEqual(len(self.submissions()), 1)

    def test_unknown_response_reconcile_live_job(self):
        self.behavior = 'ambiguous'; q.advance(self.root)
        self.active = {'80000'}
        self.assertEqual(q.advance(self.root), 'ACCEPTANCE_RECONCILED')
        self.assertEqual(len(self.submissions()), 1)
        self.assertIn('--name', self.calls[-1])

    def test_intent_without_response_holds(self):
        q.request(self.root, 'build')
        q.publish(q.directory(self.root, 'build') / 'attempt_0000_INTENT.json', {})
        self.assertEqual(q.advance(self.root), 'HOLD_ACCEPTANCE_UNKNOWN')
        self.assertFalse(self.submissions())

    def test_queue_unavailable_holds(self):
        q.advance(self.root)
        with patch.object(q, 'queue_ids', side_effect=ValueError('SQUEUE_UNAVAILABLE_HOLD')):
            with self.assertRaisesRegex(ValueError, 'SQUEUE_UNAVAILABLE'):
                q.advance(self.root)
        self.assertFalse((q.directory(self.root, 'build') / 'INTERRUPTED.json').exists())

    def test_queue_uses_base_array_id(self):
        q.queue_ids()
        self.assertIn('--format=%F', self.calls[-1])
        self.assertIn('--user', self.calls[-1])

    def test_worker_can_register_before_client_response(self):
        with patch.dict(os.environ, self.env('md_07'), clear=True):
            mapping = q.current_authorization(self.root)
        self.assertEqual(mapping, {'build_job': None, 'array_job': '70001', 'collect_job': None})
        self.assertTrue((q.directory(self.root, 'md_07') / 'ACCEPTED.json').is_file())

    def test_wrong_array_index_rejected(self):
        env = self.env('md_07'); env['SLURM_ARRAY_TASK_ID'] = '6'
        with patch.dict(os.environ, env, clear=True), self.assertRaisesRegex(ValueError, 'ARRAY_INDEX'):
            q.authorize(self.root)

    def test_wrong_job_name_rejected(self):
        env = self.env('build'); env['SLURM_JOB_NAME'] = 'other'
        with patch.dict(os.environ, env, clear=True), self.assertRaisesRegex(ValueError, 'JOB_NAME'):
            q.authorize(self.root)

    def test_wrong_partition_rejected(self):
        env = self.env('md_00'); env['SLURM_JOB_PARTITION'] = 'interact'
        with patch.dict(os.environ, env, clear=True), self.assertRaisesRegex(ValueError, 'PARTITION'):
            q.authorize(self.root)

    def test_accepted_job_id_conflict_rejected(self):
        env = self.env('build')
        with patch.dict(os.environ, env, clear=True):
            q.authorize(self.root)
            os.environ['SLURM_JOB_ID'] = '70002'
            with self.assertRaisesRegex(ValueError, 'CONFLICT'):
                q.authorize(self.root)

    def test_immutable_publish_does_not_overwrite(self):
        path = self.root / 'receipt.json'; q.publish(path, {'x': 1})
        q.publish(path, {'x': 1})
        with self.assertRaises(ValueError):
            q.publish(path, {'x': 2})
        self.assertEqual(q.read(path), {'x': 1})

    def test_static_request_mutation_rejected(self):
        q.request(self.root, 'build'); (self.root / 'FILES.sha256').write_text('mutated')
        with self.assertRaisesRegex(ValueError, 'BINDING_CHANGED'):
            q.request(self.root, 'build')

    def test_collection_interrupted_not_overwritten(self):
        q.advance(self.root); self.done('build'); q.advance(self.root)
        self.assertEqual(q.advance(self.root), 'HOLD_COLLECTION_INTERRUPTED_NO_RETRY')
        self.assertEqual(len(self.submissions()), 2)

    def test_ready_prevents_submission(self):
        q.publish(self.root / 'downloads/READY.json', {'synthetic': True})
        self.assertEqual(q.advance(self.root), 'DOWNLOAD_READY')
        self.assertFalse(self.calls)

    def test_status_before_submit_read_only(self):
        with contextlib.redirect_stdout(io.StringIO()):
            q.status(self.root)
        self.assertFalse((self.root / 'generated').exists())
        self.assertFalse(self.calls)

    def test_midchain_quota_no_repeat_md(self):
        q.advance(self.root); self.done('build', ready=True); q.advance(self.root); self.done('md_00')
        self.behavior = 'quota'; q.advance(self.root)
        self.behavior = 'ok'; q.advance(self.root)
        self.assertEqual(sum('--array=0-0' in c for c in self.submissions()), 1)
        self.assertEqual(sum('--array=1-1' in c for c in self.submissions()), 2)

    def test_submission_environment_drops_parent_allocation(self):
        seen = []
        def fake(command, **kwargs):
            seen.append(kwargs['env']); return subprocess.CompletedProcess(command, 0, '12345', '')
        with patch.dict(os.environ, {'SLURM_JOB_ID': '9', 'SBATCH_ARRAY_INX': '0-29', 'B1_STAGE_KEY': 'build'}), patch.object(q.subprocess, 'run', side_effect=fake):
            q.advance(self.root)
        self.assertFalse(any(k.startswith(('SLURM_', 'SBATCH_')) or k == 'B1_STAGE_KEY' for k in seen[0]))

    def test_real_stage_wrapper_failed_build_collects_without_md(self):
        import titan_path_benchmark_runtime as engine
        q.advance(self.root)
        with patch.dict(os.environ, self.env('build'), clear=True), patch.object(engine, 'verify'), patch.object(engine, 'prepare', side_effect=ValueError('SYNTHETIC_PREPARE_FAILURE')), patch.object(engine, 'run_trajectory') as md:
            self.assertEqual(q.execute_stage(self.root), 1)
        md.assert_not_called()
        self.assertEqual(self.submissions()[-1][-1], 'collect')
        self.assertFalse(q.read(self.root / 'generated/PREPARE_FAILURE.json')['md_executed'])

    def test_execution_claim_blocks_second_worker(self):
        import titan_path_benchmark_runtime as engine
        with patch.dict(os.environ, self.env('build'), clear=True), patch.object(engine, 'verify'), patch.object(engine, 'prepare') as prepare, patch.object(q, 'advance'):
            q.execute_stage(self.root)
            with self.assertRaises(FileExistsError):
                q.execute_stage(self.root)
        prepare.assert_called_once()

    def test_md_technical_failure_advances_only_next_stage(self):
        import titan_path_benchmark_runtime as engine
        q.advance(self.root); self.done('build', ready=True); q.advance(self.root)
        with patch.dict(os.environ, self.env('md_00', '70002'), clear=True), patch.object(engine, 'verify'), patch.object(engine, 'run_trajectory', return_value=1) as md:
            self.assertEqual(q.execute_stage(self.root), 1)
        md.assert_called_once_with(self.root, 0)
        self.assertIn('--array=1-1', self.submissions()[-1])
        self.assertEqual(q.read(q.directory(self.root, 'md_00') / 'DONE.json')['returncode'], 1)

    def test_collection_exception_preserved_no_automatic_retry(self):
        import titan_path_benchmark_runtime as engine
        with patch.dict(os.environ, self.env('collect'), clear=True), patch.object(engine, 'verify'), patch.object(engine, 'collect', side_effect=ValueError('SYNTHETIC_COLLECT_FAILURE')):
            self.assertEqual(q.execute_stage(self.root), 1)
        self.assertFalse(self.submissions())
        self.assertEqual(q.read(q.directory(self.root, 'collect') / 'DONE.json')['returncode'], 1)

    def test_shipped_authorization_hooks_match_role(self):
        import titan_path_benchmark_runtime as engine
        import titan_path_paired_loop as paired
        # Enabled only in the R2 shipped modules; source R1 remains unchanged.
        if not hasattr(engine, 'read') or engine.read is engine.replay.read:
            self.skipTest('R2 package hook integration is tested against shipped R2')
        q.publish(self.root / 'SCHEDULER_POLICY.json', {'synthetic': True})
        with patch.dict(os.environ, self.env('md_03'), clear=True):
            for module in (engine, paired):
                mapping = module.read(self.root / 'generated/SUBMISSION.json')
                self.assertEqual(mapping, {'build_job': None, 'array_job': '70001', 'collect_job': None})

    def test_broken_handoff_retains_done_without_retry(self):
        import titan_path_benchmark_runtime as engine
        with patch.dict(os.environ, self.env('md_00'), clear=True), patch.object(engine, 'verify'), patch.object(engine, 'run_trajectory', return_value=0), patch.object(q, 'advance', side_effect=ValueError('QUEUE_DOWN')):
            self.assertEqual(q.execute_stage(self.root), 0)
        place = q.directory(self.root, 'md_00')
        self.assertTrue((place / 'DONE.json').exists())
        self.assertTrue((place / 'HANDOFF_PAUSED.json').exists())


if __name__ == '__main__':
    unittest.main()
