import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_receding_dispatch as d


class DispatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
    def tearDown(self): self.tmp.cleanup()
    def test_accepted_not_resubmitted(self):
        with patch.object(d.subprocess, 'run', return_value=SimpleNamespace(returncode=0, stdout='123;cluster\n', stderr='')) as run:
            self.assertEqual(d.submit_job(self.root, 'build', []), '123')
            self.assertEqual(d.submit_job(self.root, 'build', []), '123')
            self.assertEqual(run.call_count, 1)
    def test_ambiguous_never_blind_retry(self):
        with patch.object(d.subprocess, 'run', return_value=SimpleNamespace(returncode=0, stdout='unknown', stderr='')) as run:
            with self.assertRaises(ValueError): d.submit_job(self.root, 'build', [])
            with self.assertRaises(ValueError): d.submit_job(self.root, 'build', [])
            self.assertEqual(run.call_count, 1)
    def test_failed_submission_retained(self):
        with patch.object(d.subprocess, 'run', return_value=SimpleNamespace(returncode=1, stdout='', stderr='quota')):
            with self.assertRaises(ValueError): d.submit_job(self.root, 'build', [])
        self.assertEqual(d.read(self.root/'generated/jobs/build_RESPONSE.json')['stderr'], 'quota')
    def test_prepare_failure_only_collects(self):
        with patch.object(d, 'submit_job', return_value='4') as run:
            d.dispatch_after_prepare(self.root)
            self.assertEqual(run.call_args.args[1], 'collect')
    def test_one_md_element_only(self):
        (self.root/'generated/jobs').mkdir(parents=True)
        d.write(self.root/'generated/PREPARE_READY.json', {})
        d.write(self.root/'generated/jobs/build.json', {'job_id': '1'})
        d.write(self.root/'PLAN.json', {'cells': [{'cell': 'a'}, {'cell': 'b'}, {'cell': 'c'}]})
        with patch.object(d, 'submit_job', return_value='2') as run:
            d.dispatch_after_prepare(self.root)
            self.assertEqual(run.call_args.args[1:], ('md_00', ['--array=0']))
            self.assertEqual(run.call_count, 1)
    def test_technical_failure_collects_not_next_md(self):
        d.write(self.root/'generated/PREPARE_READY.json', {})
        d.write(self.root/'PLAN.json', {'cells': [{'cell': 'a'}, {'cell': 'b'}]})
        d.write(self.root/'outputs/a/TRAJECTORY_RESULT.json', {'status': 'TECHNICAL_FAILURE'})
        with patch.object(d, 'submit_job', return_value='2') as run:
            d.dispatch_after_prepare(self.root)
            self.assertEqual(run.call_args.args[1], 'collect')
    def test_timeout_blocks_retry(self):
        with patch.object(d.subprocess, 'run', side_effect=TimeoutError):
            with self.assertRaises(TimeoutError): d.submit_job(self.root, 'build', [])
        with self.assertRaises(ValueError): d.submit_job(self.root, 'build', [])


if __name__ == '__main__': unittest.main()
