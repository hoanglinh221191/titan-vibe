"""Actual archive collection of synthetic stdout/stderr; no scheduler or MD."""
import hashlib
import json
import os
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_benchmark_runtime as runtime
import titan_path_paired_loop as paired


class DiagnosticsTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.cell={'cell':'synthetic','method':'fixed'}
        for name,value in {'PLAN.json':{},'PROTOCOL.json':{},'RUNTIME_ROUTES.json':{},'PACKAGE_LOCK.json':{},
                           'generated/SUBMISSION.json':{'collect_job':'12'}}.items():
            p=self.root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(value))
        (self.root/'FILES.sha256').write_text('synthetic\n')
        self.selection=self.root/'outputs/trajectories/synthetic/segment_081/selection'
        self.selection.mkdir(parents=True)

    def collect(self):
        with patch.object(runtime,'verify',return_value={'cells':[self.cell]}),patch.dict(os.environ,{'SLURM_JOB_ID':'12'}):
            runtime.collect(self.root);runtime.audit(self.root)
        with tarfile.open(self.root/'downloads'/runtime.BUNDLE) as archive:
            return {m.name:archive.extractfile(m).read() for m in archive}

    def test_native_stdout_and_stderr_retained(self):
        (self.selection/'native.log').write_bytes(b'partial output\n')
        (self.selection/'native.log.stderr').write_bytes(b'[ERROR] route window\n')
        data=self.collect()
        for n in ('native.log','native.log.stderr'):
            name=(self.selection/n).relative_to(self.root).as_posix()
            self.assertEqual(data[name],(self.selection/n).read_bytes())
            manifest=json.loads(data['ANALYSIS_MANIFEST.json'])
            self.assertEqual(manifest[name]['sha256'],hashlib.sha256(data[name]).hexdigest())

    def test_empty_stderr_preserved(self):
        (self.selection/'native.log.stderr').write_bytes(b'')
        data=self.collect()
        self.assertEqual(data[(self.selection/'native.log.stderr').relative_to(self.root).as_posix()],b'')

    def test_non_utf8_stderr_preserved_verbatim(self):
        (self.selection/'native.log.stderr').write_bytes(b'\xff\x00error')
        data=self.collect()
        self.assertEqual(data[(self.selection/'native.log.stderr').relative_to(self.root).as_posix()],b'\xff\x00error')

    def test_failure_status_not_changed_by_collection(self):
        data=self.collect();report=json.loads(data['generated/COLLECTION.json'])
        self.assertEqual(report['status'],'FAILED_WITH_DIAGNOSTICS')
        self.assertFalse(report['scientific_pass'])
        self.assertEqual(report['trajectories'][0]['status'],'TECHNICAL_FAILURE')

    def test_veto_and_other_stderr_but_not_large_md_log(self):
        for name in ('veto.log','veto.log.stderr','grompp.log.stderr','run.part0001.log'):
            (self.selection/name).write_text('synthetic')
        data=self.collect()
        for name in ('veto.log','veto.log.stderr','grompp.log.stderr'):
            self.assertIn((self.selection/name).relative_to(self.root).as_posix(),data)
        self.assertNotIn((self.selection/'run.part0001.log').relative_to(self.root).as_posix(),data)

    def test_runner_captures_error_channel_before_raising(self):
        log=self.selection/'native.log'
        with self.assertRaisesRegex(ValueError,'COMMAND_FAILED'):
            paired.logged([sys.executable,'-c','import sys;print("stdout");print("stderr",file=sys.stderr);sys.exit(1)'],self.selection,log)
        self.assertIn('stdout',log.read_text())
        self.assertIn('stderr',log.with_suffix('.log.stderr').read_text())
        data=self.collect()
        self.assertIn(log.with_suffix('.log.stderr').relative_to(self.root).as_posix(),data)


if __name__=='__main__':unittest.main()
