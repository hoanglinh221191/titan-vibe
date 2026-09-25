"""Independent synthetic scheduler review: never invokes a native scheduler/MD.

Copies the real runner/gate into disposable fixtures, replaces only the cluster
Python environment loader, and injects shell functions for sbatch/scontrol.
"""
from pathlib import Path
import ast
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / 'campaigns/titan_path_b1_nav5_short_smd_20260905_r1'
BASH = Path('C:/Program Files/Git/bin/bash.exe')
COMMIT = 'f2aedcc934d48559ef1979428fb7f4871253841b'
MOCK = r'''
python3() {
    if [[ "${3:-}" == record && "${MOCK_FAIL:-}" == record ]]; then return 17; fi
    "$MOCK_PYTHON" "$@"
}
sbatch() {
    local n=0
    [[ ! -f "$MOCK_STATE/count" ]] || n="$(<"$MOCK_STATE/count")"
    n=$((n+1))
    printf '%s' "$n" > "$MOCK_STATE/count"
    printf 'SBATCH %s\n' "$*" >> "$MOCK_STATE/events"
    if [[ "${MOCK_FAIL:-}" == "$n" ]]; then return 19; fi
    printf '%s\n' "$((100+n))"
}
scontrol() {
    [[ "$1" == release && -s "$MOCK_PACKAGE/generated/SUBMISSION.json" ]] || return 21
    printf 'RELEASE %s\n' "$*" >> "$MOCK_STATE/events"
    [[ "${MOCK_FAIL:-}" != release ]] || return 23
}
export -f python3 sbatch scontrol
'''


class SubmissionReview(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='nav5-mock-review-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.root = self.base / 'package with spaces'
        (self.root / 'scripts').mkdir(parents=True)
        (self.root / 'jobs').mkdir()
        shutil.copy2(SOURCE / 'run_nav5.sh', self.root / 'run_nav5.sh')
        shutil.copy2(SOURCE / 'scripts/package_gate.py', self.root / 'scripts/package_gate.py')
        for path in (SOURCE / 'jobs').glob('*.slurm'):
            shutil.copy2(path, self.root / 'jobs' / path.name)
        (self.root / 'scripts/python_env.sh').write_text(
            '# SYNTHETIC: cluster Python/module setup is not exercised here\n', encoding='utf-8')
        self.state = self.base / 'mock_state'
        self.state.mkdir()
        self.pilot = self.base / 'pilot with spaces'
        self.pilot.mkdir()
        self.mock = self.base / 'mock.sh'
        self.mock.write_text(MOCK, encoding='utf-8', newline='\n')
        self.env = {**os.environ, 'BASH_ENV': self.mock.as_posix(),
                    'MOCK_PYTHON': Path(sys.executable).as_posix(),
                    'MOCK_STATE': self.state.as_posix(), 'MOCK_PACKAGE': self.root.as_posix(),
                    'PYTHONDONTWRITEBYTECODE': '1'}
        self.gate('seal', '--source-commit', COMMIT, wanted=0)

    def gate(self, *args, wanted=None, env=None):
        result = subprocess.run([sys.executable, '-B', str(self.root / 'scripts/package_gate.py'), *args],
                                cwd=self.root, env=env or self.env, capture_output=True, text=True)
        if wanted is not None:
            self.assertEqual(result.returncode, wanted, result.stdout + result.stderr)
        return result

    def submit(self, failure=''):
        return subprocess.run([str(BASH), str(self.root / 'run_nav5.sh'), 'submit', self.pilot.as_posix()],
                              env={**self.env, 'MOCK_FAIL': failure}, capture_output=True, text=True)

    def events(self):
        path = self.state / 'events'
        return path.read_text() if path.exists() else ''

    def test_success_records_before_release_and_preserves_dependencies(self):
        result = self.submit()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        events = self.events().splitlines()
        self.assertEqual(len(events), 4)
        self.assertIn('--hold', events[0])
        self.assertIn('--dependency=afterany:101', events[1])
        self.assertIn('--array=0-4%4', events[1])
        self.assertIn('--dependency=afterany:102', events[2])
        self.assertEqual(events[3], 'RELEASE release 101')
        record = json.loads((self.root / 'generated/SUBMISSION.json').read_text())
        self.assertEqual(record['jobs'], {'prepare': '101', 'array': '102', 'collect': '103'})
        self.assertEqual(Path(record['pilot_root']).resolve(), self.pilot.resolve())
        second = self.submit()
        self.assertNotEqual(second.returncode, 0)
        self.assertEqual(self.events().splitlines(), events)

    def test_scheduler_failure_never_releases(self):
        result = self.submit('2')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn('RELEASE ', self.events())
        self.assertFalse((self.root / 'generated/SUBMISSION.json').exists())
        self.assertIn('SUBMISSION_INCOMPLETE', result.stderr)

    def test_record_failure_never_releases(self):
        result = self.submit('record')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn('RELEASE ', self.events())
        self.assertFalse((self.root / 'generated/SUBMISSION.json').exists())

    def test_release_failure_is_explicitly_uncertain(self):
        result = self.submit('release')
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue((self.root / 'generated/SUBMISSION.json').is_file())
        self.assertIn('SUBMISSION_RELEASE_UNCERTAIN', result.stderr)

    def test_job_role_and_native_id_required(self):
        self.assertEqual(self.submit().returncode, 0)
        clean = {key: value for key, value in self.env.items()
                 if key not in ('SLURM_JOB_ID', 'SLURM_ARRAY_JOB_ID')}
        self.assertNotEqual(self.gate('job-check', env=clean).returncode, 0)
        self.assertNotEqual(self.gate('job-check', '--role', 'prepare', env=clean).returncode, 0)
        self.gate('job-check', '--role', 'prepare', env={**clean, 'SLURM_JOB_ID': '101'}, wanted=0)
        self.assertNotEqual(self.gate('job-check', '--role', 'prepare', env={**clean, 'SLURM_JOB_ID': '999'}).returncode, 0)

    def test_static_drift_and_cache_rejected(self):
        cache = self.root / 'scripts/__pycache__'
        cache.mkdir()
        (cache / 'unexpected.pyc').write_bytes(b'SYNTHETIC cached data')
        self.assertNotEqual(self.gate('verify').returncode, 0)

    def test_prepared_manifest_cannot_be_replaced_under_old_marker(self):
        def sha(path):
            return hashlib.sha256(path.read_bytes()).hexdigest()
        generated = self.root / 'generated'
        generated.mkdir()
        payload = generated / 'payload.json'
        payload.write_text('{"kind":"SYNTHETIC"}\n')
        report = generated / 'NAV5_PREPARATION_REPORT.json'
        report.write_text('{"kind":"SYNTHETIC report"}\n')
        manifest = generated / 'nav5_preparation_artifacts.sha256'
        original = f'{sha(payload)}  generated/payload.json\n'
        manifest.write_text(original)
        marker = generated / 'ALL5_PREPARED.json'
        marker.write_text(json.dumps({
            'cells': ['branch_a__blue','branch_a__cyan','branch_b__cyan','branch_b__orange','branch_b__pink'],
            'status': 'PASS_ALL5_NATIVE_PREFLIGHT_NO_MD', 'md_executed': False,
            'preparation_artifacts_sha256': sha(manifest), 'report_sha256': sha(report)}))
        marker.with_suffix('.json.sha256').write_text(f'{sha(marker)}  generated/ALL5_PREPARED.json\n')
        self.gate('prepared-check', wanted=0)
        payload.write_text('{"kind":"SYNTHETIC changed"}\n')
        manifest.write_text(f'{sha(payload)}  generated/payload.json\n')
        result = self.gate('prepared-check')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('not bound', result.stderr)
        manifest.write_text(original)
        result = self.gate('prepared-check')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('artifact identity', result.stderr)

    def bootstrap(self):
        tree = ast.parse((REPO / 'tools/package_titan_path_nav5.py').read_text(encoding='utf-8'))
        return next(ast.literal_eval(node.value) for node in tree.body
                    if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name)
                    and target.id == 'BOOTSTRAP' for target in node.targets))

    def test_bootstrap_syntax_and_argument_forwarding(self):
        script = self.bootstrap()
        result = subprocess.run([str(BASH), '-n'], input=script, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        forward = next(line for line in script.splitlines() if line.startswith('exec bash '))
        mocked = 'bash() { printf "<%s>\\n" "$@"; }; NAME=synthetic; ' + forward.replace('exec bash ', 'bash ', 1)
        for args, expected in (([], ['<synthetic/run_nav5.sh>', '<verify>']),
                               (['submit', 'pilot with spaces'], ['<synthetic/run_nav5.sh>', '<submit>', '<pilot with spaces>'])):
            result = subprocess.run([str(BASH), '-c', mocked, 'mock', *args], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.splitlines(), expected)

    def test_bootstrap_safety_check_rejects_single_backslash(self):
        script = self.bootstrap()
        python = script.split("<<'PY'\n", 1)[1].split('\nPY\n', 1)[0]
        # Only the host-Python-version test is replaced; archive safety logic
        # is exactly the generated bootstrap snippet.
        python = python.replace("assert sys.version_info[:2] == (3, 12), 'Python 3.12 required'", 'pass')
        name = 'titan_path_b1_nav5_short_smd_20260905_r1'
        for index, relative in enumerate(('safe.txt', 'one' + chr(92) + 'bad.txt')):
            archive = self.base / f'bootstrap_{index}.tar.gz'
            with tarfile.open(archive, 'w:gz') as handle:
                entry = tarfile.TarInfo(name + '/' + relative)
                entry.size = 1
                handle.addfile(entry, io.BytesIO(b'x'))
            result = subprocess.run([sys.executable, '-c', python, str(archive), name], capture_output=True, text=True)
            self.assertEqual(result.returncode == 0, index == 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
