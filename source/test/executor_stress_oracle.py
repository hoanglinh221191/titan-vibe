"""CORE-08 gate: TitanExecutor with two threads submitting batches concurrently.

Builds test/titan_executor_stress_oracle.cpp with src/titan_executor.cpp and
runs it under a timeout (the old executor hung: lost wake-up on a shared
condition variable). --executor-source builds against another
titan_executor.cpp (for example the old one, to see it fail). Keep each run
in a fresh directory.

  python test/executor_stress_oracle.py --work-root test/tmp/executor_stress
"""
import argparse
import os
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--compiler', default='clang++' if os.name == 'nt' else 'g++-13')
p.add_argument('--work-root', type=Path, required=True)
p.add_argument('--executor-source', type=Path, default=root / 'src/titan_executor.cpp')
p.add_argument('--rounds', type=int, default=4000)
p.add_argument('--stall-seconds', type=int, default=20)
p.add_argument('--sanitize', choices=['none', 'thread', 'address'], default='none')
a = p.parse_args()
work = a.work_root.resolve()
work.mkdir(parents=True, exist_ok=False)
exe = work / ('executor_stress' + ('.exe' if os.name == 'nt' else ''))
flags = ['-std=c++23', '-O2', '-pthread', '-Isrc']
if a.sanitize != 'none':
    flags = ['-std=c++23', '-O1', '-g', '-pthread', '-Isrc', f'-fsanitize={a.sanitize}']
build = subprocess.run([a.compiler, *flags, 'test/titan_executor_stress_oracle.cpp',
                        str(a.executor_source.resolve()), '-o', str(exe)],
                       cwd=root, capture_output=True, text=True, encoding='utf-8', errors='replace')
(work / 'build.log').write_text(build.stdout + build.stderr, encoding='utf-8')
if build.returncode:
    print(build.stdout + build.stderr, end='')
    raise SystemExit(build.returncode)
try:
    run = subprocess.run([str(exe), str(a.rounds), str(a.stall_seconds)], cwd=work,
                         capture_output=True, text=True, encoding='utf-8', errors='replace',
                         timeout=10 * a.stall_seconds + 120)
    output, code = run.stdout + run.stderr, run.returncode
except subprocess.TimeoutExpired as e:
    output = (e.stdout or '') if isinstance(e.stdout, str) else ''
    output += '\nFAIL executor stress: killed after the outer timeout\n'
    code = 3
(work / 'run.log').write_text(output, encoding='utf-8')
print(output, end='')
raise SystemExit(code)
