"""KIN-09 gate: symmetric spectrum of reversible MSM transition matrices.

Builds test/titan_reversible_spectrum_oracle.cpp against the header
src/titan_markov_estimator.h and runs it. Keep each run in a fresh directory.

  python test/kin_spectrum_oracle.py --work-root test/tmp/kin_spectrum
"""
import argparse
import os
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--compiler', default='clang++' if os.name == 'nt' else 'g++-13')
p.add_argument('--work-root', type=Path, required=True)
a = p.parse_args()
work = a.work_root.resolve()
work.mkdir(parents=True, exist_ok=False)
exe = work / ('reversible_spectrum' + ('.exe' if os.name == 'nt' else ''))
build = subprocess.run([a.compiler, '-std=c++23', '-O2', '-Isrc', '-Isrc/eigen',
                        'test/titan_reversible_spectrum_oracle.cpp', '-o', str(exe)],
                       cwd=root, capture_output=True, text=True, encoding='utf-8', errors='replace')
(work / 'build.log').write_text(build.stdout + build.stderr, encoding='utf-8')
if build.returncode:
    print(build.stdout + build.stderr, end='')
    raise SystemExit(build.returncode)
run = subprocess.run([str(exe)], cwd=work, capture_output=True, text=True,
                     encoding='utf-8', errors='replace')
(work / 'run.log').write_text(run.stdout + run.stderr, encoding='utf-8')
print(run.stdout + run.stderr, end='')
raise SystemExit(run.returncode)
