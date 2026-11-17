"""GRID-03 gate: the Galerkin multigrid of the grid PB solvers.

Builds test/titan_galerkin_mg_oracle.cpp with src/pbe/grid/galerkin_mg.cpp
and src/titan_executor.cpp, and runs it: coarse operators against P^T A P,
CG against a dense solve, convergence on a dielectric sphere, thread-count
invariance. Keep each run in a fresh directory.

  python test/grid_galerkin_mg_oracle.py --work-root test/tmp/grid_galerkin_mg
"""
import argparse
import os
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--compiler', default='clang++' if os.name == 'nt' else 'g++-13')
p.add_argument('--work-root', type=Path, required=True)
p.add_argument('--sanitize', action='store_true')
a = p.parse_args()
work = a.work_root.resolve()
work.mkdir(parents=True, exist_ok=False)
exe = work / ('galerkin_mg' + ('.exe' if os.name == 'nt' else ''))
flags = ['-std=c++23', '-O1' if a.sanitize else '-O2', '-Isrc', '-Isrc/eigen', '-pthread']
if a.sanitize:
    flags += ['-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer']
build = subprocess.run([a.compiler] + flags + ['test/titan_galerkin_mg_oracle.cpp',
                        'src/pbe/grid/galerkin_mg.cpp', 'src/titan_executor.cpp', '-o', str(exe)],
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
