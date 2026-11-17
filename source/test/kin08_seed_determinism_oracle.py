"""KIN-08 / X-06 gate: the randomized and subspace eigensolvers are seeded.

On the H10 e2e fixture (370 frames over 738 ns):
  msm   exact TICA, 12-state MSM with the partial ITS scan (-msm-its-solver
        partial, 12 lags in parallel) and partial PCA (-pca-solver partial)
  rand  randomized TICA (-tica-backend randomized, rank 6)
Checks:
  - two runs of each with -nt 4 write identical files (timing and memory
    reports excepted);
  - msm with -nt 1 and -nt 4 write identical ITS/PCA/MSM files (each solve
    builds its own generator, so the thread layout does not matter);
  - -msm-seed 7 changes the randomized TICA eigenvalues (the seed reaches
    the sketch; at rank 6 the approximation depends on it visibly).
The old start blocks came from the process-wide, unseeded rand() behind
Eigen::MatrixXd::Random, the TICA sketch from a hard-coded 12345.

  python -I test/kin08_seed_determinism_oracle.py --binary titan.exe --work-root test/tmp/kin08_seed
"""
import argparse
import filecmp
import shutil
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
fx = (root / 'test' / 'fixtures' / 'e2e_h10').as_posix()
BASE = ['-ff=charmm36', f'-s={fx}/protein_ref.pdb', f'-n={fx}/prot.ndx',
        f'-common_topo={fx}/topol.top', '-eu=kCal/mol', '-noprint-side_chain_contact=true',
        f'-traj={fx}/h10_long.xtc', '-g=0-0', '-nocal=true',
        '-tica', '-tica-lag=5', '-tica-components=2']
MODES = {
    'msm': ['-tica-rank=8', '-msm', '-msm-lag=1', '-msm-clusters=12', '-msm-macrostates=2',
            '-msm-its-solver=partial', '-msm-its-topk=4',
            '-msm-its-lags=1,2,3,4,5,6,7,8,9,10,11,12',
            '-pca', '-pca-solver=partial', '-pca-topk=8'],
    'rand': ['-tica-backend=randomized', '-tica-rank=6'],
}
SKIP = {'arg.in', 'run.log', 'log.txt', 'cal.txt', 'memory_report.dat', 'performance_report.json'}


def run(binary, work, name, mode, nt, extra=()):
    case = work / name
    case.mkdir(parents=True)
    (case / 'arg.in').write_text('\n'.join(BASE + MODES[mode] + [f'-nt={nt}', *extra]) + '\n',
                                 encoding='ascii')
    r = subprocess.run([str(binary), '-a', 'arg.in'], cwd=case, capture_output=True, text=True,
                       encoding='utf-8', errors='replace', timeout=1200)
    (case / 'run.log').write_text(r.stdout + r.stderr, encoding='utf-8')
    if r.returncode != 0:
        raise SystemExit(f'FAIL {name}: exit {r.returncode}\n{(r.stdout + r.stderr)[-800:]}')
    return case


def compare(a, b, label, only=None):
    names = sorted(p.name for p in a.iterdir() if p.is_file() and p.name not in SKIP)
    if only:
        names = [n for n in names if n.startswith(only)]
    differ = [n for n in names if not (b / n).exists() or not filecmp.cmp(a / n, b / n, shallow=False)]
    ok = bool(names) and not differ
    print(f"[{'PASS' if ok else 'FAIL'}] {label}: {len(names) - len(differ)}/{len(names)} files identical"
          + (f'; differ: {differ}' if differ else ''))
    return ok


def eigenvalues(case):
    rows = [l.split() for l in (case / 'tica_eigenvalues.dat').read_text().splitlines()
            if l.strip() and not l.startswith('#')]
    return [float(r[1]) for r in rows]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--work-root', type=Path, required=True)
    a = p.parse_args()
    binary = a.binary.resolve()
    work = a.work_root.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    ok = True
    m1 = run(binary, work, 'msm_nt4_a', 'msm', 4)
    m2 = run(binary, work, 'msm_nt4_b', 'msm', 4)
    m3 = run(binary, work, 'msm_nt1', 'msm', 1)
    ok &= compare(m1, m2, 'partial ITS/PCA, two runs')
    ok &= compare(m1, m3, 'partial ITS/PCA, -nt 4 vs -nt 1', only=('msm_', 'pca_'))
    usage = [l for l in (m1 / 'log.txt').read_text(errors='replace').splitlines()
             if 'ITS Solver Usage' in l or 'PCA Solver' in l]
    used = any('partial(subspace)' in l for l in usage if 'ITS' in l) and \
        any('partial(subspace)' in l for l in usage if 'PCA' in l)
    print(f"[{'PASS' if used else 'FAIL'}] partial solvers were used: {usage}")
    ok &= used
    r1 = run(binary, work, 'rand_a', 'rand', 4)
    r2 = run(binary, work, 'rand_b', 'rand', 4)
    r7 = run(binary, work, 'rand_seed7', 'rand', 4, ['-msm-seed=7'])
    ok &= compare(r1, r2, 'randomized TICA, two runs')
    e1, e7 = eigenvalues(r1), eigenvalues(r7)
    changed = len(e1) == len(e7) and max(abs(x - y) for x, y in zip(e1, e7)) > 1e-6
    print(f"[{'PASS' if changed else 'FAIL'}] -msm-seed reaches the TICA sketch: seed 42 {e1}, seed 7 {e7}")
    ok &= changed
    print('RESULT:', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
