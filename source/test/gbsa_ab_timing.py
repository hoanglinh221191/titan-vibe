"""A/B timing of the MM/GBSA kernel on the H10 fixture (ENE-15, ENE-04).

Runs two titan binaries alternately (order flipped every round) on the H10
window (chain A - chain B, 100 frames, igb 8, rgbmax 25, -nt 4), each run
pinned to the P-cores on Windows, and reads gbsa_performance.dat: the fused
energy kernel and the Born-radius/SASA precompute, in ms per frame. Prints
every round, the medians and the median of the per-round ratios B/A. This
machine drifts (hybrid CPU, clock): trust the per-round ratios, and run an
A/A pair (the same binary twice) to see the noise floor.

  python test/gbsa_ab_timing.py --a old/titan.exe --b titan.exe --work-root <dir outside OneDrive>
"""
import argparse
import ctypes
import os
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'test' / 'fixtures' / 'e2e_h10'
P_CORES = (1 << 0) | (1 << 1) | (1 << 10) | (1 << 11) | (1 << 12) | (1 << 13)


def args_for(fx, threads, frames_end):
    return ['-ff=charmm36', f'-s={fx}/protein_ref.pdb', f'-n={fx}/prot.ndx',
            f'-common_topo={fx}/topol.top', '-eu=kCal/mol',
            '-noprint-side_chain_contact=true', f'-nt={threads}',
            f'-traj={fx}/h10_window.xtc', '-g=1-2', '-b=8500', f'-e={frames_end}',
            '-dt=50', '-gbsa=true', '-igb=8', '-gbsa-rgbmax=25']


def run_once(binary, run_dir, arg_lines, affinity):
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    (run_dir / 'arg.in').write_text('\n'.join(arg_lines) + '\n', encoding='utf-8')
    with open(run_dir / 'stdout.log', 'wb') as so:
        proc = subprocess.Popen([str(binary), '-a', 'arg.in'], cwd=run_dir, stdout=so,
                                stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
        if os.name == 'nt' and affinity:
            ctypes.windll.kernel32.SetProcessAffinityMask(int(proc._handle), affinity)
        rc = proc.wait()
    perf = {}
    path = run_dir / 'gbsa_performance.dat'
    if rc != 0 or not path.exists():
        raise RuntimeError(f'{binary} failed in {run_dir} (rc {rc})')
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2 and not line.startswith('#'):
            try:
                perf[parts[0]] = float(parts[1])
            except ValueError:
                pass
    frames = perf['frames']
    return (perf['gbsa_energy_kernel_ms_total'] / frames,
            perf['gbsa_precompute_ms_total'] / frames)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--a', required=True, type=Path)
    p.add_argument('--b', required=True, type=Path)
    p.add_argument('--rounds', type=int, default=10)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--end', type=int, default=13450, help='-e in ps (100 frames at 13450)')
    p.add_argument('--work-root', required=True, type=Path)
    p.add_argument('--no-pin', action='store_true')
    a = p.parse_args()
    work = a.work_root.resolve()
    work.mkdir(parents=True, exist_ok=True)
    fx = work / '_fixture'
    if not fx.exists():
        shutil.copytree(FIXTURE, fx, ignore=shutil.ignore_patterns('*.py', '__pycache__',
                                                                    'h10_long.xtc'))
    bins = {'A': work / ('A_' + a.a.name), 'B': work / ('B_' + a.b.name)}
    shutil.copy2(a.a, bins['A'])
    shutil.copy2(a.b, bins['B'])
    lines = args_for(fx.as_posix(), a.threads, a.end)
    affinity = 0 if a.no_pin else P_CORES
    res = {'A': [], 'B': []}
    for r in range(a.rounds):
        order = ('A', 'B') if r % 2 == 0 else ('B', 'A')
        for tag in order:
            res[tag].append(run_once(bins[tag], work / f'run_{tag}', lines, affinity))
        ka, kb = res['A'][-1][0], res['B'][-1][0]
        print(f'round {r + 1:2d}  kernel A {ka:8.3f}  B {kb:8.3f} ms/frame  B/A {kb / ka:6.3f}   '
              f'precompute A {res["A"][-1][1]:7.3f}  B {res["B"][-1][1]:7.3f}', flush=True)
    med = {t: statistics.median(v[0] for v in res[t]) for t in res}
    pre = {t: statistics.median(v[1] for v in res[t]) for t in res}
    ratios = [b[0] / a_[0] for a_, b in zip(res['A'], res['B'])]
    print(f'median kernel ms/frame: A {med["A"]:.3f}  B {med["B"]:.3f}; '
          f'median round ratio B/A {statistics.median(ratios):.3f} '
          f'(min {min(ratios):.3f}, max {max(ratios):.3f}); '
          f'ratio of minima {min(v[0] for v in res["B"]) / min(v[0] for v in res["A"]):.3f}')
    print(f'median precompute ms/frame: A {pre["A"]:.3f}  B {pre["B"]:.3f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
