"""Oracle for the MSD fit window and memory bound (STR-13, MDA-01/04).

Fifty independent atoms do free Brownian motion, D = 0.05 A^2/ps, dt = 1 ps,
500 frames (fixed seed), in a 300 A box; -msd-mode=atom, automatic maxlag
(min(frames - 1, 1000) = 499) and automatic fit window.

window  The automatic fit window must leave out lags with fewer than half of
        the frames as time origins: 20-80 % of the lags with lag <= frames/2
        (and at least the first three lags), here 50-201 ps. The old window,
        20-80 % of all 499 lags (100-400 ps), reached lags with 100 origins
        out of 500. D must equal the slope/6 recomputed here from the same
        coordinates over that window (1e-6 relative) and lie within 15 % of
        the true D (statistical check, fixed seed; ~7 % sampling error).
memory  -msd-memory-limit=0.05 (MB): storing every frame of 50 atoms needs
        more after ~40 frames; the run must stop with a message naming the
        limit instead of growing without bound (the old code kept every
        coordinate of every frame, 24 B per atom and frame, with no check).

Runs the energy path (-g=1-2, a fixed one-residue chain B as molecule 2),
as msd_time_block_oracle.py: -nocal re-centres molecule 1 every frame.

  python -I test/msd_window_memory_oracle.py --binary titan.exe --work-root D:/softs/titan_work/str/msd_window_oracle
"""
import argparse
import math
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 300.0
D_TRUE = 0.05
FRAMES = 500
ATOMS = 50
DT = 1.0


def walks(seed=7):
    rng = random.Random(seed)
    s = math.sqrt(2.0 * D_TRUE * DT)
    pos = [[150.0 + 3.0 * (a % 5), 150.0 + 3.0 * (a // 5), 150.0] for a in range(ATOMS)]
    frames = []
    for _ in range(FRAMES):
        frames.append([tuple(round(v, 3) for v in p) for p in pos])
        for p in pos:
            for k in range(3):
                p[k] += rng.gauss(0.0, s)
    return frames


def msd_curve(frames):
    """Atom-mode MSD for lags 1..frames-1 (origins t, targets t + lag)."""
    n = len(frames)
    out = []
    for lag in range(1, n):
        acc = 0.0
        for o in range(n - lag):
            a, b = frames[o], frames[o + lag]
            acc += sum((b[i][k] - a[i][k]) ** 2 for i in range(ATOMS) for k in range(3))
        out.append(acc / ((n - lag) * ATOMS))
    return out


def window(n_frames, lags, enforce_origins):
    eligible = [lag for lag in lags if not enforce_origins or 2 * lag <= n_frames]
    if len(eligible) < 3:
        eligible = lags[:3]
    m = len(eligible)
    return eligible[math.floor(0.2 * (m - 1))], eligible[math.ceil(0.8 * (m - 1))]


def fit_d(curve, lo, hi):
    xs = [lag * DT for lag in range(lo, hi + 1)]
    ys = [curve[lag - 1] for lag in range(lo, hi + 1)]
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx) / 6.0


def run_case(binary, work, name, frames, extra):
    case = work / name
    case.mkdir()
    traj = []
    for i, coords in enumerate(frames):
        atoms = [("A", a + 1, "MOL", "C1", *coords[a]) for a in range(ATOMS)]
        atoms += [("B", 1, "GLY", n, x, 20.0, 20.0)
                  for n, x in (("N", 20.0), ("CA", 21.4), ("C", 22.6), ("O", 23.8))]
        traj.append((DT * i, (BOX, BOX, BOX), atoms))
    fx.write_trajectory(case / "walk.pdb", traj)
    fx.write_topology(case / "walk.top", traj[0][2])
    fx.write_index(case / "walk.ndx", [("System", range(1, ATOMS + 5)),
                                       ("walkers", range(1, ATOMS + 1)),
                                       ("fixed", range(ATOMS + 1, ATOMS + 5))])
    lines = [x for x in fx.base_args(work, case / "walk.pdb", case / "walk.top",
                                     case / "walk.ndx", groups="1-2") if x != "-nocal=true"]
    lines += ["-msd=true", "-msd-mode=atom", "-msd-selection=mol1"] + extra
    return fx.run_titan(binary, case / "run", lines, timeout=1800) + (case / "run",)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    a = p.parse_args()
    binary = a.binary.resolve()
    work = a.work_root.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    checks = fx.Checks()

    frames = walks()
    curve = msd_curve(frames)
    lags = list(range(1, FRAMES))
    lo, hi = window(FRAMES, lags, True)
    old_lo, old_hi = window(FRAMES, lags, False)
    d_new, d_old = fit_d(curve, lo, hi), fit_d(curve, old_lo, old_hi)
    print(f"window {lo}-{hi} ps: D {d_new:.6f} ({100 * (d_new / D_TRUE - 1):+.1f} %); "
          f"old window {old_lo}-{old_hi} ps: D {d_old:.6f} ({100 * (d_old / D_TRUE - 1):+.1f} %)")

    code, output, run = run_case(binary, work, "window", frames, [])
    if checks.check(code == 0, "window: titan exit 0", output[-800:] if code else ""):
        files = [f for f in run.glob("msd_diffusion_summary*.dat") if "average" not in f.name]
        kv = fx.read_key_values(files[0]) if files else {}
        got = [float(v) for v in kv.get("fit_window_ps", "nan,nan").split(",")]
        checks.check(got == [lo * DT, hi * DT],
                     f"fit window {lo}-{hi} ps: no fitted lag has fewer than half the frames as origins",
                     f"titan {got}")
        d = float(kv.get("D_A2_per_ps", "nan"))
        checks.check(abs(d - d_new) <= 1e-6 * d_new, "D = slope/6 over that window, recomputed here",
                     f"titan {d}, here {d_new:.8f}")
        checks.check(abs(d / D_TRUE - 1.0) < 0.15, "D within 15 % of the true D (fixed seed)",
                     f"titan {d}, true {D_TRUE}")
    code, output, run = run_case(binary, work, "memory", frames, ["-msd-memory-limit=0.05"])
    checks.check(code != 0 and "-msd-memory-limit" in output and "MSD" in output,
                 "memory: the run stops at the MSD memory limit with a message",
                 f"exit {code}: {output[-500:]}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
