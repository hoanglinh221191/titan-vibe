"""Oracle for the MSD diffusion uncertainty (STR-25, MDA-14).

A rigid three-atom molecule whose centre of mass does free Brownian motion,
D = 0.05 A^2/ps, dt = 1 ps, 2000 frames (fixed seeds), in a 300 A box;
-msd-fit-window=10,50, -uncertainty=block, -uncertainty-blocks=10.

1. Time blocks: msd_diffusion_summary must list ten 200-ps blocks and, per
   block, the MSD slope fitted over 10-50 ps from origins inside the block;
   the slopes are recomputed here from the same coordinates (1e-6 relative).
2. msd_uncertainty.csv, D_A2_per_ps: n_blocks = 10 and sem = sd(D_b)/sqrt(10)
   of those per-block values; the row says it is a time-block estimate.
3. Calibration: the mean reported sem of eight runs (other seeds) must match
   the true sampling error of D for this trajectory length, the standard
   deviation of D over 32 independent walks fitted here the same way
   (0.6-1.6x; the time-block sem of a single run scatters by ~25 %). The old
   "block" sem split the lag points of one MSD curve (10 to 50 ps) into
   blocks and measured how straight the curve is, not the sampling error.

Runs the energy path (-g=1-2, a fixed one-residue chain B as molecule 2):
in -nocal runs the automatic -pbci 0 re-centres molecule 1 in the box every
frame, so the MSD of a lone molecule's centre of mass is 0 there.

  python -I test/msd_time_block_oracle.py --binary titan.exe --work-root test/tmp/msd_time_block_oracle
"""
import argparse
import math
import random
import shutil
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 300.0
D_TRUE = 0.05
FRAMES = 2000
DT = 1.0
BLOCKS = 10
TMIN, TMAX = 10.0, 50.0
BODY = [(0.0, 0.0, 0.0), (1.4, 0.0, 0.0), (0.0, 1.4, 0.0)]


def walk(seed, frames=FRAMES):
    rng = random.Random(seed)
    s = math.sqrt(2.0 * D_TRUE * DT)
    p = [150.0, 150.0, 150.0]
    out = []
    for _ in range(frames):
        out.append(tuple(p))
        p = [p[k] + rng.gauss(0.0, s) for k in range(3)]
    return out


def rounded_com(centre):
    atoms = [tuple(round(centre[k] + b[k], 3) for k in range(3)) for b in BODY]
    return tuple(sum(a[k] for a in atoms) / len(atoms) for k in range(3))


def slope(com, f0, f1):
    """MSD slope over TMIN..TMAX ps from origins and targets inside [f0, f1)."""
    xs, ys = [], []
    for lag in range(1, f1 - f0):
        t = lag * DT
        if t < TMIN or t > TMAX:
            continue
        acc = 0.0
        n = 0
        for o in range(f0, f1 - lag):
            acc += sum((com[o + lag][k] - com[o][k]) ** 2 for k in range(3))
            n += 1
        xs.append(t)
        ys.append(acc / n)
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def run_case(binary, work, seed):
    """Writes the walk of this seed and runs titan: (centres, run dir, exit, output)."""
    case = work / f"seed{seed}"
    case.mkdir()
    centres = walk(seed)
    frames = []
    for i, c in enumerate(centres):
        atoms = [("A", 1, "MOL", f"C{k + 1}", c[0] + b[0], c[1] + b[1], c[2] + b[2])
                 for k, b in enumerate(BODY)]
        atoms += [("B", 1, "GLY", n, x, 20.0, 20.0)
                  for n, x in (("N", 20.0), ("CA", 21.4), ("C", 22.6), ("O", 23.8))]
        frames.append((DT * i, (BOX, BOX, BOX), atoms))
    fx.write_trajectory(case / "walk.pdb", frames)
    fx.write_topology(case / "walk.top", frames[0][2])
    fx.write_index(case / "walk.ndx", [("System", range(1, 8)), ("mol", [1, 2, 3]),
                                       ("fixed", [4, 5, 6, 7])])
    lines = [x for x in fx.base_args(work, case / "walk.pdb", case / "walk.top",
                                     case / "walk.ndx", groups="1-2") if x != "-nocal=true"]
    lines += ["-msd=true", "-msd-mode=com", f"-msd-fit-window={TMIN:g},{TMAX:g}",
              "-uncertainty=block", f"-uncertainty-blocks={BLOCKS}"]
    code, output = fx.run_titan(binary, case / "run", lines, timeout=1800)
    return centres, case / "run", code, output


def d_row(run):
    sidecar = run / "msd_uncertainty.csv"
    rows = {r["observable"]: r for r in fx.read_csv(sidecar)} if sidecar.exists() else {}
    return rows.get("D_A2_per_ps", {})


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

    centres, run, code, output = run_case(binary, work, 2025)
    if not checks.check(code == 0, "titan exit 0", f"exit {code}: {output[-800:]}" if code else ""):
        return checks.exit_code()
    com = [rounded_com(c) for c in centres]
    per = FRAMES // BLOCKS
    want = [slope(com, b * per, (b + 1) * per) for b in range(BLOCKS)]
    summary = [f for f in run.glob("msd_diffusion_summary*.dat") if "average" not in f.name]
    kv = fx.read_key_values(summary[0]) if summary else {}
    got = [float(v) for v in kv.get("time_block_slope_A2_per_ps", "").split(",") if v.strip()]
    ok = len(got) == BLOCKS and all(abs(g - w) <= 1e-6 * abs(w) + 1e-9 for g, w in zip(got, want))
    checks.check(kv.get("time_blocks") == str(BLOCKS) and ok,
                 "time-block MSD slopes = slopes recomputed from the coordinates",
                 f"blocks {kv.get('time_blocks')}, length {kv.get('time_block_length_ps')} ps, "
                 f"titan {got[:3]}..., here {[round(w, 6) for w in want[:3]]}...")

    d_blocks = [w / 6.0 for w in want]
    sem_want = statistics.stdev(d_blocks) / math.sqrt(BLOCKS)
    row = d_row(run)
    sem_got = float(row.get("sem", "nan"))
    checks.check(row.get("n_blocks") == str(BLOCKS) and abs(sem_got - sem_want) <= 1e-6 * sem_want
                 and "time_block" in row.get("note", ""),
                 "D_A2_per_ps: sem over the ten time blocks",
                 f"titan n_blocks {row.get('n_blocks')}, sem {sem_got:.6g}, note '{row.get('note')}'; "
                 f"here sem {sem_want:.6g}")

    # calibration over an ensemble: the mean reported sem of eight runs against
    # the spread of the full-trajectory D over 32 independent walks
    sems = []
    for seed in range(300, 308):
        _, run, code, _ = run_case(binary, work, seed)
        sems.append(float(d_row(run).get("sem", "nan")) if code == 0 else float("nan"))
    d_ind = [slope(walk(seed), 0, FRAMES) / 6.0 for seed in range(400, 432)]
    true_sd = statistics.stdev(d_ind)
    mean_sem = statistics.mean(sems)
    checks.check(math.isfinite(mean_sem) and 0.6 * true_sd <= mean_sem <= 1.6 * true_sd,
                 "mean reported sem of 8 runs = sampling sd of D over 32 walks (0.6-1.6x)",
                 f"mean sem {mean_sem:.4g} (runs {[round(v, 4) for v in sems]}), "
                 f"sd {true_sd:.4g}, mean D {statistics.mean(d_ind):.4g} (true {D_TRUE})")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
