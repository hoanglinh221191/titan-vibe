"""Oracle for -msd-rot (STR-23, MDA-05...08).

A rigid six-atom body (+-3, +-2, +-1 A along its axes) at the centre of a
40 A box rotates by Brownian steps drawn in the BODY frame,
R(t + dt) = R(t) Exp(w), w_i ~ N(0, 2 D_i dt), with a fixed seed.

1. Anisotropic body, D = (0.004, 0.010, 0.025) rad^2/ps, dt = 1 ps, 3000
   frames: the fitted principal diffusion coefficients must match D within
   15 % and the principal axes must be the body axes (|cos| >= 0.95). A
   lab-frame increment would mix the axes and pull all three towards D_iso.
2. The same trajectory: kabsch_rmsd_A of a rigid body must be 0 (< 1e-4 A).
3. Isotropic body, D = 0.01 rad^2/ps: the exported C1 and C2 must be the
   Legendre correlations of a body vector averaged over directions,
   <(1 + 2 cos t)/3> and <(4 cos^2 t + 2 cos t - 1)/5> of the rotation angle
   t of each origin pair, computed here from the same frames, and must follow
   exp(-2 D t) and exp(-6 D t) at the first lags. The old columns were
   <cos t> and <P2(cos t)> of the rotation angle (C1 decayed 1.5x too fast).

  python -I test/rotational_diffusion_oracle.py --binary titan.exe --work-root test/tmp/rot_oracle
"""
import argparse
import math
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 40.0
BODY = [(3.0, 0.0, 0.0), (-3.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, -2.0, 0.0),
        (0.0, 0.0, 1.0), (0.0, 0.0, -1.0)]
NAMES = ["C1", "C2", "C3", "C4", "C5", "C6"]


def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + w2 * x1 + y1 * z2 - z1 * y2,
            w1 * y2 + w2 * y1 + z1 * x2 - x1 * z2,
            w1 * z2 + w2 * z1 + x1 * y2 - y1 * x2)


def qexp(w):
    angle = math.sqrt(sum(v * v for v in w))
    if angle < 1e-15:
        return (1.0, 0.0, 0.0, 0.0)
    s = math.sin(0.5 * angle) / angle
    return (math.cos(0.5 * angle), w[0] * s, w[1] * s, w[2] * s)


def qnorm(q):
    n = math.sqrt(sum(v * v for v in q))
    return tuple(v / n for v in q)


def rotate(q, v):
    w, x, y, z = q
    r = [[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
         [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
         [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]
    return [sum(r[i][k] * v[k] for k in range(3)) for i in range(3)]


def simulate(d, frames, dt, seed):
    rng = random.Random(seed)
    q = (1.0, 0.0, 0.0, 0.0)
    out = []
    for _ in range(frames):
        out.append(q)
        w = [rng.gauss(0.0, math.sqrt(2.0 * d[i] * dt)) for i in range(3)]
        q = qnorm(qmul(q, qexp(w)))  # body-frame step: right multiplication
    return out


def write_case(work, name, quats, dt):
    frames = []
    for i, q in enumerate(quats):
        atoms = []
        for n, b in zip(NAMES, BODY):
            p = rotate(q, b)
            atoms.append(("A", 1, "BOD", n, 20.0 + p[0], 20.0 + p[1], 20.0 + p[2]))
        frames.append((dt * i, (BOX, BOX, BOX), atoms))
    case = work / name
    case.mkdir(parents=True)
    fx.write_trajectory(case / "body.pdb", frames)
    fx.write_topology(case / "body.top", frames[0][2])
    fx.write_index(case / "body.ndx", [("System", range(1, len(BODY) + 1))])
    return case


def run(binary, work, case, maxlag):
    lines = fx.base_args(work, case / "body.pdb", case / "body.top", case / "body.ndx")
    lines += ["-msd=true", "-msd-rot=true", f"-msd-maxlag={maxlag}"]
    return fx.run_titan(binary, case / "run", lines, timeout=1800)


def rotational_files(case):
    run_dir = case / "run"
    summary = [f for f in run_dir.glob("*rot*summary*.dat") if "average" not in f.name]
    msd = [f for f in run_dir.glob("*rot*msd*.dat") if "average" not in f.name]
    quat = [f for f in run_dir.glob("*quaternion*.csv") if "average" not in f.name]
    return (summary[0] if summary else None, msd[0] if msd else None, quat[0] if quat else None)


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
    (work / "lib").mkdir()
    checks = fx.Checks()

    # 1-2: anisotropic body
    d_true = (0.004, 0.010, 0.025)
    quats = simulate(d_true, 3000, 1.0, 20261022)
    case = write_case(work, "anisotropic", quats, 1.0)
    code, output = run(binary, work, case, 40)
    if checks.check(code == 0, "anisotropic: titan exit 0", f"exit {code}: {output[-800:]}" if code else ""):
        summary_file, _, quat_file = rotational_files(case)
        kv = fx.read_key_values(summary_file) if summary_file else {}
        d_fit = [float(kv.get(f"D{i}_rad2_per_ps", "nan")) for i in (1, 2, 3)]
        axes = [[float(v) for v in kv.get(f"principal_axis_{i}", "nan,nan,nan").split(",")]
                for i in (1, 2, 3)]
        order = sorted(range(3), key=lambda i: d_fit[i])
        ok = all(abs(d_fit[order[k]] - d_true[k]) <= 0.15 * d_true[k] for k in range(3))
        checks.check(ok, "anisotropic: fitted D1 < D2 < D3 within 15 % of the body-frame D",
                     f"fit {sorted(d_fit)}, true {d_true}, status {kv.get('fit_status')}, "
                     f"window {kv.get('fit_window_policy')}")
        axis_ok = all(abs(axes[order[k]][k]) >= 0.95 for k in range(3))
        checks.check(axis_ok, "anisotropic: principal axes = body axes x, y, z",
                     f"axes by increasing D {[axes[i] for i in order]}")
        rmsd = [float(r["kabsch_rmsd_A"]) for r in fx.read_csv(quat_file)] if quat_file else []
        # coordinates carry 3 decimals: a rigid body fits to ~5e-4 A
        checks.check(rmsd and max(rmsd) < 2e-3, "anisotropic: kabsch_rmsd_A ~ 0 for a rigid body",
                     f"max {max(rmsd) if rmsd else None} A")

    # 3: isotropic body, C1/C2
    d_iso = 0.01
    quats = simulate((d_iso, d_iso, d_iso), 2000, 1.0, 7)
    case = write_case(work, "isotropic", quats, 1.0)
    code, output = run(binary, work, case, 20)
    if checks.check(code == 0, "isotropic: titan exit 0", f"exit {code}: {output[-800:]}" if code else ""):
        _, msd_file, quat_file = rotational_files(case)
        lines = msd_file.read_text().splitlines() if msd_file else []
        header = lines[0].split() if lines else []
        table = [dict(zip(header, l.split())) for l in lines[1:] if l.strip()]
        # the module's own orientations (8 decimals), so that the reference
        # tests the C1/C2 definition and not the 3-decimal PDB rounding
        if quat_file:
            quats = [tuple(float(r[k]) for k in ("qw", "qx", "qy", "qz"))
                     for r in fx.read_csv(quat_file)]
        worst1 = worst2 = 0.0
        model1 = model2 = 0.0
        for row in table[:20]:
            lag = int(row["lag_index"])
            c1 = c2 = 0.0
            n = 0
            for o in range(len(quats) - lag):
                q0, q1 = quats[o], quats[o + lag]
                rel = qmul((q0[0], -q0[1], -q0[2], -q0[3]), q1)
                c = math.cos(2.0 * math.acos(min(1.0, abs(rel[0]))))
                c1 += (1.0 + 2.0 * c) / 3.0
                c2 += (4.0 * c * c + 2.0 * c - 1.0) / 5.0
                n += 1
            c1 /= n
            c2 /= n
            worst1 = max(worst1, abs(float(row["C1"]) - c1))
            worst2 = max(worst2, abs(float(row["C2"]) - c2))
            if lag <= 10:
                t = float(row["lag_time_ps"])
                model1 = max(model1, abs(float(row["C1"]) - math.exp(-2 * d_iso * t)))
                model2 = max(model2, abs(float(row["C2"]) - math.exp(-6 * d_iso * t)))
        checks.check(len(table) >= 20 and worst1 < 1e-6 and worst2 < 1e-6,
                     "isotropic: C1, C2 = direction-averaged Legendre correlations",
                     f"max |C1 - ref| {worst1:.3g}, max |C2 - ref| {worst2:.3g}")
        checks.check(model1 < 0.03 and model2 < 0.05,
                     "isotropic: C1 ~ exp(-2Dt), C2 ~ exp(-6Dt) over lags 1-10",
                     f"max dev {model1:.3g}, {model2:.3g}")

    # 4: a 30 A zigzag rod in the 40 A box (MDA-07): its far atoms are more
    # than half a box from atom 0, so the old first-frame reconstruction
    # (every atom imaged next to atom 0) broke it and the fit with it
    global BODY, NAMES
    BODY = [(-15.0 + 3.0 * k, 1.0 if k % 2 else -1.0, 0.5 * (k % 3 - 1)) for k in range(11)]
    NAMES = [f"C{k + 1}" for k in range(11)]
    quats = simulate((d_iso, d_iso, d_iso), 1500, 1.0, 11)
    case = write_case(work, "long_body", quats, 1.0)
    code, output = run(binary, work, case, 20)
    if checks.check(code == 0, "long body: titan exit 0", f"exit {code}: {output[-800:]}" if code else ""):
        summary_file, _, quat_file = rotational_files(case)
        kv = fx.read_key_values(summary_file) if summary_file else {}
        rmsd = [float(r["kabsch_rmsd_A"]) for r in fx.read_csv(quat_file)] if quat_file else []
        checks.check(rmsd and max(rmsd) < 2e-3, "long body: whole in every frame (kabsch_rmsd_A ~ 0)",
                     f"max {max(rmsd) if rmsd else None} A")
        d_fit = float(kv.get("D_iso_rad2_per_ps", "nan"))
        checks.check(abs(d_fit - d_iso) <= 0.15 * d_iso, "long body: D_iso within 15 %",
                     f"D_iso {d_fit}, true {d_iso}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
