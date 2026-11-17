"""Oracle for -s2 with -s2-align=rms when the protein crosses a box face
(STR-24, MDA-09).

A rigid 20-residue helix (Calpha trace ~30 A long) in a 50 A box. Every frame
applies a random rotation (fixed seed) and moves the centre by 2.3 A along x,
so from frame ~5 on the protein straddles the x face; all atoms are wrapped
into the box one by one. The body is rigid, so after a correct fit every
NH and CA-C bond vector is constant: C(t) = 1 and S2 = 1 for every vector.

The old fit took the raw Calpha coordinates: a trace split by the face gives
a meaningless rotation and C(t) decays, S2 far below 1 with fit_status "ok".
A second run without crossing (frames kept whole) must give S2 = 1 on both
builds: the oracle checks the fit itself, not only the imaging.

The runs use the energy path (-g=1-2, a one-residue chain B far away as
molecule 2): -nocal switches on -pbci 0, which makes molecule 1 whole before
-s2 sees it and hides the bug; without it S2 reads the coordinates as given.

  python -I test/s2_pbc_oracle.py --binary titan.exe --work-root test/tmp/s2_pbc_oracle
"""
import argparse
import math
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 50.0
NRES = 20
NAMES = ["ALA", "LEU", "SER", "GLU", "LYS"]


def helix():
    atoms = []
    for i in range(NRES):
        a = math.radians(100.0 * i)
        z = 1.5 * i - 14.0
        ca = (2.3 * math.cos(a), 2.3 * math.sin(a), z)
        n = (1.6 * math.cos(a - 0.45), 1.6 * math.sin(a - 0.45), z - 0.7)
        c = (1.9 * math.cos(a + 0.5), 1.9 * math.sin(a + 0.5), z + 0.8)
        h = (n[0] * 0.5, n[1] * 0.5, n[2] - 0.6)
        o = (c[0] * 1.5, c[1] * 1.5, c[2] + 0.3)
        cb = (3.6 * math.cos(a), 3.6 * math.sin(a), z + 0.2)
        name = NAMES[i % len(NAMES)]
        for atom, p in (("N", n), ("H", h), ("CA", ca), ("CB", cb), ("C", c), ("O", o)):
            atoms.append(("A", i + 1, name, atom, p))
    return atoms


def random_rotation(rng):
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    q = (math.sqrt(1 - u1) * math.sin(2 * math.pi * u2), math.sqrt(1 - u1) * math.cos(2 * math.pi * u2),
         math.sqrt(u1) * math.sin(2 * math.pi * u3), math.sqrt(u1) * math.cos(2 * math.pi * u3))
    w, x, y, z = q
    return [[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]


def frames(crossing, count=40):
    rng = random.Random(24)
    body = helix()
    out = []
    for f in range(count):
        r = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] if f == 0 else random_rotation(rng)
        cx = 25.0 + (2.3 * f if crossing else 0.0)
        atoms = []
        for chain, resnr, resname, name, p in body:
            q = [sum(r[i][k] * p[k] for k in range(3)) for i in range(3)]
            x, y, z = cx + q[0], 25.0 + q[1], 25.0 + q[2]
            atoms.append((chain, resnr, resname, name, fx.wrap(x, BOX), fx.wrap(y, BOX), fx.wrap(z, BOX)))
        # molecule 2: one GLY residue, fixed, away from the helix's path
        for name, (x, y, z) in (("N", (25.0, 5.0, 5.0)), ("CA", (26.4, 5.2, 5.0)),
                                ("C", (27.0, 4.0, 5.1)), ("O", (26.4, 3.0, 5.3))):
            atoms.append(("B", 1, "GLY", name, x, y, z))
        out.append((10.0 * f, (BOX, BOX, BOX), atoms))
    return out


def order_rows(case):
    files = [f for f in sorted(case.glob("s2_order_parameters*.csv")) if "average" not in f.name]
    return fx.read_csv(files[0]) if files else []


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
    for name, crossing in (("whole", False), ("crossing", True)):
        case = work / name
        case.mkdir()
        fr = frames(crossing)
        fx.write_trajectory(case / "helix.pdb", fr)
        fx.write_topology(case / "helix.top", fr[0][2])
        n_a = 6 * NRES
        fx.write_index(case / "helix.ndx", [("System", range(1, len(fr[0][2]) + 1)),
                                            ("chA", range(1, n_a + 1)),
                                            ("chB", range(n_a + 1, len(fr[0][2]) + 1))])
        lines = [a for a in fx.base_args(work, case / "helix.pdb", case / "helix.top",
                                         case / "helix.ndx", groups="1-2") if a != "-nocal=true"]
        lines += ["-s2=true", "-s2-align=rms", "-s2-vector=both"]
        code, output = fx.run_titan(binary, case / "run", lines)
        if not checks.check(code == 0, f"{name}: titan exit 0", f"exit {code}: {output[-800:]}" if code else ""):
            continue
        rows = order_rows(case / "run")
        s2 = [float(r["S2_clipped"]) for r in rows]
        plateau = [float(r["plateau_C_raw"]) for r in rows]
        checks.check(len(rows) >= 2 * (NRES - 1),
                     f"{name}: one row per NH and CA-C vector", f"{len(rows)} rows")
        # coordinates carry 3 decimals: a constant unit vector up to ~1e-3
        checks.check(s2 and min(s2) > 0.995 and min(plateau) > 0.995,
                     f"{name}: S2 = 1 for every vector of a rigid body",
                     f"min S2 {min(s2) if s2 else None}, min plateau {min(plateau) if plateau else None}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
