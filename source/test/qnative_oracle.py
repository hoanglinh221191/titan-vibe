"""Oracle for Q-native with several points per residue (STR-20, MDB-11/12).

A ten-residue LYS beta-hairpin in a 60 A box, five frames:
  0  the reference itself;
  1  strand 2 moved 1.6 A away from strand 1 (some contacts break);
  2  frame 0 translated across the x face of the box, coordinates wrapped
     atom by atom (minimum image must keep every contact);
  3  frame 1 wrapped the same way;
  4  strand 2 moved 3 A away (most contacts break).

For -qnative-mode sidechain-heavy, bb and ca the native pairs and Q(t) are
computed here: a pair (|i - j| >= 3 in the chain) is native when the minimum
minimum-image distance between the residues' contact atoms in frame 0 is at
most 6.5 A, and Q(t) is the fraction of native pairs within 6.5 A at frame t.
Q(0) must be 1. The old kernel indexed the per-atom point list with residue
indices, so in sidechain-heavy and bb mode it measured unrelated atoms.

The -qnative-ref case gives a reference without residue 1: the native
residues must be found in the trajectory by chain, number and name (the old
code applied reference list positions to the trajectory, one residue off).

  python -I test/qnative_oracle.py --binary titan.exe --work-root test/tmp/qnative_oracle
"""
import argparse
import itertools
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 60.0
CUT = 6.5
SIDE = [("CB", 1.53), ("CG", 3.05), ("CD", 4.55), ("CE", 6.05), ("NZ", 7.50)]
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def residue_atoms(resnr, x, y, up, strand2_shift=0.0):
    """LYS with its side chain along +z (up) or -z, backbone in the plane."""
    s = 1.0 if up else -1.0
    y = y + strand2_shift
    atoms = [("N", x - 1.2, y - 0.4, 0.0), ("H", x - 1.5, y - 1.3, 0.0),
             ("CA", x, y, 0.0), ("HA", x + 0.3, y + 0.9, -0.5 * s),
             ("C", x + 1.2, y - 0.5, 0.0), ("O", x + 1.4, y - 1.7, 0.0)]
    for k, (name, dz) in enumerate(SIDE):
        atoms.append((name, x + 0.3 * (k % 2), y + 0.2, s * dz))
    atoms += [("HB1", x + 0.9, y + 0.8, s * 1.5), ("HZ1", x + 0.2, y + 0.9, s * 8.0)]
    if resnr == 10:
        atoms.append(("OXT", x + 2.2, y - 0.2, 0.0))
    return atoms


def hairpin(shift2):
    """Residues 1-5 along +x at y = 0, residues 6-10 back along x at y = 5."""
    atoms = []
    for resnr in range(1, 11):
        k = resnr - 1 if resnr <= 5 else 10 - resnr
        x = 20.0 + 3.6 * k
        y = 30.0 if resnr <= 5 else 35.0
        up = k % 2 == 0
        for name, ax, ay, az in residue_atoms(resnr, x, y, up, shift2 if resnr > 5 else 0.0):
            atoms.append(("A", resnr, "LYS", name, ax, ay, 30.0 + az))
    return atoms


def translate_wrap(atoms, dx):
    return [(c, r, n, a, fx.wrap(x + dx, BOX), y, z) for c, r, n, a, x, y, z in atoms]


def frames():
    f0 = hairpin(0.0)
    f1 = hairpin(1.6)
    return [(0.0, (BOX, BOX, BOX), f0), (10.0, (BOX, BOX, BOX), f1),
            (20.0, (BOX, BOX, BOX), translate_wrap(f0, 33.0)),
            (30.0, (BOX, BOX, BOX), translate_wrap(f1, 31.0)),
            (40.0, (BOX, BOX, BOX), hairpin(3.0))]


def contact_points(atoms, mode):
    out = {}
    for chain, resnr, resname, name, x, y, z in atoms:
        heavy = not name.startswith("H")
        if mode == "sidechain-heavy":
            keep = heavy and name not in BACKBONE
        elif mode == "bb":
            keep = heavy and name in BACKBONE
        else:
            keep = name == "CA"
        if keep:
            out.setdefault(resnr, []).append((x, y, z))
    return out


def mic_distance(p, q):
    d = [p[k] - q[k] for k in range(3)]
    d = [v - BOX * round(v / BOX) for v in d]
    return math.sqrt(sum(v * v for v in d))


def min_distance(points, i, j):
    return min(mic_distance(p, q) for p in points[i] for q in points[j])


def expected(frames_, mode, ref_atoms):
    ref = contact_points(ref_atoms, mode)
    natives = [(i, j) for i, j in itertools.combinations(sorted(ref), 2)
               if j - i >= 3 and min_distance(ref, i, j) <= CUT]
    q = []
    for _, _, atoms in frames_:
        pts = contact_points(atoms, mode)
        q.append(sum(1 for i, j in natives if min_distance(pts, i, j) <= CUT) / len(natives))
    return natives, q


def read_q(case):
    files = sorted(case.glob("qnative_timeseries*.csv"))
    files = [f for f in files if "average" not in f.name]
    if not files:
        return None
    return [float(row["Q"]) for row in fx.read_csv(files[0])]


def read_natives(case):
    files = [f for f in sorted(case.glob("qnative_native_contacts*.csv")) if "average" not in f.name]
    if not files:
        return None
    pairs = []
    for row in fx.read_csv(files[0]):
        i = int(row["residue_i"].split(":")[1])
        j = int(row["residue_j"].split(":")[1])
        pairs.append(tuple(sorted((i, j))))
    return sorted(pairs)


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

    fr = frames()
    atoms0 = fr[0][2]
    fx.write_trajectory(work / "hairpin.pdb", fr)
    fx.write_topology(work / "hairpin.top", atoms0)
    fx.write_index(work / "hairpin.ndx", [("System", range(1, len(atoms0) + 1))])
    # reference without residue 1, same coordinates for residues 2-10
    ref_atoms = [atom for atom in atoms0 if atom[1] != 1]
    fx.write_trajectory(work / "ref_missing_res1.pdb", [(0.0, (BOX, BOX, BOX), ref_atoms)])

    checks = fx.Checks()
    base = fx.base_args(work, work / "hairpin.pdb", work / "hairpin.top", work / "hairpin.ndx")
    cases = [("sidechain_heavy", "sidechain-heavy", None, atoms0),
             ("bb", "bb", None, atoms0),
             ("ca", "ca", None, atoms0),
             ("ref_missing_res1", "sidechain-heavy", work / "ref_missing_res1.pdb", ref_atoms),
             # one point per residue: only the residue matching (MDB-12) is tested
             ("ref_missing_res1_ca", "ca", work / "ref_missing_res1.pdb", ref_atoms)]
    for name, mode, ref, ref_coords in cases:
        natives, q_expected = expected(fr, mode, ref_coords)
        lines = base + ["-qnative=true", f"-qnative-mode={mode}", f"-qnative-cut={CUT}",
                        "-qnative-min-seqsep=3", "-qnative-domain=both"]
        if ref is not None:
            lines.append("-qnative-ref=" + ref.as_posix())
        code, output = fx.run_titan(binary, work / name, lines)
        if not checks.check(code == 0, f"{name}: titan exit 0", f"exit {code}: {output[-600:]}" if code else ""):
            continue
        got_natives = read_natives(work / name)
        checks.check(got_natives == sorted(natives), f"{name}: native pairs",
                     f"titan {got_natives}, expected {sorted(natives)}")
        q = read_q(work / name)
        ok = q is not None and len(q) == len(q_expected) and all(
            abs(g - e) <= 2e-6 for g, e in zip(q, q_expected))
        checks.check(ok, f"{name}: Q(t) for {len(natives)} native pairs",
                     f"titan {q}, expected {[round(v, 6) for v in q_expected]}")
        if q:
            checks.check(abs(q[0] - 1.0) <= 1e-6, f"{name}: Q(0) = 1 for the reference itself",
                         f"Q(0) = {q[0]}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
