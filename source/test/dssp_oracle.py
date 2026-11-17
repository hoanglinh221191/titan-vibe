"""Oracle for the DSSP module (STR-10, MDB-06/07; STR-26 MDB-08/09).

Backbones are built here from ideal internal coordinates (NeRF) and the
secondary structure is computed by an independent Kabsch-Sander / mkdssp
implementation in this file, then compared residue by residue with titan's
dssp_mol1.dat and with its bridge list:

  helix_gap       an ideal alpha helix of 20 ALA, heavy atoms only, with
                  residue 11 deleted (numbering kept). C(10)-N(12) is 3.5 A,
                  a chain break: no helix, turn or rebuilt H may span it. The
                  old code tested breaks by a sequential index, so the helix
                  ran through the gap (i -> i+4 bonds became 3-turns).
  marginal_helix  a uniform helix whose i -> i+4 bond energy straddles the
                  -0.5 kcal/mol cutoff: below it with the old N + (C - O) H
                  (N-H 1.23 A), above it with mkdssp's N + (C - O)/|C - O|
                  (N-H 1.0 A).
  explicit_h      marginal_helix with the H atoms written at N + unit(C - O):
                  the rebuilt H must give the same answer as the explicit H.
  tight_turn      a 3-10-like segment in which the old bridge search from
                  j = i + 1 finds bridges (i, i + 2); mkdssp starts at i + 3.
  bulge_sign      a strand paired antiparallel with two partner strands where
                  the second partner lies after the first in sequence: the old
                  signed bulge test (jbi - jej < 6 with jbi - jej < 0) merged
                  the two ladders and the partner residues lost their E state;
                  mkdssp computes the differences unsigned and keeps them apart.

For each case the oracle also checks that the old rules give a different
answer, so the case keeps testing the fix.

  python -I test/dssp_oracle.py --binary titan.exe --work-root D:/softs/titan_work/str/dssp_oracle
"""
import argparse
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 200.0

# ---------------------------------------------------------------- geometry


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    return scale(a, 1.0 / n)


def dist(a, b):
    return norm(sub(a, b))


def place(a, b, c, bond, angle_deg, torsion_deg):
    """NeRF: D with |CD| = bond, angle BCD, torsion ABCD."""
    angle = math.radians(angle_deg)
    torsion = math.radians(torsion_deg)
    bc = unit(sub(c, b))
    n = unit(cross(sub(b, a), bc))
    m = cross(n, bc)
    d = (-bond * math.cos(angle), bond * math.sin(angle) * math.cos(torsion),
         bond * math.sin(angle) * math.sin(torsion))
    return add(c, add(scale(bc, d[0]), add(scale(m, d[1]), scale(n, d[2]))))


def build_backbone(phi_psi, omega=180.0):
    """Residues {N, CA, C, O} from (phi, psi) per residue, ideal bonds/angles."""
    n0 = (0.0, 0.0, 0.0)
    ca0 = (1.458, 0.0, 0.0)
    ang = math.radians(111.2)
    c0 = add(ca0, (-1.525 * math.cos(ang), 1.525 * math.sin(ang), 0.0))
    atoms = [n0, ca0, c0]
    for k in range(1, len(phi_psi)):
        psi_prev = phi_psi[k - 1][1]
        phi = phi_psi[k][0]
        n = place(atoms[-3], atoms[-2], atoms[-1], 1.329, 116.2, psi_prev)
        ca = place(atoms[-2], atoms[-1], n, 1.458, 121.7, omega)
        c = place(atoms[-1], n, ca, 1.525, 111.2, phi)
        atoms += [n, ca, c]
    residues = []
    for k in range(len(phi_psi)):
        n, ca, c = atoms[3 * k:3 * k + 3]
        if k + 1 < len(phi_psi):
            # O in the peptide plane, opposite the next N
            o = place(atoms[3 * k + 3], ca, c, 1.231, 120.5, 180.0)
        else:
            o = place(n, ca, c, 1.231, 120.5, phi_psi[k][1] + 180.0)
        residues.append({"N": n, "CA": ca, "C": c, "O": o})
    return residues


def rotation(axis, angle_deg):
    x, y, z = unit(axis)
    a = math.radians(angle_deg)
    c, s, t = math.cos(a), math.sin(a), 1.0 - math.cos(a)
    return ((t * x * x + c, t * x * y - s * z, t * x * z + s * y),
            (t * x * y + s * z, t * y * y + c, t * y * z - s * x),
            (t * x * z - s * y, t * y * z + s * x, t * z * z + c))


def transform(residues, rot=None, shift=(0.0, 0.0, 0.0)):
    out = []
    for res in residues:
        new = {}
        for name, p in res.items():
            if rot is not None:
                p = (dot(rot[0], p), dot(rot[1], p), dot(rot[2], p))
            new[name] = add(p, shift)
        out.append(new)
    return out


# ------------------------------------------------------- reference DSSP

CUTOFF = -0.5
STATES = "HBEGITSC"


def round_half_away(x):
    return math.copysign(math.floor(abs(x) * 1000.0 + 0.5), x) / 1000.0


def hbond_energy(donor, acceptor):
    """Kabsch-Sander energy (kcal/mol) as mkdssp: rounded to 1e-3, >= -9.9."""
    n, h = donor["N"], donor["H"]
    c, o = acceptor["C"], acceptor["O"]
    d_ho, d_hc, d_nc, d_no = dist(h, o), dist(h, c), dist(n, c), dist(n, o)
    if min(d_ho, d_hc, d_nc, d_no) < 0.5:
        return -9.9
    e = -27.888 / d_ho + 27.888 / d_hc - 27.888 / d_nc + 27.888 / d_no
    return max(-9.9, round_half_away(e))


LEGACY_RULES = ("h", "breaks", "bridge_start", "bulge")


class Dssp:
    """mkdssp rules as titan implements them. legacy names the old titan rules
    to use instead (True = all): "h" unnormalised H = N + (C - O), also from
    the previous chain; "breaks" chain breaks by chain only, not by
    C(i)-N(i+1) > 2.5 A; "bridge_start" bridge search from j = i + 1;
    "bulge" signed bulge differences."""

    def __init__(self, residues, legacy=()):
        self.res = [dict(r) for r in residues]
        self.legacy = set(LEGACY_RULES) if legacy is True else set(legacy)
        n = len(self.res)
        self.n = n
        self.breaks = [False] * max(0, n - 1)
        for i in range(n - 1):
            a, b = self.res[i], self.res[i + 1]
            if a["chain"] != b["chain"]:
                self.breaks[i] = True
            elif "breaks" not in self.legacy:
                self.breaks[i] = ("C" not in a or "N" not in b or
                                  dist(a["C"], b["N"]) > 2.5)
        self.prefix = [0]
        for b in self.breaks:
            self.prefix.append(self.prefix[-1] + (1 if b else 0))
        for i in range(1, n):
            r, p = self.res[i], self.res[i - 1]
            if "H" in r or "N" not in r or r["resname"] == "PRO":
                continue
            if "C" not in p or "O" not in p:
                continue
            if "h" in self.legacy:
                r["H"] = add(r["N"], sub(p["C"], p["O"]))
            elif not self.breaks[i - 1]:
                r["H"] = add(r["N"], unit(sub(p["C"], p["O"])))
        self.best = [[] for _ in range(n)]
        for d in range(n):
            rd = self.res[d]
            if "N" not in rd or "H" not in rd or "CA" not in rd or rd["resname"] == "PRO":
                continue
            best = [(-1, 0.0), (-1, 0.0)]
            for a in range(n):
                ra = self.res[a]
                if "CA" not in ra or "C" not in ra or "O" not in ra:
                    continue
                if a == d or a + 1 == d:
                    continue
                if dist(rd["CA"], ra["CA"]) > 9.0:
                    continue
                e = hbond_energy(rd, ra)
                if e < best[0][1]:
                    best = [(a, e), best[0]]
                elif e < best[1][1]:
                    best[1] = (a, e)
            self.best[d] = best

    def no_break(self, i, j):
        if i > j:
            i, j = j, i
        if i < 0 or j >= self.n:
            return False
        return self.prefix[j] == self.prefix[i]

    def bond(self, donor, acceptor):
        if donor < 0 or donor >= self.n:
            return False
        return any(a == acceptor and e < CUTOFF for a, e in self.best[donor])

    def bridge_type(self, r1, r2):
        n = self.n
        if r1 == 0 or r2 == 0 or r1 + 1 >= n or r2 + 1 >= n:
            return None
        if not self.no_break(r1 - 1, r1 + 1) or not self.no_break(r2 - 1, r2 + 1):
            return None
        b = self.bond
        if (b(r1 + 1, r2) and b(r2, r1 - 1)) or (b(r2 + 1, r1) and b(r1, r2 - 1)):
            return "parallel"
        if (b(r1 + 1, r2 - 1) and b(r2 + 1, r1 - 1)) or (b(r2, r1) and b(r1, r2)):
            return "antiparallel"
        return None

    def assign(self):
        n = self.n
        states = ["C"] * n
        pairs = []
        bridges = []
        first_j = 1 if "bridge_start" in self.legacy else 3
        for i in range(1, n - 1):
            if "CA" not in self.res[i]:
                continue
            for j in range(i + first_j, n - 1):
                if "CA" not in self.res[j]:
                    continue
                if dist(self.res[i]["CA"], self.res[j]["CA"]) > 9.0:
                    continue
                t = self.bridge_type(i, j)
                if t is None:
                    continue
                pairs.append((i, j, t))
                for br in bridges:
                    if br["type"] != t or i != br["i"][-1] + 1:
                        continue
                    if t == "parallel" and br["j"][-1] + 1 == j:
                        br["i"].append(i)
                        br["j"].append(j)
                        break
                    if t == "antiparallel" and br["j"][0] == j + 1:
                        br["i"].append(i)
                        br["j"].insert(0, j)
                        break
                else:
                    bridges.append({"type": t, "i": [i], "j": [j]})
        bridges.sort(key=lambda br: br["i"][0])

        def diff(a, b):
            # mkdssp: uint32 arithmetic, a negative difference is huge
            if "bulge" in self.legacy:
                return a - b
            return a - b if a >= b else 1 << 40

        k = 0
        while k < len(bridges):
            m = k + 1
            while m < len(bridges):
                bi, bj = bridges[k], bridges[m]
                ibi, iei, jbi, jei = bi["i"][0], bi["i"][-1], bi["j"][0], bi["j"][-1]
                ibj, iej, jbj, jej = bj["i"][0], bj["i"][-1], bj["j"][0], bj["j"][-1]
                if (bi["type"] != bj["type"] or
                        not self.no_break(min(ibi, ibj), max(iei, iej)) or
                        not self.no_break(min(jbi, jbj), max(jei, jej)) or
                        diff(ibj, iei) >= 6 or (iei >= ibj and ibi <= iej)):
                    m += 1
                    continue
                if bi["type"] == "parallel":
                    bulge = (diff(jbj, jei) < 6 and diff(ibj, iei) < 3) or diff(jbj, jei) < 3
                else:
                    bulge = (diff(jbi, jej) < 6 and diff(ibj, iei) < 3) or diff(jbi, jej) < 3
                if not bulge:
                    m += 1
                    continue
                bi["i"] += bj["i"]
                if bi["type"] == "parallel":
                    bi["j"] += bj["j"]
                else:
                    bi["j"] = bj["j"] + bi["j"]
                del bridges[m]
            k += 1
        for br in bridges:
            ss = "E" if len(br["i"]) > 1 else "B"
            for idx in range(br["i"][0], br["i"][-1] + 1):
                if states[idx] != "E":
                    states[idx] = ss
            for idx in range(br["j"][0], br["j"][-1] + 1):
                if states[idx] != "E":
                    states[idx] = ss

        flags = [[None, None, None] for _ in range(n)]
        for h, stride in enumerate((3, 4, 5)):
            for i in range(n - stride):
                if not self.no_break(i, i + stride) or not self.bond(i + stride, i):
                    continue
                flags[i + stride][h] = "end"
                for j in range(i + 1, i + stride):
                    if flags[j][h] is None:
                        flags[j][h] = "middle"
                flags[i][h] = "startend" if flags[i][h] == "end" else "start"

        def start(i, h):
            return flags[i][h] in ("start", "startend")

        for i in range(1, n - 4):
            if start(i, 1) and start(i - 1, 1):
                for j in range(i, i + 4):
                    states[j] = "H"
        for i in range(1, n - 3):
            if start(i, 0) and start(i - 1, 0):
                if all(states[j] in "CG" for j in range(i, i + 3)):
                    for j in range(i, i + 3):
                        states[j] = "G"
        for i in range(1, n - 5):
            if start(i, 2) and start(i - 1, 2):
                if all(states[j] in "CIH" for j in range(i, i + 5)):
                    for j in range(i, i + 5):
                        states[j] = "I"
        for i in range(1, n - 1):
            if states[i] != "C":
                continue
            turn = any(i >= k and start(i - k, h)
                       for h, stride in enumerate((3, 4, 5)) for k in range(1, stride))
            if turn:
                states[i] = "T"
            elif 2 <= i < n - 2 and self.kappa(i) > 70.0:
                states[i] = "S"
        return "".join(states), sorted(pairs)

    def kappa(self, i):
        r = self.res
        if not all("CA" in r[k] for k in (i - 2, i, i + 2)) or not self.no_break(i - 2, i + 2):
            return 0.0
        v1 = sub(r[i]["CA"], r[i - 2]["CA"])
        v2 = sub(r[i + 2]["CA"], r[i]["CA"])
        c = max(-1.0, min(1.0, dot(v1, v2) / (norm(v1) * norm(v2))))
        return math.degrees(math.acos(c))


# ------------------------------------------------------------ test cases


def label(residues, chain="A", first=1, skip=()):
    out = []
    number = first
    for res in residues:
        while number in skip:
            number += 1
        item = dict(res)
        item.update(chain=chain, resnr=number, resname="ALA")
        out.append(item)
        number += 1
    return out


def helix(count, phi=-57.8, psi=-47.0):
    return build_backbone([(phi, psi)] * count)


def case_helix_gap():
    res = label(helix(20))
    del res[10]  # residue 11
    return res


def case_marginal_helix():
    return label(helix(16, MARGINAL[0], MARGINAL[1]))


def case_explicit_h():
    res = case_marginal_helix()
    for i in range(1, len(res)):
        res[i]["H"] = add(res[i]["N"], unit(sub(res[i - 1]["C"], res[i - 1]["O"])))
    return res


def case_tight_turn():
    return label(build_backbone(TIGHT))


def strand(count, phi=-139.0, psi=135.0):
    return build_backbone([(phi, psi)] * count)


def aligned_strand(count):
    """Strand with CA(0) at the origin, the CA axis along x and the C=O of
    residue 2 along +y (sheet in the xy plane)."""
    res = strand(count)
    x = unit(sub(res[-1]["CA"], res[0]["CA"]))
    co = sub(res[2]["O"], res[2]["C"])
    y = unit(sub(co, scale(x, dot(co, x))))
    z = cross(x, y)
    origin = res[0]["CA"]
    return [{k: (dot(sub(p, origin), x), dot(sub(p, origin), y), dot(sub(p, origin), z))
             for k, p in r.items()} for r in res]


def antiparallel_partner(segment, tx, ty):
    """The segment turned by 180 degrees about z through its centre, then
    shifted by (tx, ty): an antiparallel partner strand in the same sheet."""
    centre = scale(add(segment[0]["CA"], segment[-1]["CA"]), 0.5)
    moved = transform(segment, None, scale(centre, -1.0))
    return transform(moved, rotation((0.0, 0.0, 1.0), 180.0), add(centre, (tx, ty, 0.0)))


def path_point(points, s):
    for a, b in zip(points, points[1:]):
        length = dist(a, b)
        if s <= length:
            return add(a, scale(unit(sub(b, a)), s))
        s -= length
    return points[-1]


def loop_between(last_c, next_n, waypoints):
    """Backbone-like residues along the polyline last_c -> waypoints -> next_n:
    every C(k)-N(k+1) link 1.1-1.6 A, so no chain break is detected."""
    points = [last_c] + list(waypoints) + [next_n]
    total = sum(dist(a, b) for a, b in zip(points, points[1:]))
    count = max(1, round((total - 1.33) / 3.76))
    f = (total - 1.33) / (count * 3.76)
    residues = []
    for k in range(count):
        s0 = 1.33 * f + k * 3.76 * f
        n = path_point(points, s0)
        c = path_point(points, s0 + 2.43 * f)
        axis = unit(sub(c, n))
        side = unit(cross(axis, (0.3, 0.5, 0.8)))
        ca = add(scale(add(n, c), 0.5), scale(side, 0.6))
        o = add(c, scale(side, -1.23))
        residues.append({"N": n, "CA": ca, "C": c, "O": o})
    return residues


def case_bulge_sign():
    """Chain A: a 9-residue strand. Chain B: a 4-residue antiparallel partner
    of A2-A3 (P1), a loop above the sheet, then a 4-residue antiparallel
    partner of A6-A7 (P2): two ladders whose j ranges do not overlap, P2
    after P1 in sequence."""
    s1 = aligned_strand(9)
    p1 = antiparallel_partner(s1[1:5], *BULGE_P1)
    p2 = antiparallel_partner(s1[5:9], *BULGE_P2)
    lift = (0.0, 0.0, 10.0)
    loop = loop_between(p1[-1]["C"], p2[0]["N"],
                        [add(p1[-1]["C"], lift), add(p2[0]["N"], lift)])
    return label(s1, "A") + label(p1 + loop + p2, "B")


# Found by scanning with the reference above: helix angles whose i -> i+4
# energy straddles the cutoff (-0.454 kcal/mol with an N-H of 1.0 A, -0.548
# with 1.23 A); a peptide (random search) in which a search from j = i + 1
# finds a parallel bridge (2, 4) from the bonds 3 -> 4 (-0.663) and 4 -> 1
# (-0.660); the rigid shifts that pair the partner strands with chain A.
MARGINAL = (-78.0, -22.0)
TIGHT = [(-49.0, -27.0), (-157.0, 11.0), (-60.0, -27.0), (-160.0, 72.0),
         (66.0, -155.0), (-146.0, -36.0), (42.0, 72.0), (138.0, 75.0)]
BULGE_P1 = (-5.8, 4.6)
BULGE_P2 = (-5.8, 4.6)


def atoms_of(residues):
    atoms = []
    for res in residues:
        for name in ("N", "H", "CA", "C", "O"):
            if name in res:
                x, y, z = res[name]
                atoms.append((res["chain"], res["resnr"], res["resname"], name,
                              x + BOX / 2, y + BOX / 2, z + BOX / 2))
    return atoms


def read_titan(case):
    states = {}
    path = case / "dssp_mol1.dat"
    if not path.exists():
        return None, None
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        states[(parts[1], int(parts[2]))] = parts[4]
    pairs = []
    bridges = case / "dssp_mol1_bridges.csv"
    if bridges.exists():
        rows = [line for line in bridges.read_text().splitlines()
                if line.strip() and not line.startswith("#")]
        header = rows[0].split(",")
        for line in rows[1:]:
            row = dict(zip(header, line.split(",")))
            pairs.append((row["chain_i"], int(row["resnr_i"]), row["chain_j"],
                          int(row["resnr_j"]), row["orientation"]))
    return states, sorted(pairs)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    a = p.parse_args()
    binary = a.binary.resolve()
    work = a.work_root.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    checks = fx.Checks()
    # (name, residues, the old rule the case must separate from the new one)
    cases = [("helix_gap", case_helix_gap(), "breaks"),
             ("marginal_helix", case_marginal_helix(), "h"),
             ("explicit_h", case_explicit_h(), None),
             ("tight_turn", case_tight_turn(), "bridge_start"),
             ("bulge_sign", case_bulge_sign(), "bulge")]
    for name, residues, rule in cases:
        new_states, new_pairs = Dssp(residues).assign()
        old_states, old_pairs = Dssp(residues, legacy=True).assign()
        print(f"{name}: reference {new_states} (all old rules {old_states})")
        if rule is not None:
            one_states, one_pairs = Dssp(residues, legacy=[rule]).assign()
            checks.check((new_states, new_pairs) != (one_states, one_pairs),
                         f"{name}: the old '{rule}' rule alone changes the answer",
                         f"old '{rule}' rule: {one_states}")
        case = work / name
        case.mkdir()
        atoms = atoms_of(residues)
        fx.write_trajectory(case / "s.pdb", [(0.0, (BOX, BOX, BOX), atoms)])
        fx.write_topology(case / "s.top", atoms)
        fx.write_index(case / "s.ndx", [("System", range(1, len(atoms) + 1))])
        lines = fx.base_args(case, case / "s.pdb", case / "s.top", case / "s.ndx") + ["-dssp=true"]
        code, output = fx.run_titan(binary, case, lines)
        if not checks.check(code == 0, f"{name}: titan exit 0", output[-600:] if code else ""):
            continue
        states, pairs = read_titan(case)
        got = "".join(states.get((r["chain"], r["resnr"]), "?") for r in residues) if states else None
        checks.check(got == new_states, f"{name}: per-residue states",
                     f"titan {got}, expected {new_states}")
        expected_pairs = sorted((residues[i]["chain"], residues[i]["resnr"], residues[j]["chain"],
                                 residues[j]["resnr"], t) for i, j, t in new_pairs)
        checks.check(pairs == expected_pairs, f"{name}: bridge pairs",
                     f"titan {pairs}, expected {expected_pairs}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
