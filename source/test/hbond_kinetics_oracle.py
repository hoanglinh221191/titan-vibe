"""Oracle for -hbk donor hydrogens and periodic images (STR-21, MDB-01...04).

Chain A (molecule 1) and chain B (molecule 2) in a 60 A box, three identical
frames, four designed inter-chain H-bonds (D...A 2.8-2.9 A, D-H...A 180 deg):

  A1 ARG NH1 -> B1 ASP OD1 and A1 ARG NH2 -> B1 ASP OD2   (the old name rule
     took NH1/NH2 for hydrogens, so each donor's "hydrogen" was itself);
  A2 TYR OH  -> B2 GLU OE1                                  (OH as hydrogen);
  A3 LYS NZ  -> B3 GLY O through HZ2, not HZ1               (one H per donor);
  A4 SER OG  -> B4 ALA O across the x face of the box       (no minimum image).

The expected bond list is computed here from the module's rules: donors and
acceptors by residue and atom name, hydrogens = atoms whose name starts with H
within 1.2 A of the donor, D...A < 3.5 A and D-H...A >= 135 deg for any
hydrogen of the donor, minimum image in the box. Every listed inter-chain bond
must be present with 100 % occupancy and nothing else may be reported.

Runs the energy path (-g=1-2 without -nocal) so that the old binary, whose
-nocal loop never called -hbk, is tested on the same frames; a second run
with -nocal must report the same intra-chain table (it was empty before).

  python -I test/hbond_kinetics_oracle.py --binary titan.exe --work-root test/tmp/hbk_oracle
"""
import argparse
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 60.0
Z = 30.0
DONORS = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("TYR", "OH"), ("LYS", "NZ"),
          ("SER", "OG")}
ACCEPTORS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"), ("TYR", "OH"),
             ("SER", "OG")}


def add(atoms, chain, resnr, resname, name, x, y, z=Z):
    atoms.append((chain, resnr, resname, name, round(fx.wrap(x, BOX), 3), round(y, 3), round(z, 3)))


def backbone(atoms, chain, resnr, resname, x, y, sign):
    """Backbone around (x, y); sign = +1 places it to the right of x."""
    s = sign
    add(atoms, chain, resnr, resname, "N", x + s * 1.2, y - 1.4)
    add(atoms, chain, resnr, resname, "H", x + s * 1.2, y - 2.4)
    add(atoms, chain, resnr, resname, "CA", x, y)
    add(atoms, chain, resnr, resname, "HA", x, y + 0.6, Z + 0.9)
    add(atoms, chain, resnr, resname, "C", x + s * 1.3, y + 1.0)
    add(atoms, chain, resnr, resname, "O", x + s * 2.4, y + 0.6)


def build():
    a, b = [], []
    # region 0, y = 8: ARG NH1/NH2 -> ASP OD1/OD2
    y = 8.0
    backbone(a, "A", 1, "ARG", 10.0, y, -1)
    for name, x in (("CB", 11.5), ("CG", 13.0), ("CD", 14.5), ("NE", 16.0)):
        add(a, "A", 1, "ARG", name, x, y)
    add(a, "A", 1, "ARG", "HE", 15.5, y + 0.87)
    add(a, "A", 1, "ARG", "CZ", 17.3, y)
    add(a, "A", 1, "ARG", "NH1", 18.0, y + 1.15)
    add(a, "A", 1, "ARG", "HH11", 19.0, y + 1.15)
    add(a, "A", 1, "ARG", "HH12", 17.5, y + 2.02)
    add(a, "A", 1, "ARG", "NH2", 18.0, y - 1.15)
    add(a, "A", 1, "ARG", "HH21", 19.0, y - 1.15)
    add(a, "A", 1, "ARG", "HH22", 17.5, y - 2.02)
    backbone(b, "B", 1, "ASP", 25.5, y, 1)
    add(b, "B", 1, "ASP", "CB", 23.1, y)
    add(b, "B", 1, "ASP", "CG", 21.6, y)
    add(b, "B", 1, "ASP", "OD1", 20.9, y + 1.15)
    add(b, "B", 1, "ASP", "OD2", 20.9, y - 1.15)
    # region 1, y = 22: TYR OH -> GLU OE1
    y = 22.0
    backbone(a, "A", 2, "TYR", 10.0, y, -1)
    for name, x, dy in (("CB", 11.5, 0.0), ("CG", 13.0, 0.0), ("CD1", 13.7, 1.2),
                        ("HD1", 13.2, 2.1), ("CD2", 13.7, -1.2), ("HD2", 13.2, -2.1),
                        ("CE1", 15.1, 1.2), ("HE1", 15.6, 2.1), ("CE2", 15.1, -1.2),
                        ("HE2", 15.6, -2.1), ("CZ", 16.0, 0.0), ("OH", 17.4, 0.0),
                        ("HH", 18.36, 0.0)):
        add(a, "A", 2, "TYR", name, x, y + dy)
    backbone(b, "B", 2, "GLU", 25.5, y - 3.0, 1)
    add(b, "B", 2, "GLU", "CB", 22.0, y - 3.4)
    add(b, "B", 2, "GLU", "CG", 20.6, y - 2.5)
    add(b, "B", 2, "GLU", "CD", 20.9, y - 1.1)
    add(b, "B", 2, "GLU", "OE1", 20.2, y)
    add(b, "B", 2, "GLU", "OE2", 22.1, y - 0.9)
    # region 2, y = 36: LYS NZ -> GLY O through the second hydrogen
    y = 36.0
    backbone(a, "A", 3, "LYS", 10.0, y, -1)
    for name, x in (("CB", 11.5), ("CG", 13.0), ("CD", 14.5), ("CE", 16.0), ("NZ", 17.5)):
        add(a, "A", 3, "LYS", name, x, y)
    nz = (17.5, y, Z)
    dirs = [(0.333, 0.943, 0.0), (0.333, -0.471, 0.816), (0.333, -0.471, -0.816)]
    for k, d in enumerate(dirs, 1):
        add(a, "A", 3, "LYS", f"HZ{k}", nz[0] + d[0], nz[1] + d[1], nz[2] + d[2])
    o = [nz[i] + 2.9 * dirs[1][i] for i in range(3)]
    c = [o[i] + 1.23 * dirs[1][i] for i in range(3)]
    add(b, "B", 3, "GLY", "N", c[0] + 2.0, c[1] - 1.0, c[2] + 0.8)
    add(b, "B", 3, "GLY", "H", c[0] + 2.5, c[1] - 1.8, c[2] + 0.8)
    add(b, "B", 3, "GLY", "CA", c[0] + 1.4, c[1] + 0.2, c[2] + 0.6)
    add(b, "B", 3, "GLY", "HA1", c[0] + 1.8, c[1] + 1.1, c[2] + 0.2)
    add(b, "B", 3, "GLY", "C", c[0], c[1], c[2])
    add(b, "B", 3, "GLY", "O", o[0], o[1], o[2])
    # region 3, y = 50: SER OG at x = 1.2 -> ALA O at x = -1.6 (58.4)
    y = 50.0
    add(a, "A", 4, "SER", "N", 4.8, y - 1.4)
    add(a, "A", 4, "SER", "H", 4.8, y - 2.4)
    add(a, "A", 4, "SER", "CA", 3.6, y)
    add(a, "A", 4, "SER", "HA", 3.6, y + 0.6, Z + 0.9)
    add(a, "A", 4, "SER", "CB", 2.6, y)
    add(a, "A", 4, "SER", "OG", 1.2, y)
    add(a, "A", 4, "SER", "HG", 0.24, y)
    add(a, "A", 4, "SER", "C", 4.4, y + 1.2)
    add(a, "A", 4, "SER", "O", 5.6, y + 1.0)
    add(b, "B", 4, "ALA", "N", -5.0, y - 1.0)
    add(b, "B", 4, "ALA", "H", -5.0, y - 2.0)
    add(b, "B", 4, "ALA", "CA", -4.0, y)
    add(b, "B", 4, "ALA", "HA", -4.0, y + 0.6, Z + 0.9)
    add(b, "B", 4, "ALA", "CB", -4.5, y + 1.4)
    add(b, "B", 4, "ALA", "C", -2.8, y)
    add(b, "B", 4, "ALA", "O", -1.6, y)
    # intra-chain: A1 ARG NE-HE -> A5 GLU OE1 (the -nocal run reports it)
    y = 8.0
    add(a, "A", 5, "GLU", "OE1", 14.55, y + 2.52)
    add(a, "A", 5, "GLU", "OE2", 12.9, y + 3.9)
    add(a, "A", 5, "GLU", "CD", 14.0, y + 3.6)
    add(a, "A", 5, "GLU", "CG", 14.6, y + 4.9)
    add(a, "A", 5, "GLU", "CB", 14.6, y + 6.4)
    backbone(a, "A", 5, "GLU", 14.6, y + 7.9, 1)
    return a, b


def mic(p, q):
    d = [q[k] - p[k] for k in range(3)]
    return [v - BOX * round(v / BOX) for v in d]


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def expected_bonds(atoms, scopes):
    """(scope, donor chain, donor resnr, donor atom, acceptor chain, acceptor resnr, acceptor atom)."""
    def is_donor(r, n):
        return n == "N" or (r, n) in DONORS

    def is_acceptor(r, n):
        return n == "O" or (r, n) in ACCEPTORS

    pos = {i: (a[4], a[5], a[6]) for i, a in enumerate(atoms)}
    bonds = set()
    for i, d in enumerate(atoms):
        if not is_donor(d[2], d[3]):
            continue
        hyd = [j for j, h in enumerate(atoms) if j != i and h[0] == d[0] and h[1] == d[1]
               and h[3].lstrip("0123456789").startswith("H") and norm(mic(pos[i], pos[j])) < 1.2]
        for k, acc in enumerate(atoms):
            if k == i or not is_acceptor(acc[2], acc[3]):
                continue
            same_chain = acc[0] == d[0]
            scope = "intra" if same_chain else "inter"
            if scope not in scopes or (same_chain and acc[1] == d[1]):
                continue
            if not norm(mic(pos[i], pos[k])) < 3.5:
                continue
            for j in hyd:
                hd = mic(pos[j], pos[i])
                ha = mic(pos[j], pos[k])
                cos = sum(x * y for x, y in zip(hd, ha)) / (norm(hd) * norm(ha))
                if cos <= math.cos(math.radians(135.0)):
                    bonds.add((scope, d[0], d[1], d[3], acc[0], acc[1], acc[3]))
                    break
    return bonds


def read_pairs(case):
    files = [f for f in sorted(case.glob("hbond_kinetics_hbond_pairs*.csv")) if "average" not in f.name]
    if not files:
        return None
    out = {}
    for row in fx.read_csv(files[0]):
        scope = "inter" if row["Scope"] == "inter" else "intra"
        key = (scope, row["DonorChain"], int(row["DonorResID"]), row["DonorAtom"],
               row["AcceptorChain"], int(row["AcceptorResID"]), row["AcceptorAtom"])
        out[key] = float(row["Occupancy_%"])
    return out


def frames_seen(case):
    files = [f for f in sorted(case.glob("hbond_kinetics_summary*.dat")) if "average" not in f.name]
    return fx.read_key_values(files[0]).get("frames_seen") if files else None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    args = p.parse_args()
    binary = args.binary.resolve()
    work = args.work_root.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    chain_a, chain_b = build()
    atoms = chain_a + chain_b
    fr = [(t, (BOX, BOX, BOX), atoms) for t in (0.0, 10.0, 20.0)]
    fx.write_trajectory(work / "pair.pdb", fr)
    fx.write_topology(work / "pair.top", atoms)
    n_a = len(chain_a)
    fx.write_index(work / "pair.ndx", [("System", range(1, len(atoms) + 1)),
                                       ("chA", range(1, n_a + 1)),
                                       ("chB", range(n_a + 1, len(atoms) + 1))])
    designed = {("inter", "A", 1, "NH1", "B", 1, "OD1"), ("inter", "A", 1, "NH2", "B", 1, "OD2"),
                ("inter", "A", 2, "OH", "B", 2, "OE1"), ("inter", "A", 3, "NZ", "B", 3, "O"),
                ("inter", "A", 4, "OG", "B", 4, "O")}
    checks = fx.Checks()
    want = expected_bonds(atoms, {"inter", "intra"})
    checks.check(designed <= want, "fixture: the designed bonds satisfy the criteria",
                 str(sorted(designed - want)))

    base = [a for a in fx.base_args(work, work / "pair.pdb", work / "pair.top", work / "pair.ndx",
                                    groups="1-2") if a != "-nocal=true"]
    hbk = ["-hbk=true", "-hbk-domain=both"]
    code, output = fx.run_titan(binary, work / "energy", base + hbk)
    if checks.check(code == 0, "energy path: titan exit 0", f"exit {code}: {output[-800:]}" if code else ""):
        got = read_pairs(work / "energy") or {}
        for bond in sorted(designed):
            label = f"{bond[1]}{bond[2]} {bond[3]} -> {bond[4]}{bond[5]} {bond[6]}"
            checks.check(got.get(bond) == 100.0, f"energy path: {label} at 100 %",
                         f"occupancy {got.get(bond)}")
        checks.check(set(got) == want, "energy path: reported bonds = computed bonds",
                     f"extra {sorted(set(got) - want)}, missing {sorted(want - set(got))}")
    nocal = base + ["-nocal=true"] + hbk
    code, output = fx.run_titan(binary, work / "nocal", nocal)
    if checks.check(code == 0, "-nocal: titan exit 0", f"exit {code}: {output[-800:]}" if code else ""):
        checks.check(frames_seen(work / "nocal") == "3", "-nocal: -hbk saw the three frames",
                     f"frames_seen {frames_seen(work / 'nocal')}")
        got = read_pairs(work / "nocal") or {}
        intra = {b for b in want if b[0] == "intra"}
        checks.check(("intra", "A", 1, "NE", "A", 5, "OE1") in intra,
                     "fixture: A1 ARG NE -> A5 GLU OE1 is an intra-chain bond")
        checks.check(set(got) == intra, "-nocal: intra-chain bonds = computed (inter is off in -nocal)",
                     f"got {sorted(got)}, want {sorted(intra)}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
