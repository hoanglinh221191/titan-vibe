"""Oracle for what the membrane order-parameter and "curvature" outputs are
(STR-11, MDB-19).

martini_bump  the 128-lipid MARTINI bilayer of membrane_pbc_oracle.py
              (85 x 85 x 70 A box, centre z = 35), with the 32 lipids of one
              2 x 2 grid cell (x, y < 40 A) raised by 2 A, -membrane-normal=z,
              -membrane-grid=2,2, -membrane-curvature=grid. The map must be
              the height deviation (mean head height in the cell - leaflet
              mean height)/sqrt(cell area), averaged over the two leaflets:
              +1.5/42.5 in the raised cell, -0.5/42.5 elsewhere, in a column
              named height_dev_over_cell_edge (the old header said
              curvature_proxy_Ainv, a curvature in 1/A).
              Vertical tails: P2 of every bead bond = 1, column mean_P2_bond.
aa_tilted     a bilayer of lipids with all-atom names (head P, tails
              C22-C23-C24 and C32-C33-C34, -membrane-*=atomname:...) whose
              C-C bonds are tilted 30 degrees from the normal: the output is
              S_CC = P2(cos 30) = 0.625 for every bond, labelled P2_bond, and
              the summary note says S_CC, not S_CD (the old output called it
              mean_S with a MARTINI-only note).

  python -I test/membrane_labels_oracle.py --binary titan.exe --work-root D:/softs/titan_work/str/membrane_labels_oracle
"""
import argparse
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402
import membrane_pbc_oracle as mp  # noqa: E402

ZC = 35.0
BUMP = 2.0


def martini_frame():
    atoms = []
    resnr = 0
    for side in (1, -1):
        for ix in range(mp.N_SIDE):
            for iy in range(mp.N_SIDE):
                resnr += 1
                x, y = 5.0 + 10.0 * ix, 5.0 + 10.0 * iy
                lift = BUMP if (x < 40.0 and y < 40.0) else 0.0
                for name, bx, by, bz in mp.lipid_beads(x, y, ZC + lift, side):
                    atoms.append(("L", resnr, "POPC", name, bx, by, bz))
    return atoms


def aa_frame():
    atoms = []
    resnr = 0
    step = 1.53
    tilt = math.radians(30.0)
    for side in (1, -1):
        for ix in range(mp.N_SIDE):
            for iy in range(mp.N_SIDE):
                resnr += 1
                x, y = 5.0 + 10.0 * ix, 5.0 + 10.0 * iy
                zh = ZC + side * 20.0
                atoms.append(("L", resnr, "POPC", "P", x, y, zh))
                for tail, dy in (("2", -1.0), ("3", 1.0)):
                    for k in range(3):
                        d = step * (k + 1)
                        atoms.append(("L", resnr, "POPC", f"C{tail}{k + 2}",
                                      x + d * math.sin(tilt), y + dy,
                                      zh - side * (2.0 + d * math.cos(tilt))))
    return atoms


def rows(case, stem):
    files = [f for f in sorted(case.glob(stem + "*.csv")) if "average" not in f.name]
    return fx.read_csv(files[0]) if files else None


def summary_text(case):
    files = [f for f in sorted(case.glob("membrane_summary*.dat")) if "average" not in f.name]
    return files[0].read_text() if files else ""


def run_case(binary, work, name, atoms, extra):
    frames = [(0.0, (mp.LX, mp.LY, mp.LZ), atoms), (10.0, (mp.LX, mp.LY, mp.LZ), atoms)]
    fx.write_trajectory(work / f"{name}.pdb", frames)
    fx.write_topology(work / f"{name}.top", atoms)
    fx.write_index(work / f"{name}.ndx", [("lipid", range(1, len(atoms) + 1))])
    lines = fx.base_args(work, work / f"{name}.pdb", work / f"{name}.top", work / f"{name}.ndx")
    lines += ["-membrane=true", "-membrane-group=auto", "-membrane-normal=z",
              "-membrane-grid=2,2"] + extra
    return fx.run_titan(binary, work / name, lines)


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

    code, output = run_case(binary, work, "martini_bump", martini_frame(), ["-membrane-curvature=grid"])
    if checks.check(code == 0, "martini_bump: titan exit 0", output[-800:] if code else ""):
        case = work / "martini_bump"
        grid = rows(case, "membrane_curvature_map") or []
        column = "height_dev_over_cell_edge"
        has = bool(grid) and column in grid[0]
        checks.check(has, f"martini_bump: curvature map column {column}",
                     f"header {list(grid[0]) if grid else None}")
        if has:
            # cells: u = (0,-1,0), v = (1,0,0) for the normal +z, so the
            # raised lipids (x, y < 40 A) are in cell (grid_x 1, grid_y 0)
            edge = mp.LX / 2.0
            expected = {(ix, iy): (1.5 if (ix, iy) == (1, 0) else -0.5) / edge
                        for ix in range(2) for iy in range(2)}
            got = {(int(r["grid_x"]), int(r["grid_y"])): float(r[column]) for r in grid}
            checks.check(all(abs(got.get(k, 1e9) - v) < 1e-6 for k, v in expected.items()),
                         "martini_bump: height deviation (h_cell - h_leaflet)/sqrt(A_cell)",
                         f"titan {got}, expected {expected}")
        order = rows(case, "membrane_order_parameters") or []
        checks.check(bool(order) and "mean_P2_bond" in order[0] and
                     all(abs(float(r["mean_P2_bond"]) - 1.0) < 1e-6 for r in order),
                     "martini_bump: mean_P2_bond = 1 for vertical bead bonds",
                     f"header {list(order[0]) if order else None}")
        text = summary_text(case)
        checks.check("not a mean or Gaussian curvature" in text,
                     "martini_bump: summary says the map is not a curvature")

    aa_extra = ["-membrane-lipids=resname:POPC", "-membrane-head=atomname:P",
                "-membrane-tail=atomname:C22,C23,C24;C32,C33,C34"]
    code, output = run_case(binary, work, "aa_tilted", aa_frame(), aa_extra)
    if checks.check(code == 0, "aa_tilted: titan exit 0", output[-800:] if code else ""):
        case = work / "aa_tilted"
        order = rows(case, "membrane_order_parameters") or []
        # P2 of the C-C vectors as written to the PDB (3 decimals), averaged
        # over the lipids: 0.625 up to the rounding of the coordinates
        written = [(a[1], a[3], tuple(float(f"{v:8.3f}") for v in a[4:7])) for a in aa_frame()]
        pos = {(r, n): xyz for r, n, xyz in written}
        expected = {}
        for tail in ("2", "3"):
            for bond in (1, 2):
                values = []
                for resnr in {r for r, _, _ in written}:
                    a, b = pos[(resnr, f"C{tail}{bond + 1}")], pos[(resnr, f"C{tail}{bond + 2}")]
                    d = [b[k] - a[k] for k in range(3)]
                    values.append(0.5 * (3.0 * d[2] ** 2 / sum(v * v for v in d) - 1.0))
                expected[(f"tail{int(tail) - 1}", f"bond{bond}")] = sum(values) / len(values)
        p2 = 0.5 * (3.0 * math.cos(math.radians(30.0)) ** 2 - 1.0)
        ok = bool(order) and "mean_P2_bond" in order[0] and len(order) == 4 and all(
            abs(float(r["mean_P2_bond"]) - expected[(r["tail"], r["bond"])]) < 2e-6 and
            abs(float(r["mean_P2_bond"]) - p2) < 1e-4 for r in order)
        checks.check(ok, f"aa_tilted: S_CC = P2(cos 30) = {p2:.4f} for the 4 C-C bonds",
                     f"{order}")
        text = summary_text(case)
        checks.check("S_CC" in text and "not the NMR S_CD" in text,
                     "aa_tilted: summary names the quantity S_CC, not S_CD")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
