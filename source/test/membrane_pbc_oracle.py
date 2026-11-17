"""Oracle for a bilayer that straddles the box along its normal (STR-22).

128 MARTINI-like POPC (8 x 8 per leaflet, 10 A apart) in an 80 x 80 x 70 A box,
head beads (PO4) 20 A above and below the bilayer centre, tails pointing to
the centre and ending 4 A from it. Six frames, every bead wrapped into the box:
the centre is at z = 35, 0, 68, 3, 66, 35, so in frames 1-4 the bilayer is
split by the z face and the lower leaflet's heads sit at the top of the box.
The bilayer (40 A) is thicker than the water layer (30 A): the head beads
alone cannot tell the bilayer centre from the water centre, the jumps of the
centre between frames are larger than any previous-frame scheme tolerates.

For -membrane-normal=z and for the automatic normal, every frame must give:
thickness 40 A, every lipid in the leaflet of frame 0 (no flip-flop event),
64 lipids per leaflet, mean APL = box area / 64, tail penetration 0 and the tail
order parameter S = 1 for every bond. The old code used the raw head
coordinates: centre in the water, leaflets swapped (128 flip-flop events),
thickness 60 A, and for the automatic normal an in-plane "normal".

  python -I test/membrane_pbc_oracle.py --binary titan.exe --work-root test/tmp/membrane_pbc_oracle
"""
import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

N_SIDE = 8  # lipids per side and leaflet; 80 A patch: in-plane head spread > normal spread
# a 5 A wider box than the lattice: no two heads are exactly half a box apart
# in the plane (a tie of the minimum image, which real coordinates never hit)
LX, LY, LZ = 10.0 * N_SIDE + 5.0, 10.0 * N_SIDE + 5.0, 70.0
APL = LX * LY / (N_SIDE * N_SIDE)
CENTRES = [35.0, 0.0, 68.0, 3.0, 66.0, 35.0]
HALF = 20.0


def lipid_beads(x, y, zc, side):
    """side +1: upper leaflet (head above the centre), -1: lower."""
    zh = zc + side * HALF
    beads = [("NC3", x, y, zh + side * 2.0), ("PO4", x, y, zh),
             ("GL1", x - 0.5, y, zh - side * 2.0), ("GL2", x + 0.5, y, zh - side * 2.0)]
    for k in range(4):
        z = zh - side * (4.0 + 4.0 * k)
        beads.append((f"C{k + 1}A", x - 0.5, y, z))
        beads.append((f"C{k + 1}B", x + 0.5, y, z))
    return beads


def frame_atoms(zc):
    atoms = []
    resnr = 0
    for side in (1, -1):
        for ix in range(N_SIDE):
            for iy in range(N_SIDE):
                resnr += 1
                for name, x, y, z in lipid_beads(5.0 + 10.0 * ix, 5.0 + 10.0 * iy, zc, side):
                    atoms.append(("L", resnr, "POPC", name, fx.wrap(x, LX), fx.wrap(y, LY),
                                  fx.wrap(z, LZ)))
    return atoms


def rows(case, stem):
    files = [f for f in sorted(case.glob(stem + "*.csv")) if "average" not in f.name]
    return fx.read_csv(files[0]) if files else None


def interdigitation(case):
    files = [f for f in sorted(case.glob("membrane_interdigitation*.dat")) if "average" not in f.name]
    if not files:
        return None
    lines = files[0].read_text().split("\n")
    return [[float(v) for v in l.split()] for l in lines[1:] if l.strip()]


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

    frames = [(10.0 * i, (LX, LY, LZ), frame_atoms(zc)) for i, zc in enumerate(CENTRES)]
    atoms0 = frames[0][2]
    fx.write_trajectory(work / "bilayer.pdb", frames)
    fx.write_topology(work / "bilayer.top", atoms0)
    fx.write_index(work / "bilayer.ndx", [("lipid", range(1, len(atoms0) + 1))])
    per_leaflet = N_SIDE * N_SIDE
    upper_ids = {f"L:{r}:POPC" for r in range(1, per_leaflet + 1)}

    checks = fx.Checks()
    base = fx.base_args(work, work / "bilayer.pdb", work / "bilayer.top", work / "bilayer.ndx")
    for name, normal in (("normal_z", "z"), ("normal_auto", "auto")):
        lines = base + ["-membrane=true", "-membrane-group=auto", f"-membrane-normal={normal}",
                        "-membrane-grid=2,2", "-membrane-flipflop-persistence=3"]
        case = work / name
        code, output = fx.run_titan(binary, case, lines)
        if not checks.check(code == 0, f"{name}: titan exit 0", f"exit {code}: {output[-800:]}" if code else ""):
            continue
        thick = rows(case, "membrane_thickness_timeseries") or []
        values = [float(r["mean_thickness_A"]) for r in thick]
        checks.check(len(values) == len(CENTRES) and all(abs(v - 40.0) < 1e-6 for v in values),
                     f"{name}: thickness 40 A in every frame", f"{values}")
        leaf = rows(case, "membrane_leaflets") or []
        first = {r["lipid_id"]: r["leaflet"] for r in leaf if r["frame"] == "0"}
        upper_label = {first[i] for i in upper_ids if i in first}
        consistent = len(first) == 2 * per_leaflet and len(upper_label) == 1 and all(
            r["leaflet"] == first.get(r["lipid_id"]) for r in leaf)
        changed = sorted({(r["frame"], r["lipid_id"]) for r in leaf
                          if r["leaflet"] != first.get(r["lipid_id"])})[:4]
        checks.check(consistent, f"{name}: every lipid keeps its frame-0 leaflet",
                     f"{len(leaf)} rows, first changes {changed}")
        heights = [abs(float(r["signed_distance_A"])) for r in leaf]
        checks.check(heights and all(abs(h - 20.0) < 1e-6 for h in heights),
                     f"{name}: |signed distance| of every head = 20 A",
                     f"range {min(heights) if heights else None}..{max(heights) if heights else None}")
        flips = rows(case, "membrane_flipflop_events") or []
        checks.check(len(flips) == 0, f"{name}: no flip-flop event", f"{len(flips)} events")
        apl = rows(case, "membrane_apl_timeseries") or []
        checks.check(len(apl) == len(CENTRES) and all(
            abs(float(r["mean_apl_A2"]) - APL) < 1e-3 and r["upper_lipids"] == str(per_leaflet)
            and r["lower_lipids"] == str(per_leaflet) for r in apl),
            f"{name}: APL {APL:.4f} A^2 and {per_leaflet} + {per_leaflet} lipids in every frame",
            f"{[(r['mean_apl_A2'], r['upper_lipids'], r['lower_lipids']) for r in apl]}")
        inter = interdigitation(case) or []
        checks.check(len(inter) == len(CENTRES) and all(abs(r[3]) < 1e-6 for r in inter),
                     f"{name}: no tail penetration", f"{[r[3] for r in inter]}")
        order = rows(case, "membrane_order_parameters") or []
        column = "mean_P2_bond" if order and "mean_P2_bond" in order[0] else "mean_S"
        checks.check(order and all(abs(float(r[column]) - 1.0) < 1e-6 for r in order),
                     f"{name}: order parameter P2 = 1 for every tail bond",
                     f"{[r[column] for r in order][:6]}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
