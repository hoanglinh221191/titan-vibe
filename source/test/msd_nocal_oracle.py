"""STR-28 oracle: -msd in a -nocal run must follow the molecule's own motion.

Rigid bodies drift at constant velocity v = (0.30, 0.20, 0.10) A/ps through a
40 A box, 200 frames, dt = 1 ps, run with -nocal -g=0-0 -msd-mode=com; the
centre-of-mass MSD is then exactly |v|^2 t^2 (ballistic, no noise).

  small  three atoms, whole in every frame. Before STR-28 the automatic
         -pbci 0 of -nocal moved molecule 1's centre to the box centre every
         frame (align_to_box), so the MSD was 0.
  rod    30 atoms over 29 A (longer than half the box), each atom folded into
         the box on its own as in a raw trajectory. Before STR-28 the body was
         rebuilt by imaging every atom next to atom 0, which breaks a body
         larger than half the box and makes its centre jump.

  python -I test/msd_nocal_oracle.py --binary titan.exe --work-root test/tmp/msd_nocal_oracle
"""
import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 40.0
FRAMES = 200
VELOCITY = (0.30, 0.20, 0.10)
SMALL = [(0.0, 0.0, 0.0), (1.4, 0.0, 0.0), (0.0, 1.4, 0.0)]
ROD = [(1.0 * k, 0.3 * (k % 2), 0.0) for k in range(30)]
LAGS = [1, 5, 20, 60]


def run_case(binary, work, name, body, fold_atoms):
    """Writes and runs one case; returns failure strings."""
    case = work / name
    case.mkdir()
    frames = []
    for i in range(FRAMES):
        centre = [10.0 + VELOCITY[k] * i for k in range(3)]
        shift = [fx.wrap(centre[k], BOX) - centre[k] for k in range(3)]
        atoms = []
        for n, b in enumerate(body):
            pos = [centre[k] + b[k] for k in range(3)]
            if fold_atoms:
                pos = [fx.wrap(p, BOX) for p in pos]
            else:
                pos = [pos[k] + shift[k] for k in range(3)]
            atoms.append(("A", 1, "MOL", f"C{n + 1}", pos[0], pos[1], pos[2]))
        frames.append((float(i), (BOX, BOX, BOX), atoms))
    fx.write_trajectory(case / "drift.pdb", frames)
    fx.write_topology(case / "drift.top", frames[0][2])
    fx.write_index(case / "drift.ndx", [("System", list(range(1, len(body) + 1)))])
    lines = fx.base_args(case, case / "drift.pdb", case / "drift.top", case / "drift.ndx")
    lines += ["-msd=true", "-msd-mode=com"]
    code, _ = fx.run_titan(binary, case / "run", lines)

    failures = []
    found = sorted((case / "run").glob("msd_timeseries*.dat"))  # no _<n> suffix for PDB input
    if code != 0 or not found:
        return [f"{name}: titan exit {code} or no msd_timeseries*.dat (see {case / 'run'})"]
    rows = {}
    for line in found[0].read_text().splitlines()[1:]:
        cols = line.split()
        rows[int(cols[0])] = float(cols[2])
    v2 = sum(v * v for v in VELOCITY)
    for lag in LAGS:
        want = v2 * lag * lag
        got = rows.get(lag)
        ok = got is not None and abs(got - want) <= 1e-3 * want + 1e-4
        print(f"  {name:5s} lag {lag:3d} ps: MSD {got if got is not None else 'missing'} A^2, "
              f"expected {want:.4f}  {'ok' if ok else 'FAIL'}")
        if not ok:
            failures.append(f"{name} lag {lag}: MSD {got}, expected {want:.4f}")
    return failures


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    a = p.parse_args()
    work = a.work_root.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    binary = a.binary.resolve()
    failures = run_case(binary, work, "small", SMALL, fold_atoms=False)
    failures += run_case(binary, work, "rod", ROD, fold_atoms=True)
    if failures:
        print("FAIL")
        for f in failures:
            print("  -", f)
        return 1
    print("PASS: -nocal -msd follows the molecule's own drift")
    return 0


if __name__ == "__main__":
    sys.exit(main())
