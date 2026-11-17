"""Oracle for the LLPS degree distribution P(k) (STR-04).

Six one-bead chains A-F (the LLPS macro nodes) in a 150 A box, contact cutoff
7 A, four frames with a designed contact graph:

  frame 0  path A-B-C (5 A spacing), D E F apart   degrees 1 2 1 0 0 0
  frame 1  pairs A-B, C-D, E-F                     degrees 1 1 1 1 1 1
  frame 2  star around B (A, C, D at 5 A from B,
           7.07 or 10 A from each other), E F apart degrees 1 3 1 1 0 0
  frame 3  all apart                               degrees 0 0 0 0 0 0

P_f(k) is the fraction of nodes with degree k in frame f. llps_degree_dist.dat
must hold, for every k, the mean of P_f(k) over the frames and its standard
deviation over the frames (population, 1/N), computed here with fractions.
The old code accumulated (1/n)^2 once per NODE instead of P_f(k)^2 once per
frame, which gives sqrt(sum_f c_f(k)/n^2 / N - mean^2): 0 for k = 0 and 1
here, where the exact value is 0.360844.

  python -I test/llps_degree_oracle.py --binary titan.exe --work-root D:/softs/titan_work/str/llps_degree_oracle
"""
import argparse
import math
import shutil
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 150.0
CHAINS = "ABCDEF"


def frame_positions():
    """Bead positions (A) of the six chains per frame."""
    far = [(20.0 + 20.0 * k, 120.0, 120.0) for k in range(6)]
    f0 = [(50.0, 50.0, 50.0), (55.0, 50.0, 50.0), (60.0, 50.0, 50.0), far[3], far[4], far[5]]
    f1 = [(30.0, 30.0, 30.0), (35.0, 30.0, 30.0), (60.0, 30.0, 30.0), (65.0, 30.0, 30.0),
          (90.0, 30.0, 30.0), (95.0, 30.0, 30.0)]
    f2 = [(55.0, 50.0, 50.0), (50.0, 50.0, 50.0), (45.0, 50.0, 50.0), (50.0, 55.0, 50.0),
          far[4], far[5]]
    f3 = far
    return [f0, f1, f2, f3]


def degrees(positions, cutoff=7.0):
    deg = [0] * len(positions)
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            if math.dist(positions[i], positions[j]) < cutoff:
                deg[i] += 1
                deg[j] += 1
    return deg


def expected_distribution(all_degrees):
    n = len(all_degrees[0])
    kmax = max(max(d) for d in all_degrees)
    per_frame = [[Fraction(d.count(k), n) for k in range(kmax + 1)] for d in all_degrees]
    frames = len(per_frame)
    mean = [sum(p[k] for p in per_frame) / frames for k in range(kmax + 1)]
    var = [sum((p[k] - mean[k]) ** 2 for p in per_frame) / frames for k in range(kmax + 1)]
    old = []
    for k in range(kmax + 1):
        sq = sum(Fraction(d.count(k), n * n) for d in all_degrees)
        old.append(math.sqrt(max(0.0, float(sq / frames - mean[k] ** 2))))
    return [float(m) for m in mean], [math.sqrt(float(v)) for v in var], old


def read_distribution(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        k, mean, std = line.split()
        rows.append((int(k), float(mean), float(std)))
    return rows


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

    positions = frame_positions()
    all_degrees = [degrees(pos) for pos in positions]
    checks.check(all_degrees == [[1, 2, 1, 0, 0, 0], [1, 1, 1, 1, 1, 1], [1, 3, 1, 1, 0, 0],
                                 [0, 0, 0, 0, 0, 0]], "designed degree sequences",
                 str(all_degrees))
    mean, std, old = expected_distribution(all_degrees)
    print("expected mean", [round(v, 6) for v in mean])
    print("expected std ", [round(v, 6) for v in std])
    print("old estimator", [round(v, 6) for v in old])
    checks.check(max(abs(x - y) for x, y in zip(std, old)) > 0.1,
                 "the old estimator differs from the frame standard deviation on this input")

    frames = []
    for f, pos in enumerate(positions):
        atoms = [(CHAINS[i], 1, "ALA", "CA", *pos[i]) for i in range(6)]
        # group 2: one far bead (the LLPS core runs in the energy frame loop)
        atoms.append(("K", 1, "LIG", "X1", 140.0, 140.0, 140.0))
        frames.append((1000.0 * f, (BOX, BOX, BOX), atoms))
    atoms0 = frames[0][2]
    fx.write_trajectory(work / "llps.pdb", frames)
    fx.write_topology(work / "llps.top", atoms0)
    n = len(atoms0)
    fx.write_index(work / "llps.ndx", [("System", range(1, n + 1)), ("Chains", range(1, n)),
                                       ("Far", [n])])
    lines = [x for x in fx.base_args(work, work / "llps.pdb", work / "llps.top", work / "llps.ndx",
                                     groups="1-2") if not x.startswith("-nocal")]
    lines += ["-mol-1-type=protein", "-mol-2-type=lig", "-llps=true", "-llps-dbscan=false",
              "-llps-density=off", "-llps-rcut=0.7"]
    code, output = fx.run_titan(binary, work / "case", lines, timeout=600)
    if checks.check(code == 0, "titan exit 0", output[-800:] if code else ""):
        path = work / "case" / "llps_degree_dist.dat"
        got = read_distribution(path) if path.exists() else []
        checks.check(len(got) == len(mean), "one row per degree 0..kmax",
                     f"{len(got)} rows, expected {len(mean)}")
        checks.check(all(abs(g[1] - e) <= 2e-6 for g, e in zip(got, mean)), "P(k) mean",
                     f"titan {[g[1] for g in got]}")
        checks.check(all(abs(g[2] - e) <= 2e-6 for g, e in zip(got, std)),
                     "P(k) standard deviation over frames",
                     f"titan {[g[2] for g in got]}, expected {[round(v, 6) for v in std]}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
