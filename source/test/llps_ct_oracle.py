"""Oracle for the LLPS contact-persistence function C(t) (STR-03).

Ten one-bead chains in a 150 A box, 1600 frames:
  - chains A-F form a line; their spacing is 10 A (no contact at the default
    7 A cutoff) for frames 0-499 and 6 A (in contact) from frame 500 on, so
    the giant component reaches 6/10 and the system percolates at frame 500;
  - chains G-J touch A-D (5.5 A) or sit far away. Their contacts follow a
    two-state Markov chain from frame 500 on (k_on 0.05, k_off 0.10 per
    frame: occupancy 1/3, lambda = 0.85); G touches A in every frame before
    percolation as well;
  - the frame times are 0, 500, 2000, 3000, ... ps: the median step is
    1 ns, the first one 0.5 ns.

C(t) is computed here as titan defines it: over the post-percolation frames
P, the persistent pattern M marks the pairs in contact in at least half of P,
d(t) = A(t) xor M, C(k) = sum_t |d(t) and d(t+k)| / sum_t |d(t)| over all
origins t of P, lags k < min(maxlag, |P|/2), lag time k * median(dt); tau is
the 1/e crossing. The old code counted M over every frame (G-A became
persistent through its pre-percolation contacts), used only the last
2*500 frames of the ring as origins and took dt from the first two frames.
The four Markov pairs also give C(k) close to 1/3 + 2/3 * 0.85^k.

  python -I test/llps_ct_oracle.py --binary titan.exe --work-root D:/softs/titan_work/str/llps_ct_oracle
"""
import argparse
import math
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 150.0
FRAMES = 1600
PERCOLATION = 500
MAXLAG = 50
CHAINS = "ABCDEFGHIJ"
K_ON, K_OFF = 0.05, 0.10
OFFSETS = [(0.0, 5.5, 0.0), (0.0, 0.0, 5.5), (0.0, -5.5, 0.0), (0.0, 0.0, -5.5)]


def times_ps():
    return [0.0, 500.0] + [1000.0 * k for k in range(2, FRAMES)]


def contact_series():
    """on[c][f] for the flicker pair (G+c, A+c), c = 0..3."""
    rng = random.Random(20261022)
    series = []
    for c in range(4):
        state = c == 0
        values = []
        for f in range(FRAMES):
            if f < PERCOLATION:
                state = c == 0
            elif f == PERCOLATION:
                state = rng.random() < K_ON / (K_ON + K_OFF)
            else:
                state = (rng.random() >= K_OFF) if state else (rng.random() < K_ON)
            values.append(state)
        series.append(values)
    return series


def frames_of(series):
    out = []
    for f, t in enumerate(times_ps()):
        spacing = 10.0 if f < PERCOLATION else 6.0
        line = [(20.0 + spacing * k, 50.0, 50.0) for k in range(6)]
        atoms = []
        for k, (x, y, z) in enumerate(line):
            atoms.append((CHAINS[k], 1, "ALA", "CA", x, y, z))
        for c in range(4):
            x, y, z = line[c]
            if series[c][f]:
                dx, dy, dz = OFFSETS[c]
                pos = (x + dx, y + dy, z + dz)
            else:
                pos = (x, 100.0 + 10.0 * c, 100.0)
            atoms.append((CHAINS[6 + c], 1, "ALA", "CA", *pos))
        # group 2: one far bead (the LLPS core runs in the energy frame loop,
        # which needs two groups)
        atoms.append(("K", 1, "LIG", "X1", 140.0, 140.0, 140.0))
        out.append((t, (BOX, BOX, BOX), atoms))
    return out


def adjacency(frames):
    """Set of chain pairs within 7 A per frame (one bead per chain)."""
    result = []
    for _, _, atoms in frames:
        pos = [a[4:7] for a in atoms if a[2] == "ALA"]
        pairs = set()
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                if math.dist(pos[i], pos[j]) < 7.0:
                    pairs.add((i, j))
        result.append(pairs)
    return result


def ct(adj, first, origins_from, mask_frames, dt, maxlag):
    window = adj[first:]
    counts = {}
    for frame in mask_frames:
        for p in frame:
            counts[p] = counts.get(p, 0) + 1
    mask = {p for p, n in counts.items() if 2 * n >= len(window)}
    dev = [frame ^ mask for frame in adj[origins_from:]]
    n_lags = min(maxlag, len(dev) // 2)
    values = []
    for k in range(n_lags):
        num = sum(len(dev[t] & dev[t + k]) for t in range(len(dev) - k))
        den = sum(len(dev[t]) for t in range(len(dev) - k))
        values.append(num / den if den else 0.0)
    values = [v / values[0] for v in values] if values and values[0] > 1e-12 else values
    lags = [k * dt for k in range(n_lags)]
    return lags, values


def fit_tau(values, lags):
    inv_e = 1.0 / math.e
    for i in range(1, len(values)):
        if values[i] <= inv_e and values[i - 1] > values[i]:
            frac = (values[i - 1] - inv_e) / (values[i - 1] - values[i])
            return lags[i - 1] + frac * (lags[i] - lags[i - 1])
    return sum(0.5 * (values[i - 1] + values[i]) * (lags[i] - lags[i - 1]) for i in range(1, len(values)))


def read_ct(case):
    path = case / "llps_autocorr.dat"
    if not path.exists():
        return None, None
    lags, values = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        a, b = line.split()
        lags.append(float(a))
        values.append(float(b))
    return lags, values


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

    frames = frames_of(contact_series())
    adj = adjacency(frames)
    times = [t / 1000.0 for t in times_ps()]
    steps = sorted(b - a for a, b in zip(times, times[1:]) if b > a)
    dt = steps[len(steps) // 2] if len(steps) % 2 else 0.5 * (steps[len(steps) // 2 - 1] + steps[len(steps) // 2])
    lags, expected = ct(adj, PERCOLATION, PERCOLATION, adj[PERCOLATION:], dt, MAXLAG)
    tau = fit_tau(expected, lags)
    # the old estimator: mask over every frame against |P|/2, origins = the
    # last 1000 frames, dt from the first two frames
    old_counts_window = adj
    old_lags, old = ct(adj, PERCOLATION, FRAMES - 1000, old_counts_window, times[1] - times[0], MAXLAG)
    markov = [1.0 / 3.0 + 2.0 / 3.0 * 0.85 ** k for k in range(MAXLAG)]
    print(f"expected C(1..5) {[round(v, 4) for v in expected[1:6]]}, tau {tau:.4f} ns; "
          f"old C(1..5) {[round(v, 4) for v in old[1:6]]}, tau {fit_tau(old, old_lags):.4f} ns")
    checks.check(max(abs(x - y) for x, y in zip(expected, old)) > 0.05,
                 "the old estimator differs from the definition on this input")
    checks.check(max(abs(x - y) for x, y in zip(expected, markov)) < 0.1,
                 "C(k) of the input follows the Markov value 1/3 + 2/3 0.85^k within 0.1",
                 f"max deviation {max(abs(x - y) for x, y in zip(expected, markov)):.4f}")

    atoms0 = frames[0][2]
    fx.write_trajectory(work / "llps.pdb", frames)
    fx.write_topology(work / "llps.top", atoms0)
    n = len(atoms0)
    fx.write_index(work / "llps.ndx", [("System", range(1, n + 1)), ("Chains", range(1, n)),
                                       ("Far", [n])])
    lines = [a for a in fx.base_args(work, work / "llps.pdb", work / "llps.top", work / "llps.ndx",
                                     groups="1-2") if not a.startswith("-nocal")]
    lines += ["-mol-1-type=protein", "-mol-2-type=lig", "-llps=true", "-llps-dbscan=false",
              "-llps-density=off", "-llps-rcut=0.7", f"-llps-maxlag={MAXLAG}"]
    code, output = fx.run_titan(binary, work / "case", lines, timeout=1200)
    if checks.check(code == 0, "titan exit 0", output[-800:] if code else ""):
        got_lags, got = read_ct(work / "case")
        ok = got is not None and len(got) == len(expected) and all(
            abs(g - e) <= 2e-6 for g, e in zip(got, expected))
        checks.check(ok, "C(t) values", f"titan {got[:6] if got else got}..., expected "
                     f"{[round(v, 6) for v in expected[:6]]}... ({len(expected)} lags)")
        ok = got_lags is not None and len(got_lags) == len(lags) and all(
            abs(g - e) <= 1e-6 for g, e in zip(got_lags, lags))
        checks.check(ok, "lag times use the median step (1 ns)",
                     f"titan {got_lags[:3] if got_lags else got_lags}")
        summary = fx.read_key_values(work / "case" / "llps_summary.dat")
        checks.check(summary.get("percolation_time_ns") is not None and
                     abs(float(summary["percolation_time_ns"]) - times[PERCOLATION]) < 1e-3,
                     "percolation at frame 500", str(summary.get("percolation_time_ns")))
        checks.check(summary.get("tau_ns") is not None and abs(float(summary["tau_ns"]) - tau) < 1e-4,
                     "tau", f"titan {summary.get('tau_ns')}, expected {tau:.6f}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
