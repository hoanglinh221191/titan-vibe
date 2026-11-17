"""Oracle for LLPS droplet fusion/fission counting (STR-09).

105 one-bead chains in a 200 A box, DBSCAN (eps 15 A, minpts 3), nine frames:
  frames 0-2  a 100-chain droplet (5 x 5 x 4 grid, 4 A spacing) and a 5-chain
              droplet 50 A above it: two clusters;
  frames 3-5  the small droplet 6 A above the large one: one cluster of 105;
  frames 6-8  the small droplet back at 50 A: two clusters.

A previous cluster belongs to a new cluster when at least half of its
members are in it; two previous clusters in one new cluster are a fusion,
one previous cluster whose members make up at least half of each of two new
clusters is a fission. Expected: one fusion (frame 3), one fission (frame 6),
the 5-chain cluster of frame 6 continues (stable), it is not born again.
The old rule (Jaccard >= 0.30) gave J = 5/105 for the small droplet: no
fusion, no fission, and the small droplet "born" at frame 6.

  python -I test/llps_lineage_oracle.py --binary titan.exe --work-root D:/softs/titan_work/str/llps_lineage_oracle
"""
import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402

BOX = 200.0
LABELS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def frame_atoms(near):
    big = [(48.0 + 4 * i, 48.0 + 4 * j, 54.0 + 4 * k)
           for i in range(5) for j in range(5) for k in range(4)]
    top = 72.0 if near else 116.0
    small = [(56.0, 56.0, top), (52.0, 56.0, top), (60.0, 56.0, top),
             (56.0, 52.0, top), (56.0, 60.0, top)]
    atoms = []
    for n, (x, y, z) in enumerate(big + small):
        # consecutive beads in different chain IDs: one chain per bead
        atoms.append((LABELS[n % len(LABELS)], n + 1, "ALA", "CA", x, y, z))
    # group 2: one far bead (LLPS runs in the energy frame loop, two groups)
    atoms.append(("Z", 999, "LIG", "X1", 190.0, 190.0, 190.0))
    return atoms


def write_topology(path, count):
    """One moleculetype per bead (each its own chain), then the far bead."""
    top = ["[ defaults ]", "1 2 yes 0.5 0.8333333333", "[ atomtypes ]",
           "TC 12.011 0.0 A 0.3 0.1", "TX 72.0 0.0 A 0.3 0.1"]
    for n in range(1, count):
        top += ["[ moleculetype ]", f"B{n} 1", "[ atoms ]", f"1 TC {n} ALA CA 1 0.0 12.011"]
    top += ["[ moleculetype ]", "FAR 1", "[ atoms ]", "1 TX 999 LIG X1 1 0.0 72.0",
            "[ system ]", "lineage oracle", "[ molecules ]"]
    top += [f"B{n} 1" for n in range(1, count)] + ["FAR 1"]
    Path(path).write_text("\n".join(top) + "\n")


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

    near = [False] * 3 + [True] * 3 + [False] * 3
    frames = [(1000.0 * f, (BOX, BOX, BOX), frame_atoms(n)) for f, n in enumerate(near)]
    n_atoms = len(frames[0][2])
    fx.write_trajectory(work / "drops.pdb", frames)
    write_topology(work / "drops.top", n_atoms)
    fx.write_index(work / "drops.ndx", [("System", range(1, n_atoms + 1)),
                                        ("Beads", range(1, n_atoms)), ("Far", [n_atoms])])
    lines = [a for a in fx.base_args(work, work / "drops.pdb", work / "drops.top",
                                     work / "drops.ndx", groups="1-2")
             if not a.startswith("-nocal")]
    lines += ["-mol-1-type=protein", "-mol-2-type=lig", "-llps=true", "-llps-dbscan=true",
              "-llps-dynamics=true", "-llps-dbscan-eps=15", "-llps-dbscan-minpts=3",
              "-llps-density=off"]
    code, output = fx.run_titan(binary, work / "case", lines)
    if not checks.check(code == 0, "titan exit 0", output[-800:] if code else ""):
        return checks.exit_code()
    case = work / "case"
    rows = fx.read_csv(case / "llps_cluster_lineage.csv")
    sizes = {}
    for row in rows:
        sizes.setdefault(round(float(row["Time_ns"])), []).append(int(row["ProteinChains"]))
    expected_sizes = {f: ([105] if n else [100, 5]) for f, n in enumerate(near)}
    checks.check({f: sorted(v, reverse=True) for f, v in sizes.items()} == expected_sizes,
                 "clusters per frame (DBSCAN input as designed)", str(sizes))
    events = {}
    for row in rows:
        events.setdefault(round(float(row["Time_ns"])), []).append(row["Event"])
    counts = fx.read_key_values(case / "llps_fusion_fission.dat")
    checks.check(counts.get("fusion_events") == "1", "one fusion (5-chain droplet into the 100-chain one)",
                 f"fusion_events {counts.get('fusion_events')}")
    checks.check(counts.get("fission_events") == "1", "one fission (the 5-chain droplet leaves)",
                 f"fission_events {counts.get('fission_events')}")
    checks.check(events.get(3) == ["fusion"], "frame 3 event is fusion", str(events.get(3)))
    checks.check(sorted(events.get(6, [])) == ["stable", "stable"],
                 "frame 6: both clusters continue a previous cluster (none born)", str(events.get(6)))
    checks.check(sorted(events.get(0, [])) == ["born", "born"], "frame 0: both born", str(events.get(0)))
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
