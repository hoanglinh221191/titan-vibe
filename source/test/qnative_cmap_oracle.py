"""Oracle for Q-native and the contact map without a periodic box, the soft
Q switch, and the glycine contact point (STR-26: MDB-13, MDB-14, MDB-16).

The ten-residue LYS hairpin of qnative_oracle.py, three frames: the
reference, strand 2 moved 1.6 A away, strand 2 moved 3 A away.

  no_box_*     the frames written without CRYST1, and with a CRYST1 whose
               angles are degenerate (90 90 0): Q(t) and the contact counts
               must be those of plain distances. The old code dropped every
               frame with a <= 0 (empty Q, zero contact probabilities) and
               passed degenerate angles to a NaN cell (Q = 0).
  soft_*       -qnative-soft=true: Q(t) = mean over native pairs of
               1/(1 + exp(beta (r - lambda r0))), beta = 5 1/A, r0 the pair's
               reference distance, lambda 1.5 (ca) or 1.8 (sidechain-heavy).
               The old switch sat at the cutoff for every pair.
  gly_*        residues 3 and 8 made glycine, -qnative-mode/-cmap-mode
               sidechain-com: the glycine point is its CA (the old code took
               the COM of the whole residue). The cutoff is chosen here so
               that the two definitions give different native pairs.

  python -I test/qnative_cmap_oracle.py --binary titan.exe --work-root D:/softs/titan_work/str/qnative_cmap_oracle
"""
import argparse
import itertools
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import md_oracle_fixture as fx  # noqa: E402
import qnative_oracle as qo  # noqa: E402

BOX = qo.BOX
BETA = 5.0
GLY = (3, 8)
BACKBONE = {"N", "CA", "C", "O", "H", "HN", "HA", "H1", "H2", "H3", "OXT"}


def write_trajectory(path, frames, cryst):
    """As md_oracle_fixture.write_trajectory; cryst is a CRYST1 line or None."""
    lines = []
    for index, (time_ps, atoms) in enumerate(frames):
        if cryst:
            lines.append(cryst)
        lines.append(f"MODEL     {index + 1:4d}")
        lines.append(f"TITLE     t= {time_ps:.5f} step= {index}")
        for serial, (chain, resnr, resname, atom, x, y, z) in enumerate(atoms, 1):
            atom_field = f" {atom:<3s}" if len(atom) < 4 else atom
            lines.append(f"ATOM  {serial:5d} {atom_field} {resname:<4s}{chain}{resnr:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {fx.element(atom)}")
        lines.append("TER")
        lines.append("ENDMDL")
    Path(path).write_text("\n".join(lines) + "\nEND\n")


def hairpin_frames(gly=False):
    frames = []
    for time_ps, shift in ((0.0, 0.0), (10.0, 1.6), (20.0, 3.0)):
        atoms = qo.hairpin(shift)
        if gly:
            atoms = glycine_variant(atoms)
        frames.append((time_ps, atoms))
    return frames


def glycine_variant(atoms):
    out = []
    for chain, resnr, resname, name, x, y, z in atoms:
        if resnr not in GLY:
            out.append((chain, resnr, resname, name, x, y, z))
            continue
        if name not in BACKBONE:
            continue
        if name == "HA":
            out.append((chain, resnr, "GLY", "HA2", x, y, z))
            out.append((chain, resnr, "GLY", "HA3", x - 0.6, y + 0.4, z + 0.8))
        else:
            out.append((chain, resnr, "GLY", name, x, y, z))
    return out


def points(atoms, mode, gly_rule="ca"):
    """One or more contact points per residue, as titan's contact geometry."""
    residues = {}
    for chain, resnr, resname, name, x, y, z in atoms:
        residues.setdefault(resnr, []).append((resname, name, (x, y, z)))
    out = {}
    for resnr, members in residues.items():
        gly = members[0][0] == "GLY"
        if mode == "ca":
            out[resnr] = [p for _, n, p in members if n == "CA"]
        elif mode == "sidechain-heavy":
            out[resnr] = [p for _, n, p in members
                          if (n == "CA" if gly else (not n.startswith("H") and n not in BACKBONE))]
        else:  # sidechain-com
            if gly and gly_rule == "ca":
                out[resnr] = [p for _, n, p in members if n == "CA"]
                continue
            chosen = [(n, p) for _, n, p in members if gly or n not in BACKBONE]
            mass = [fx.MASS[fx.element(n)] for n, _ in chosen]
            total = sum(mass)
            out[resnr] = [tuple(sum(m * p[k] for m, (_, p) in zip(mass, chosen)) / total
                                for k in range(3))]
    return out


def residue_distance(pts, i, j):
    return min(math.dist(p, q) for p in pts[i] for q in pts[j])


def natives_of(ref_pts, cutoff):
    return {(i, j): residue_distance(ref_pts, i, j)
            for i, j in itertools.combinations(sorted(ref_pts), 2)
            if j - i >= 3 and residue_distance(ref_pts, i, j) <= cutoff}


def q_series(frames, mode, cutoff, soft=None, gly_rule="ca"):
    natives = natives_of(points(frames[0][1], mode, gly_rule), cutoff)
    series = []
    for _, atoms in frames:
        pts = points(atoms, mode, gly_rule)
        total = 0.0
        for (i, j), r0 in natives.items():
            r = residue_distance(pts, i, j)
            if soft is None:
                total += 1.0 if r <= cutoff else 0.0
            else:
                total += 1.0 / (1.0 + math.exp(BETA * (r - soft * r0)))
        series.append(total / len(natives))
    return sorted(natives), series


def contact_counts(frames, mode, cutoff, gly_rule="ca"):
    counts = []
    for _, atoms in frames:
        pts = points(atoms, mode, gly_rule)
        counts.append(sum(1 for i, j in itertools.combinations(sorted(pts), 2)
                          if j - i >= 3 and residue_distance(pts, i, j) <= cutoff))
    return counts


def gly_cutoff(frames):
    """A sidechain-com cutoff at which CA and whole-residue points for the
    glycines give different native sets, every reference distance at least
    0.05 A from it."""
    new = points(frames[0][1], "sidechain-com", "ca")
    old = points(frames[0][1], "sidechain-com", "com")
    distances = []
    for i, j in itertools.combinations(sorted(new), 2):
        if j - i >= 3:
            distances += [residue_distance(new, i, j), residue_distance(old, i, j)]
    best = None
    for step in range(400, 900):
        cut = step / 100.0
        margin = min(abs(d - cut) for d in distances)
        if margin < 0.05 or natives_of(new, cut).keys() == natives_of(old, cut).keys():
            continue
        if best is None or margin > best[1]:
            best = (cut, margin)
    return best[0]


def read_q(case):
    files = [f for f in sorted(case.glob("qnative_timeseries*.csv")) if "average" not in f.name]
    return [float(r["Q"]) for r in fx.read_csv(files[0])] if files else None


def read_natives(case):
    files = [f for f in sorted(case.glob("qnative_native_contacts*.csv")) if "average" not in f.name]
    if not files:
        return None
    return sorted(tuple(sorted((int(r["residue_i"].split(":")[1]), int(r["residue_j"].split(":")[1]))))
                  for r in fx.read_csv(files[0]))


def read_counts(case):
    files = [f for f in sorted(case.glob("contact_count_timeseries*.csv")) if "average" not in f.name]
    return [int(float(r["contact_pairs"])) for r in fx.read_csv(files[0])] if files else None


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

    plain = hairpin_frames()
    gly = hairpin_frames(gly=True)
    cubic = f"CRYST1{BOX:9.3f}{BOX:9.3f}{BOX:9.3f}  90.00  90.00  90.00 P 1           1"
    degenerate = f"CRYST1{BOX:9.3f}{BOX:9.3f}{BOX:9.3f}  90.00  90.00   0.00 P 1           1"
    for name, frames, cryst in (("plain_nobox", plain, None), ("plain_degenerate", plain, degenerate),
                                ("plain_box", plain, cubic), ("gly_box", gly, cubic)):
        write_trajectory(work / f"{name}.pdb", frames, cryst)
    fx.write_topology(work / "plain.top", plain[0][1])
    fx.write_index(work / "plain.ndx", [("System", range(1, len(plain[0][1]) + 1))])
    fx.write_topology(work / "gly.top", gly[0][1])
    fx.write_index(work / "gly.ndx", [("System", range(1, len(gly[0][1]) + 1))])
    cut_gly = gly_cutoff(gly)
    print(f"gly cutoff {cut_gly:.2f} A")

    cases = []
    for traj in ("plain_nobox", "plain_degenerate"):
        for mode, cut in (("ca", 6.5), ("sidechain-heavy", 6.5)):
            cases.append((f"{traj}_{mode}", traj, "plain", mode, cut, None))
    cases.append(("soft_ca", "plain_box", "plain", "ca", 6.5, 1.5))
    cases.append(("soft_sidechain_heavy", "plain_box", "plain", "sidechain-heavy", 6.5, 1.8))
    cases.append(("gly_sidechain_com", "gly_box", "gly", "sidechain-com", cut_gly, None))

    for name, traj, system, mode, cut, soft in cases:
        frames = gly if system == "gly" else plain
        natives, q_expected = q_series(frames, mode, cut, soft)
        counts_expected = contact_counts(frames, mode, cut)
        if soft is not None:
            _, q_old = q_series(frames, mode, cut, None)
            old = [sum(1.0 / (1.0 + math.exp(BETA * (residue_distance(points(atoms, mode), i, j) - cut)))
                       for i, j in natives) / len(natives) for _, atoms in frames]
            checks.check(max(abs(x - y) for x, y in zip(old, q_expected)) > 0.05,
                         f"{name}: switch at lambda*r0 differs from the switch at the cutoff",
                         f"{[round(v, 4) for v in q_expected]} vs {[round(v, 4) for v in old]}")
        if system == "gly":
            old_natives, _ = q_series(frames, mode, cut, None, gly_rule="com")
            checks.check(old_natives != natives, f"{name}: CA and whole-residue COM give other natives",
                         f"{natives} vs {old_natives}")
        lines = fx.base_args(work, work / f"{traj}.pdb", work / f"{system}.top",
                             work / f"{system}.ndx")
        lines += ["-qnative=true", f"-qnative-mode={mode}", f"-qnative-cut={cut}",
                  "-qnative-min-seqsep=3", "-qnative-domain=both",
                  "-cmap=true", f"-cmap-mode={mode}", f"-cmap-cut={cut}", "-cmap-min-seqsep=3"]
        if soft is not None:
            lines.append("-qnative-soft=true")
        code, output = fx.run_titan(binary, work / name, lines)
        if not checks.check(code == 0, f"{name}: titan exit 0", output[-600:] if code else ""):
            continue
        got_natives = read_natives(work / name)
        checks.check(got_natives == natives, f"{name}: native pairs", f"titan {got_natives}, expected {natives}")
        q = read_q(work / name)
        ok = q is not None and len(q) == len(q_expected) and all(
            abs(g - e) <= 2e-6 for g, e in zip(q, q_expected))
        checks.check(ok, f"{name}: Q(t)", f"titan {q}, expected {[round(v, 6) for v in q_expected]}")
        counts = read_counts(work / name)
        checks.check(counts == counts_expected, f"{name}: contact-map counts per frame",
                     f"titan {counts}, expected {counts_expected}")
    return checks.exit_code()


if __name__ == "__main__":
    sys.exit(main())
