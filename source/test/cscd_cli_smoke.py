"""End-to-end gate for -cscd (residue centre-of-mass distances).

titan.exe crashed (exit 139) on the first frame of every -cscd run that did
not also ask for -opdbsep; -common_topo runs took the old ligand branch on
top of that (12 Oct 2026, H10 dimer). This runs the same route on a
synthetic two-chain fixture: chain A is molecule 1, chain B molecule 2, five
frames in one multi-model PDB, masses from a -common_topo topology.

1. -cscd: exit 0, one file per chain-A residue, one row per frame at the
   frame's time, and each value equal to the nearest chain-B centre computed
   here from the PDB coordinates and topology masses.
2. -cscd -act: one column per chain-B residue, every pair checked.
3. Two replica trajectories (a list file; chain B 1 A further in the
   second): one row per time, each the mean of the two replicas.
4. -opdbsep with the absolute -f path: the trajectory's own files in ./dimer.
5. -tu=ns: the same rows with times in ns.
6. -b=10: energy rows, -cscd rows and the -nocal H-bond table start at
   10 ps, the first analysed frame, not at the reference frame's 0 ps.

The centres must be side-chain centres of mass: -common_topo marked no side
chain until 12 Oct 2026, so every centre was the whole residue. Chain B uses
CHARMM names (amide HN, C-terminal OT1/OT2) and chain A ends in OXT; all of
them are backbone, and counting any as side chain moves a centre by far more
than the 0.001 A tolerance.
"""
import argparse
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MASS = {"N": 14.007, "C": 12.011, "O": 15.999, "H": 1.008}
BACKBONE = {"N", "CA", "C", "O", "H", "HN", "HA", "H1", "H2", "H3", "HT1", "HT2",
            "HT3", "OXT", "OT1", "OT2", "OC1", "OC2"}
# (residue, [(atom, dx, dy, dz)]) with offsets from the residue origin, A.
RESIDUES = {
    "LYS": [("N", 0.0, 0.0, 0.0), ("H", -0.5, 0.8, 0.0), ("CA", 1.4, 0.2, 0.0),
            ("HA", 1.6, 0.8, 0.9), ("CB", 2.0, 1.0, -1.2), ("CG", 3.4, 1.4, -1.1),
            ("CD", 4.0, 2.2, -2.2), ("CE", 5.4, 2.6, -2.1), ("NZ", 6.0, 3.4, -3.1),
            ("C", 2.0, -1.2, 0.1), ("O", 1.4, -2.2, 0.3)],
    "GLY": [("N", 0.0, 0.0, 0.0), ("H", -0.5, 0.8, 0.0), ("CA", 1.4, 0.2, 0.0),
            ("HA1", 1.6, 0.8, 0.9), ("HA2", 1.7, 0.8, -0.9), ("C", 2.0, -1.2, 0.1),
            ("O", 1.4, -2.2, 0.3)],
    "GLU": [("N", 0.0, 0.0, 0.0), ("H", -0.5, 0.8, 0.0), ("CA", 1.4, 0.2, 0.0),
            ("HA", 1.6, 0.8, 0.9), ("CB", 2.0, 1.0, -1.2), ("CG", 3.4, 1.5, -1.0),
            ("CD", 4.1, 2.3, -2.0), ("OE1", 3.6, 2.6, -3.1), ("OE2", 5.2, 2.7, -1.7),
            ("C", 2.0, -1.2, 0.1), ("O", 1.4, -2.2, 0.3)],
    "ASP": [("N", 0.0, 0.0, 0.0), ("H", -0.5, 0.8, 0.0), ("CA", 1.4, 0.2, 0.0),
            ("HA", 1.6, 0.8, 0.9), ("CB", 2.0, 1.0, -1.2), ("CG", 3.4, 1.5, -1.0),
            ("OD1", 4.0, 2.2, -1.9), ("OD2", 4.0, 1.2, 0.1), ("C", 2.0, -1.2, 0.1),
            ("O", 1.4, -2.2, 0.3)],
    "ALA": [("N", 0.0, 0.0, 0.0), ("H", -0.5, 0.8, 0.0), ("CA", 1.4, 0.2, 0.0),
            ("HA", 1.6, 0.8, 0.9), ("CB", 2.0, 1.0, -1.2), ("C", 2.0, -1.2, 0.1),
            ("O", 1.4, -2.2, 0.3)],
}
CHAINS = [("A", ["LYS", "GLY", "GLU"]), ("B", ["ASP", "ALA", "LYS"])]
TIMES = [0.0, 10.0, 20.0, 30.0, 40.0]


def element(atom):
    return atom[0]


def residue_atoms(chain, r, name):
    """Atoms of residue r (0-based) with this chain's naming: chain B is
    CHARMM-style (HN, C-terminal OT1/OT2), chain A ends in O + OXT."""
    atoms = list(RESIDUES[name])
    last = r == len(dict(CHAINS)[chain]) - 1
    if chain == "B":
        atoms = [("HN",) + a[1:] if a[0] == "H" else a for a in atoms]
        if last:
            atoms = [("OT1",) + a[1:] if a[0] == "O" else a for a in atoms]
            atoms.append(("OT2", 3.2, -1.0, 0.0))
    elif last:
        atoms.append(("OXT", 3.2, -1.0, 0.0))
    return atoms


def frame_atoms(frame, replica=1):
    """Atoms of one frame: chain B slides a little every frame; in replica 2
    it sits 1 A further along x."""
    atoms = []
    for chain, names in CHAINS:
        for r, name in enumerate(names):
            if chain == "A":
                ox, oy, oz = 20.0 + 3.8 * r, 25.0, 25.0
            else:
                ox = 20.0 + 3.8 * (2 - r) + 0.7 * frame + (replica - 1)
                oy, oz = 33.0 - 0.9 * frame, 25.0 + 0.5 * frame
            for atom, dx, dy, dz in residue_atoms(chain, r, name):
                y_sign = 1.0 if chain == "A" else -1.0
                atoms.append((chain, r + 1, name, atom,
                              round(ox + dx, 3), round(oy + y_sign * dy, 3), round(oz + dz, 3)))
    return atoms


def write_inputs(work):
    for replica, name in ((1, "dimer.pdb"), (2, "dimer_rep2.pdb")):
        lines = []
        for f, t in enumerate(TIMES):
            lines.append("CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1")
            lines.append(f"MODEL     {f + 1:4d}")
            lines.append(f"TITLE     t= {t:.5f} step= {f}")
            atoms = frame_atoms(f, replica)
            for serial, (chain, resnr, resname, atom, x, y, z) in enumerate(atoms, 1):
                atom_field = f" {atom:<3s}" if len(atom) < 4 else atom
                lines.append(f"ATOM  {serial:5d} {atom_field} {resname:3s} {chain}{resnr:4d}    "
                             f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element(atom)}")
            lines.append("TER")
            lines.append("ENDMDL")
        (work / name).write_text("\n".join(lines) + "\nEND\n")
    (work / "replicas.txt").write_text(
        (work / "dimer.pdb").as_posix() + "\n" + (work / "dimer_rep2.pdb").as_posix() + "\n")

    top = ["[ defaults ]", "1 2 yes 0.5 0.8333333333", "[ atomtypes ]"]
    for el, mass in MASS.items():
        top.append(f"T{el} {mass} 0.0 A 0.3 0.1")
    serial = 0
    for chain, names in CHAINS:
        top += ["[ moleculetype ]", f"CH{chain} 3", "[ atoms ]"]
        local = 0
        for r, name in enumerate(names):
            for atom, *_ in residue_atoms(chain, r, name):
                local += 1
                el = element(atom)
                top.append(f"{local} T{el} {r + 1} {name} {atom} {local} 0.0 {MASS[el]}")
        serial += local
    top += ["[ system ]", "Two-chain -cscd fixture", "[ molecules ]", "CHA 1", "CHB 1"]
    (work / "dimer.top").write_text("\n".join(top) + "\n")

    n_a = sum(1 for a in frame_atoms(0) if a[0] == "A")
    n = len(frame_atoms(0))
    ndx = ["[ Protein ]", " ".join(str(i) for i in range(1, n + 1)),
           "[ chA ]", " ".join(str(i) for i in range(1, n_a + 1)),
           "[ chB ]", " ".join(str(i) for i in range(n_a + 1, n + 1))]
    (work / "dimer.ndx").write_text("\n".join(ndx) + "\n")
    (work / "lib").mkdir(exist_ok=True)


def centres(frame, side_chain, replica=1):
    out = {}
    atoms = frame_atoms(frame, replica)
    for chain, names in CHAINS:
        for r, name in enumerate(names):
            members = [a for a in atoms if a[0] == chain and a[1] == r + 1]
            if side_chain and name != "GLY":
                side = [a for a in members if a[3] not in BACKBONE]
                members = side or members
            w = [MASS[element(a[3])] for a in members]
            out[(chain, r + 1)] = tuple(
                sum(wi * a[4 + k] for wi, a in zip(w, members)) / sum(w) for k in range(3))
    return out


def run_case(binary, work, name, extra, source=None):
    case = work / name
    case.mkdir(exist_ok=True)
    for old in case.glob("out_cal_dis_com_*.dat"):
        old.unlink()
    lines = (source or ["-f=" + (work / "dimer.pdb").as_posix()]) + [
        "-1pdb=true",
        "-common_topo=" + (work / "dimer.top").as_posix(),
        "-n=" + (work / "dimer.ndx").as_posix(),
        "-ff=AMBER99SB",
        "-g=1-2",
        "-path-a=" + (work / "lib").as_posix() + "/",
        "-path-n=" + (work / "lib").as_posix() + "/",
        "-nt=2",
        "-cscd=true",
    ] + extra
    (case / "arg.in").write_text("\n".join(lines) + "\n")
    result = subprocess.run([str(binary), "-a", "arg.in"], cwd=case, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=600)
    output = result.stdout + result.stderr
    (case / "result.log").write_text(output, encoding="utf-8")
    require(result.returncode == 0, f"{name}: exit {result.returncode}\n{output[-3000:]}")
    require("TITAN ENGINE TERMINATED SUCCESSFULLY" in output, f"{name}: no clean finish")
    return case


def read_table(path):
    header, rows = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            header.append(line)
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    return header, rows


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def require_times(path, times):
    """The first column of every data row (comment lines start with # or @)."""
    require(path.exists(), f"missing {path.name}")
    got = [float(line.split()[0]) for line in path.read_text().splitlines()
           if line.strip() and not line.lstrip().startswith(("#", "@"))]
    require(len(got) == len(times) and all(abs(g - t) < 1e-6 for g, t in zip(got, times)),
            f"{path.parent.name}/{path.name}: times {got}, want {times}")


def check_case(case, all_partners, replicas=(1,), time_scale=1.0):
    for r in range(1, 4):
        path = case / f"out_cal_dis_com_A{r}.dat"
        require(path.exists(), f"missing {path.name}")
        header, rows = read_table(path)
        centre = [line for line in header if line.startswith("# centres: ")]
        require(len(centre) == 1, f"{path.name}: no '# centres:' header line")
        kind = centre[0][len("# centres: "):]
        side_chain = kind.startswith("side chain")
        require(side_chain, f"{path.name}: centres are not side chains: {kind}")
        if len(replicas) > 1:
            require(any(f"# rows: {len(replicas)} trajectories merged by time; "
                        f"{len(TIMES)} of {len(TIMES)} rows average" in line for line in header),
                    f"{path.name}: replicas not merged by time")
        require(len(rows) == len(TIMES), f"{path.name}: {len(rows)} rows, want {len(TIMES)}")
        for f, row in enumerate(rows):
            require(abs(row[0] - TIMES[f] * time_scale) < 1e-6,
                    f"{path.name}: row {f} time {row[0]}")
            # Each replica's distances, then their mean at this time.
            per_replica = []
            for replica in replicas:
                c = centres(f, side_chain, replica)
                d = [math.dist(c[("A", r)], c[("B", j)]) for j in range(1, 4)]
                per_replica.append(d if all_partners else [min(d)])
            want = [sum(v) / len(v) for v in zip(*per_replica)]
            got = row[1:]
            require(len(got) == len(want), f"{path.name}: row {f} has {len(got)} values")
            for g, w in zip(got, want):
                require(abs(g - w) <= 1e-3, f"{path.name}: row {f} {g} vs {w}")
    return kind


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, default=ROOT / "test/tmp/cscd_cli")
    args = parser.parse_args()
    binary = args.binary.resolve()
    work = args.work_root.resolve()
    work.mkdir(parents=True, exist_ok=True)
    write_inputs(work)

    nearest = check_case(run_case(binary, work, "nearest", []), all_partners=False)
    print(f"PASS -cscd nearest ({nearest} centres, 3 residues x {len(TIMES)} frames)")
    every = check_case(run_case(binary, work, "act", ["-act=true"]), all_partners=True)
    print(f"PASS -cscd -act ({every} centres, 3 x 3 pairs x {len(TIMES)} frames)")
    replicas = run_case(binary, work, "replicas", ["-act=true"],
                        source=["-i=" + (work / "replicas.txt").as_posix()])
    check_case(replicas, all_partners=True, replicas=(1, 2))
    print(f"PASS -cscd two replicas averaged at equal times ({len(TIMES)} rows)")

    # -opdbsep with an absolute -f path: the result folder used to be
    # "./" + the absolute path, which cannot be created; now ./dimer.
    separate = run_case(binary, work, "separate", ["-opdbsep=true"])
    for r in range(1, 4):
        own = separate / "dimer" / f"out_cal_dis_com_A{r}.dat"
        require(own.exists(), f"-opdbsep: missing {own.relative_to(work)}")
        require(read_table(own)[1] == read_table(separate / own.name)[1],
                f"-opdbsep: {own.name} differs from the merged file")
    check_case(separate, all_partners=False)
    print("PASS -cscd -opdbsep with an absolute input path (./dimer/)")

    # -tu=ns: frame times written in ns.
    nanoseconds = run_case(binary, work, "nanoseconds", ["-tu=ns"])
    check_case(nanoseconds, all_partners=False, time_scale=1e-3)
    print("PASS -tu=ns writes frame times in ns")

    # -b skips the first model. The first energy row used to carry the time
    # of the reference frame (0 ps) instead of its own, and -nocal H-bond
    # statistics read times past the end of their table.
    late = run_case(binary, work, "late_begin", ["-b=10"])
    require_times(late / "total_energy.dat", TIMES[1:])
    require_times(late / "out_cal_dis_com_A1.dat", TIMES[1:])
    late_nocal = run_case(binary, work, "late_begin_nocal",
                          ["-b=10", "-nocal=true", "-hb=true"])
    require_times(late_nocal / "stat_hydrogenbonds.dat", TIMES[1:])
    print("PASS -b=10: energy, -cscd and -nocal H-bond rows start at 10 ps")
    print("PASS cscd_cli_smoke")
    return 0


if __name__ == "__main__":
    # TITAN's banner is not ASCII; a failure message carries its log tail.
    sys.stdout.reconfigure(errors="backslashreplace")
    try:
        sys.exit(main())
    except RuntimeError as error:
        print("FAIL", error)
        sys.exit(1)
