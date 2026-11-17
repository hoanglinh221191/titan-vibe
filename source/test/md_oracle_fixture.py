"""Shared writer for the small synthetic systems of the MD-analysis oracles
(test/*_oracle.py for STR-20...STR-25): a multi-model PDB trajectory with a
CRYST1 box and a "t=" title per frame, a GROMACS topology that carries the
masses (-common_topo), an index file, and one titan run per case directory.

Atoms are tuples (chain, resnr, resname, atom, x, y, z) in Angstrom. The
element is the first letter of the atom name; masses from MASS.
"""
import math
import subprocess
from pathlib import Path

MASS = {"N": 14.007, "C": 12.011, "O": 15.999, "H": 1.008, "S": 32.06, "P": 30.974,
        "X": 72.0}


def element(atom):
    """First letter of the name; coarse-grained beads (GL1, ROH, ...) whose
    letter is not an element here are X, mass 72."""
    name = atom.lstrip("0123456789")
    return name[0] if name[0] in MASS else "X"


def wrap(value, length):
    """Coordinate folded into [0, length)."""
    return value - math.floor(value / length) * length


def write_trajectory(path, frames):
    """frames: list of (time_ps, (a, b, c), atoms)."""
    lines = []
    for index, (time_ps, box, atoms) in enumerate(frames):
        a, b, c = box
        lines.append(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}  90.00  90.00  90.00 P 1           1")
        lines.append(f"MODEL     {index + 1:4d}")
        lines.append(f"TITLE     t= {time_ps:.5f} step= {index}")
        for serial, (chain, resnr, resname, atom, x, y, z) in enumerate(atoms, 1):
            atom_field = f" {atom:<3s}" if len(atom) < 4 else atom
            lines.append(f"ATOM  {serial:5d} {atom_field} {resname:<4s}{chain}{resnr:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element(atom)}")
        lines.append("TER")
        lines.append("ENDMDL")
    Path(path).write_text("\n".join(lines) + "\nEND\n")


def write_topology(path, atoms, title="MD oracle fixture"):
    """One moleculetype per chain, in the order of the atoms; zero charges."""
    top = ["[ defaults ]", "1 2 yes 0.5 0.8333333333", "[ atomtypes ]"]
    for el, mass in MASS.items():
        top.append(f"T{el} {mass} 0.0 A 0.3 0.1")
    chains = []
    for atom in atoms:
        if not chains or chains[-1][0] != atom[0]:
            chains.append((atom[0], []))
        chains[-1][1].append(atom)
    for chain, members in chains:
        top += ["[ moleculetype ]", f"CH{chain} 3", "[ atoms ]"]
        residues = []
        for local, (_, resnr, resname, name, *_) in enumerate(members, 1):
            if not residues or residues[-1] != resnr:
                residues.append(resnr)
            el = element(name)
            top.append(f"{local} T{el} {resnr} {resname} {name} {local} 0.0 {MASS[el]}")
    top += ["[ system ]", title, "[ molecules ]"]
    top += [f"CH{chain} 1" for chain, _ in chains]
    Path(path).write_text("\n".join(top) + "\n")


def write_index(path, groups):
    """groups: list of (name, [1-based atom numbers])."""
    lines = []
    for name, members in groups:
        lines.append(f"[ {name} ]")
        lines.append(" ".join(str(i) for i in members))
    Path(path).write_text("\n".join(lines) + "\n")


def run_titan(binary, case, lines, timeout=600):
    """Writes case/arg.in and runs `binary -a arg.in` in case. Returns
    (exit code, stdout + stderr)."""
    case = Path(case)
    case.mkdir(parents=True, exist_ok=True)
    (case / "arg.in").write_text("\n".join(lines) + "\n")
    result = subprocess.run([str(binary), "-a", "arg.in"], cwd=case, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=timeout)
    output = result.stdout + result.stderr
    (case / "result.log").write_text(output, encoding="utf-8")
    return result.returncode, output


def base_args(work, traj, top, ndx, groups="0-0"):
    work = Path(work)
    (work / "lib").mkdir(exist_ok=True)
    return ["-f=" + Path(traj).as_posix(), "-1pdb=true",
            "-common_topo=" + Path(top).as_posix(), "-n=" + Path(ndx).as_posix(),
            "-ff=AMBER99SB", "-g=" + groups,
            "-path-a=" + (work / "lib").as_posix() + "/",
            "-path-n=" + (work / "lib").as_posix() + "/", "-nt=1", "-nocal=true"]


def read_csv(path):
    rows = Path(path).read_text().splitlines()
    header = rows[0].split(",")
    return [dict(zip(header, row.split(","))) for row in rows[1:] if row.strip()]


def read_key_values(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            out[key.strip()] = value.strip()
    return out


class Checks:
    """Collects PASS/FAIL lines; exit code 1 if anything failed."""

    def __init__(self):
        self.failed = 0

    def check(self, condition, label, detail=""):
        status = "PASS" if condition else "FAIL"
        if not condition:
            self.failed += 1
        print(f"[{status}] {label}" + (f": {detail}" if detail else ""), flush=True)
        return condition

    def exit_code(self):
        print("RESULT:", "FAIL" if self.failed else "PASS", f"({self.failed} failed)")
        return 1 if self.failed else 0
