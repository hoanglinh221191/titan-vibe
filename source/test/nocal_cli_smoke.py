"""End-to-end gate for the -nocal analyses of src/nocal/.

The two-chain fixture of cscd_cli_smoke.py (chain A = molecule 1, chain B =
molecule 2, models at 0-40 ps), with charges on the LYS, GLU and ASP side
chains, run with -nocal and -b=10 so the first model is skipped:

- -com-dis: one row per frame at the frame's own time (the old code wrote
  the reference frame's time on every row), d_COM and the per-axis minima
  against centres of mass computed here;
- -dipole-moment: one row per frame, sum q (r - geometric centre);
- -charged_atom_dis: three distributions written once per input, each bin's
  share equal to the histogram computed here;
- -resdis2com: mean and SD per residue over the analysed frames;
- -cal-volume: one row per frame, 216000 A^3, and the per-folder statistics;
- -angle_dipole_efield: no water in the fixture, so a header and no rows;
- two inputs through a list file: each input gets its own dipole file and
  the COM-distance file holds both, in input order.
"""
import argparse
import importlib.util
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cscd_cli_smoke",
                                              Path(__file__).with_name("cscd_cli_smoke.py"))
fixture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixture)

CHARGES = {("LYS", "NZ"): 1.0, ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
           ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5}
ANALYSED = fixture.TIMES[1:]  # -b=10


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def write_charged_topology(work):
    top = ["[ defaults ]", "1 2 yes 0.5 0.8333333333", "[ atomtypes ]"]
    for el, mass in fixture.MASS.items():
        top.append(f"T{el} {mass} 0.0 A 0.3 0.1")
    for chain, names in fixture.CHAINS:
        top += ["[ moleculetype ]", f"CH{chain} 3", "[ atoms ]"]
        local = 0
        for r, name in enumerate(names):
            for atom, *_ in fixture.residue_atoms(chain, r, name):
                local += 1
                el = fixture.element(atom)
                q = CHARGES.get((name, atom), 0.0)
                top.append(f"{local} T{el} {r + 1} {name} {atom} {local} {q} {fixture.MASS[el]}")
    top += ["[ system ]", "Two-chain -nocal fixture", "[ molecules ]", "CHA 1", "CHB 1"]
    (work / "dimer.top").write_text("\n".join(top) + "\n")


def run(binary, work, name, source, nocal=True):
    case = work / name
    case.mkdir(exist_ok=True)
    lines = source + [
        "-1pdb=true",
        "-common_topo=" + (work / "dimer.top").as_posix(),
        "-n=" + (work / "dimer.ndx").as_posix(),
        "-ff=AMBER99SB",
        "-g=1-2",
        "-path-a=" + (work / "lib").as_posix() + "/",
        "-path-n=" + (work / "lib").as_posix() + "/",
        "-nt=2",
        "-b=10",
        "-com-dis=xyz",
    ]
    if nocal:
        lines += [
            "-nocal=true",
            "-dipole-moment=true",
            "-charged_atom_dis=true",
            "-resdis2com=true",
            "-cal-volume=true",
            "-angle_dipole_efield=true",
        ]
    (case / "arg.in").write_text("\n".join(lines) + "\n")
    result = subprocess.run([str(binary), "-a", "arg.in"], cwd=case, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=600)
    output = result.stdout + result.stderr
    (case / "result.log").write_text(output, encoding="utf-8")
    require(result.returncode == 0, f"{name}: exit {result.returncode}\n{output[-3000:]}")
    require("TITAN ENGINE TERMINATED SUCCESSFULLY" in output, f"{name}: no clean finish")
    return case


def rows(path):
    require(path.exists(), f"missing {path.name}")
    return [[float(v) for v in line.split()] for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]


def atoms_of(frame, chain, replica=1):
    return [a for a in fixture.frame_atoms(frame, replica) if a[0] == chain]


def mass_centre(atoms):
    w = [fixture.MASS[fixture.element(a[3])] for a in atoms]
    return [sum(wi * a[4 + k] for wi, a in zip(w, atoms)) / sum(w) for k in range(3)]


def close(a, b, tol, what):
    require(abs(a - b) <= tol, f"{what}: {a} vs {b}")


def check_input(case, stem, replica=1):
    frames = [f for f, t in enumerate(fixture.TIMES) if t in ANALYSED]
    # Dipole: sum q (r - geometric centre) over chain A.
    raw = rows(case / f"{stem}_notnormalized__dipole_moment.dat")
    unit = rows(case / f"{stem}_dipole_moment.dat")
    require(len(raw) == len(frames) and len(unit) == len(frames), "dipole: one row per frame")
    for row, f in zip(raw, frames):
        a = atoms_of(f, "A", replica)
        centre = [sum(x[4 + k] for x in a) / len(a) for k in range(3)]
        m = [sum(CHARGES.get((x[2], x[3]), 0.0) * (x[4 + k] - centre[k]) for x in a)
             for k in range(3)]
        close(row[0], fixture.TIMES[f], 1e-6, "dipole time")
        for k in range(3):
            close(row[1 + k], m[k], 1e-3, f"dipole M{k}")
        close(row[4], math.sqrt(sum(c * c for c in m)), 1e-3, "dipole length")

    # Charged-atom distances: 400 bins over [0, 50) A, percent per bin.
    for sign, keep in (("", lambda q: q != 0), ("_positive", lambda q: q > 0),
                       ("_negative", lambda q: q < 0)):
        hist = [0] * 400
        for f in frames:
            a = atoms_of(f, "A", replica)
            com = mass_centre(a)
            for x in a:
                q = CHARGES.get((x[2], x[3]), 0.0)
                if q != 0 and keep(q):
                    d = math.dist(com, x[4:7])
                    hist[min(int(d / 0.125), 399)] += 1
        total = sum(hist)
        got = rows(case / f"{stem}{sign}_distributionChargedAtoms.dat")
        require(len(got) == 400, f"charged{sign}: {len(got)} bins")
        for i, (row, count) in enumerate(zip(got, hist)):
            close(row[1], 100.0 * count / total, 2e-5, f"charged{sign} bin {i}")

    # Residue-to-COM distances of chain A, mean and SD over the frames.
    got = {line.split()[0]: [float(v) for v in line.split()[1:]]
           for line in (case / f"{stem}_distances_res_2COM.dat").read_text().splitlines()
           if line.strip() and not line.startswith("#")}
    names = fixture.CHAINS[0][1]
    for r in range(len(names)):
        d = []
        for f in frames:
            a = atoms_of(f, "A", replica)
            d.append(math.dist(mass_centre([x for x in a if x[1] == r + 1]), mass_centre(a)))
        mean = sum(d) / len(d)
        sd = math.sqrt(sum((v - mean) ** 2 for v in d) / len(d))
        label = f"A{r + 1}"
        require(label in got, f"resdis2com: no row {label}")
        close(got[label][0], mean, 1e-3, f"resdis2com {label} mean")
        close(got[label][1], sd, 1e-3, f"resdis2com {label} SD")

    efield = case / f"{stem}_distributionDipoleEfield.dat"
    require(efield.exists() and not rows(efield), "angle_dipole_efield: header only")


def check_com(case, replicas):
    got = rows(case / "out_dis_com.dat")
    frames = [f for f, t in enumerate(fixture.TIMES) if t in ANALYSED]
    want = [(r, f) for r in replicas for f in frames]
    require(len(got) == len(want), f"com-dis: {len(got)} rows, want {len(want)}")
    for row, (r, f) in zip(got, want):
        a, b = atoms_of(f, "A", r), atoms_of(f, "B", r)
        ca, cb = mass_centre(a), mass_centre(b)
        close(row[0], fixture.TIMES[f], 1e-6, "com-dis time")
        close(row[1], math.dist(ca, cb), 1e-3, "com-dis d_COM")
        for k in range(3):
            close(row[2 + k], min(abs(x[4 + k] - ca[k]) for x in b), 1e-3, f"com-dis min{k}")


def check_volume(case, n_rows):
    got = rows(case / "volume.dat")
    require(len(got) == n_rows and all(abs(v[1] - 216000.0) < 1e-6 for v in got),
            f"volume.dat: {got}")
    stats = (case / "stat_volumes.dat").read_text().splitlines()
    require(len(stats) == 2 and abs(float(stats[1].split()[1]) - 216000.0) < 1e-6
            and float(stats[1].split()[2]) == 0.0, f"stat_volumes.dat: {stats}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, default=ROOT / "test/tmp/nocal_cli")
    args = parser.parse_args()
    binary = args.binary.resolve()
    work = args.work_root.resolve()
    work.mkdir(parents=True, exist_ok=True)
    fixture.write_inputs(work)
    write_charged_topology(work)

    one = run(binary, work, "single", ["-f=" + (work / "dimer.pdb").as_posix()])
    check_input(one, "dimer")
    check_com(one, (1,))
    check_volume(one, len(ANALYSED))
    print(f"PASS -nocal analyses on one input ({len(ANALYSED)} frames from 10 ps)")

    two = run(binary, work, "replicas", ["-i=" + (work / "replicas.txt").as_posix()])
    check_input(two, "dimer", 1)
    check_input(two, "dimer_rep2", 2)
    check_com(two, (1, 2))
    check_volume(two, 2 * len(ANALYSED))
    print("PASS -nocal analyses on two inputs (own files per input, shared COM file)")

    # The energy mode runs -com-dis through the same module.
    energy = run(binary, work, "energy", ["-f=" + (work / "dimer.pdb").as_posix()],
                 nocal=False)
    check_com(energy, (1,))
    print("PASS -com-dis in the energy mode")
    print("PASS nocal_cli_smoke")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(errors="backslashreplace")
    try:
        sys.exit(main())
    except RuntimeError as error:
        print("FAIL", error)
        sys.exit(1)
