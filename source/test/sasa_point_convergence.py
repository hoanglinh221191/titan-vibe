"""Per-residue SASA error of the Shrake-Rupley point count (MDB-15).

Builds test/titan_sasa_point_convergence.cpp against src/titan_surface.cpp,
runs it on the H10 reference structure (test/fixtures/e2e_h10/protein_ref.pdb,
1254 atoms, radii r_min/2 = 2^(1/6) sigma / 2 from the CHARMM36 subset, as
the GBSA path uses) and reports, per point count, the per-residue SASA error
against the largest count and the time of one call (one thread).

  python -I test/sasa_point_convergence.py --work-root D:/softs/titan_work/str/sasa_convergence
"""
import argparse
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "test" / "fixtures" / "e2e_h10"
COUNTS = (96, 192, 384, 768, 1536, 3072, 6144, 12288)
GAMMA = 0.0072  # kcal/mol/A^2, the GBSA default surface tension


def sections(path):
    current = None
    for raw in Path(path).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            current = line.strip("[] ").lower()
            continue
        yield current, line.split()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work-root", type=Path, required=True)
    p.add_argument("--compiler", default="clang++" if os.name == "nt" else "g++")
    a = p.parse_args()
    work = a.work_root.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    sigma = {f[0]: float(f[5]) for sec, f in sections(FIXTURE / "charmm36_subset.ff" / "ffnonbonded.itp")
             if sec == "atomtypes" and len(f) >= 7}
    atoms = []  # (chain, resnr, resname, radius)
    for chain, itp in (("A", "topol_Protein_chain_A.itp"), ("B", "topol_Protein_chain_A2.itp")):
        for sec, f in sections(FIXTURE / itp):
            if sec == "atoms" and len(f) >= 8:
                atoms.append((chain, int(f[2]), f[3], 10.0 * sigma[f[1]] * 2.0 ** (1.0 / 6.0) / 2.0))
    coords = [(float(l[30:38]), float(l[38:46]), float(l[46:54]))
              for l in (FIXTURE / "protein_ref.pdb").read_text().splitlines()
              if l.startswith(("ATOM", "HETATM"))]
    if len(coords) != len(atoms):
        raise SystemExit(f"{len(coords)} coordinates for {len(atoms)} topology atoms")
    (work / "atoms.txt").write_text("".join(f"{x} {y} {z} {r}\n" for (x, y, z), (_, _, _, r)
                                            in zip(coords, atoms)))

    exe = work / ("titan_sasa_point_convergence" + (".exe" if os.name == "nt" else ""))
    subprocess.run([a.compiler, "-std=c++23", "-O2", "-mavx2", "-mfma",
                    "-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1", "-I" + str(ROOT / "src"),
                    str(ROOT / "src" / "titan_surface.cpp"), str(ROOT / "src" / "titan_sasa.cpp"),
                    str(ROOT / "src" / "titan_executor.cpp"),
                    str(ROOT / "test" / "titan_sasa_point_convergence.cpp"), "-o", str(exe)],
                   check=True)
    run = subprocess.run([str(exe), str(work / "atoms.txt"), str(work / "area.csv")]
                         + [str(c) for c in COUNTS], capture_output=True, text=True, check=True)
    timing = {int(r.split(",")[0]): float(r.split(",")[1]) for r in run.stdout.split()[1:]}
    rows = (work / "area.csv").read_text().splitlines()
    per_atom = [[float(v) for v in r.split(",")] for r in rows[1:]]

    residues = {}
    for (chain, resnr, resname, _), values in zip(atoms, per_atom):
        residues.setdefault((chain, resnr, resname), [0.0] * len(COUNTS))
        for k, v in enumerate(values):
            residues[(chain, resnr, resname)][k] += v
    ref = len(COUNTS) - 1
    total_ref = sum(v[ref] for v in residues.values())
    lines = [f"structure: H10 protein_ref.pdb, {len(atoms)} atoms, {len(residues)} residues",
             f"reference: {COUNTS[ref]} points, total SASA {total_ref:.2f} A^2",
             "points,ms_per_call,total_err_A2,res_mean_abs_err_A2,res_p95_abs_err_A2,"
             "res_max_abs_err_A2,res_mean_rel_err_pct(SASA>10),res_max_rel_err_pct(SASA>10),"
             "res_max_gamma_err_kcal"]
    for k, count in enumerate(COUNTS[:-1]):
        errors = sorted(abs(v[k] - v[ref]) for v in residues.values())
        rel = [abs(v[k] - v[ref]) / v[ref] * 100.0 for v in residues.values() if v[ref] > 10.0]
        p95 = errors[min(len(errors) - 1, math.ceil(0.95 * len(errors)) - 1)]
        total = sum(v[k] for v in residues.values())
        lines.append(f"{count},{timing[count]:.2f},{total - total_ref:.3f},"
                     f"{sum(errors) / len(errors):.3f},{p95:.3f},{errors[-1]:.3f},"
                     f"{sum(rel) / len(rel):.2f},{max(rel):.2f},{GAMMA * errors[-1]:.4f}")
    lines.append(f"{COUNTS[ref]},{timing[COUNTS[ref]]:.2f},reference")
    report = "\n".join(lines) + "\n"
    (work / "convergence.csv").write_text(report)
    print(report, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
