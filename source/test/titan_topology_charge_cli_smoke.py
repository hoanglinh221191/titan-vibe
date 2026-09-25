#!/usr/bin/env python3
"""Real TITAN CLI regression; synthetic, deliberately distinct topology/lib charges.

No MD. Every run uses a new directory. The tested binary, commands, raw audit
tables and logs are retained. This is a software test, not B1 replay evidence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(root: Path) -> dict[str, Path]:
    root.mkdir(parents=True)
    # (atom, residue, lib type, topo type, lib charge, topo charge, topo mass)
    atoms = [
        ("R1", "RCP", "LR1", "TR1", 0.11, 1.25, 21.25),
        ("R2", "RCP", "LR2", "TR2", -0.22, -0.15, 22.50),
        ("L1", "LIG", "LL1", "TL1", 0.33, -0.40, 31.20),
        ("L2", "LIG", "LL2", "TL2", -0.44, 0.30, 32.30),
    ]
    table = "[ atomtypes ]\n; name atnum mass charge ptype sigma epsilon\n"
    for i, atom in enumerate(atoms):
        table += f"{atom[2]} 6 {12.0+i} 0.0 A {0.30+i*0.01:.3f} {0.20+i*0.02:.3f}\n"
        table += f"{atom[3]} 6 {40.0+i} 0.0 A {0.40+i*0.01:.3f} {0.60+i*0.02:.3f}\n"
    ff = root / "forcefield" / "cli-test.ff" / "forcefield.itp"
    write(ff, "[ defaults ]\n1 2 yes 1.0 1.0\n" + table)
    lib = root / "lib_ff"
    write(lib / "ffCHARMM36nonbonded.dat", table)
    rtp = "[ bondedtypes ]\n1 1 1 2 1 3 1 0\n"
    for residue in ("RCP", "LIG"):
        rtp += f"[ {residue} ]\n[ atoms ]\n"
        for i, atom in enumerate(atoms):
            if atom[1] == residue:
                rtp += f"{atom[0]} {atom[2]} {atom[4]} {i}\n"
    write(lib / "CHARMM36_aminoacids.rtp", rtp)
    pdb = "TITLE     Synthetic topology mapping regression t=0\n"
    pdb += "CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1\n"
    for i, atom in enumerate(atoms):
        res = 1 if i < 2 else 2
        x = 5.0 + 1.2 * (i % 2) + (10.0 if i >= 2 else 0.0)
        pdb += (f"ATOM  {i+1:5d} {atom[0]:<4s} {atom[1]:3s} A{res:4d}    "
                f"{x:8.3f}{5.0:8.3f}{5.0:8.3f}{1.0:6.2f}{0.0:6.2f}           C\n")
    write(root / "fixture.pdb", pdb + "END\n")
    write(root / "index.ndx", "[ RCP ]\n1 2\n[ LIG ]\n3 4\n")
    blocks = {}
    for residue in ("RCP", "LIG"):
        block = f"[ moleculetype ]\n{residue} 3\n[ atoms ]\n"
        for k, atom in enumerate(a for a in atoms if a[1] == residue):
            block += f"{k+1} {atom[3]} 1 {residue} {atom[0]} {k+1} {atom[5]} {atom[6]}\n"
        blocks[residue] = block + "[ bonds ]\n1 2 1\n"
    include = '#include "cli-test.ff/forcefield.itp"\n'
    for residue in blocks:
        write(root / f"{residue}.top", include + blocks[residue] +
              f"[ system ]\nSynthetic {residue}\n[ molecules ]\n{residue} 1\n")
    write(root / "common.top", include + blocks["RCP"] + blocks["LIG"] +
          "[ system ]\nSynthetic pair\n[ molecules ]\nRCP 1\nLIG 1\n")
    return {"pdb": root / "fixture.pdb", "ndx": root / "index.ndx",
            "lib": lib, "gmxlib": root / "forcefield", "common": root / "common.top",
            "rec": root / "RCP.top", "lig": root / "LIG.top"}


def check_rows(path: Path, topology: tuple[bool, bool], single: bool = False) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != (2 if single else 4):
        raise AssertionError(f"Expected complete selected parameters, got {len(rows)}")
    for i, row in enumerate(rows):
        source = topology[i // 2]
        expected_type = (("TR1", "TR2", "TL1", "TL2") if source else
                         ("LR1", "LR2", "LL1", "LL2"))[i]
        exact = {
            "molecule": "mol1" if i < 2 else "mol2",
            "selection_index_0based": str(i % 2), "source_atom_id_0based": str(i),
            "has_source_id": "1", "atom": ("R1", "R2", "L1", "L2")[i],
            "atom_type": expected_type, "assignment_policy": "topology" if source else "lib_ff",
        }
        for key, value in exact.items():
            if row[key] != value:
                raise AssertionError(f"atom {i} {key}: {row[key]!r} != {value!r}")
        sigma = (0.40 if source else 0.30) + 0.01 * i
        numeric = {
            "charge_e": ((1.25, -0.15, -0.4, 0.3) if source else (0.11, -0.22, 0.33, -0.44))[i],
            "mass_u": (21.25, 22.50, 31.20, 32.30)[i] if source else 12.0 + i,
            "sigma_nm": sigma,
            "epsilon_kjmol": (0.60 if source else 0.20) + 0.02 * i,
            "radius_A": 5 * 2**(1/6) * sigma,
        }
        for key, value in numeric.items():
            actual = float(row[key])
            if not math.isfinite(actual) or not math.isclose(actual, value, rel_tol=1e-10, abs_tol=1e-10):
                raise AssertionError(f"atom {i} {key}: {actual} != {value}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True, type=Path)
    ap.add_argument("--work-root", required=True, type=Path)
    args = ap.parse_args()
    binary, root = args.binary.resolve(), args.work_root.resolve()
    if root.exists():
        raise RuntimeError(f"Refusing an existing regression directory: {root}")
    f = fixture(root / "fixture")
    cases = [
        ("library_only", [], (False, False), False),
        ("topo_ligand", ["-topo", str(f["lig"])], (False, True), False),
        ("topo2_receptor", ["-topo2", str(f["rec"])], (True, False), False),
        ("separate_both", ["-topo", str(f["lig"]), "-topo2", str(f["rec"])], (True, True), False),
        ("common_both", ["-common_topo", str(f["common"])], (True, True), False),
        ("common_alias", ["-common-topo", str(f["common"])], (True, True), False),
        ("implicit_common", ["-topo", str(f["common"]), "-topo2", str(f["common"])], (True, True), False),
        ("common_single", ["-common_topo", str(f["common"]), "-nocal"], (True, False), True),
    ]
    results = []
    for name, extra, topology, single in cases:
        work = root / name
        work.mkdir()
        audit = work / "parameters.tsv"
        command = [str(binary), "-f", str(f["pdb"]), "-ff", "CHARMM36", "-nt", "1",
                   "-g", "0-0" if single else "0-1", "-n", str(f["ndx"]),
                   "-path-a", str(f["lib"]), "-path-n", str(f["lib"]),
                   "-fflib-hydrogens", "none", "-fflib-missing", "strict", "-no-dssp", *extra]
        env = os.environ.copy()
        env.update(TITAN_PARAMETER_AUDIT=str(audit), GMXLIB=str(f["gmxlib"]),
                   TITAN_NO_PAUSE="1", OMP_NUM_THREADS="1", TITAN_THREADS="1")
        run = subprocess.run(command, cwd=work, env=env, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                             timeout=60, check=False)
        write(work / "command.json", json.dumps(command, indent=2) + "\n")
        write(work / "stdout_stderr.log", run.stdout)
        if run.returncode != 0:
            raise RuntimeError(f"{name}: CLI exit={run.returncode}; see {work}")
        check_rows(audit, topology, single)
        results.append({"case": name, "status": "PASS", "audit_sha256": digest(audit)})
        print(f"PASS topology_charge_cli {name}", flush=True)
    # The diagnostic must never overwrite an existing evidence file.
    before = digest(audit)
    run = subprocess.run(command, cwd=work, env=env, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                         timeout=60, check=False)
    write(work / "refuse_existing_audit.log", run.stdout)
    if run.returncode == 0 or digest(audit) != before:
        raise AssertionError("Parameter audit overwrote existing evidence or did not fail closed")
    results.append({"case": "refuse_existing_audit", "status": "PASS"})
    report = {"status": "PASS", "kind": "SYNTHETIC_SOFTWARE_REGRESSION_NOT_B1_REPLAY",
              "binary_sha256": digest(binary), "cases": results, "md_executed": False}
    write(root / "report.json", json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"TITAN_TOPOLOGY_CHARGE_CLI_SMOKE_PASS cases={len(results)} md_executed=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
