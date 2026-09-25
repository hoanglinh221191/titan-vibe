#!/usr/bin/env python3
"""Verify the literature-seeded 4W52 H/J route from pinned coordinates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


INNER_RESIDUES = (102, 111)
OUTER_RESIDUES = (126, 153, 154)
EXTENSIONS_A = (8.0, 18.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return digest


def mean(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    return tuple(sum(point[i] for point in points) / len(points) for i in range(3))  # type: ignore[return-value]


def add(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, float, float]:
    return tuple(a[i] + b[i] for i in range(3))  # type: ignore[return-value]


def scale(a: tuple[float, ...], factor: float) -> tuple[float, float, float]:
    return tuple(value * factor for value in a)  # type: ignore[return-value]


def subtract(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, float, float]:
    return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]


def normalize(a: tuple[float, ...]) -> tuple[float, float, float]:
    norm = math.sqrt(sum(value * value for value in a))
    if norm <= 0.0:
        raise RuntimeError("cannot normalize a zero route direction")
    return scale(a, 1.0 / norm)


def parse_pdb(path: Path) -> tuple[dict[int, tuple[float, float, float]], list[tuple[str, tuple[float, float, float]]]]:
    ca: dict[int, tuple[float, float, float]] = {}
    ligand: list[tuple[str, tuple[float, float, float]]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        residue = line[17:21].strip()
        residue_number = int(line[22:26])
        position = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if atom_name == "CA":
            ca[residue_number] = position
        if residue == "BNZ":
            ligand.append((atom_name, position))
    if len(ligand) != 12:
        raise RuntimeError(f"expected 12 BNZ atoms in {path}, found {len(ligand)}")
    return ca, ligand


def parse_ligand_masses(path: Path) -> dict[str, float]:
    masses: dict[str, float] = {}
    in_atoms = False
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("["):
            in_atoms = line.strip("[] ").lower() == "atoms"
            continue
        if not in_atoms or not line:
            continue
        fields = line.split()
        if len(fields) >= 8 and fields[0].isdigit():
            masses[fields[4]] = float(fields[7])
    if len(masses) != 12:
        raise RuntimeError(f"expected 12 ligand masses in {path}, found {len(masses)}")
    return masses


def route_points(path: Path) -> list[tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        fields = raw.strip().split()
        if fields and fields[0] == "waypoint" and len(fields) == 7:
            points.append(tuple(float(value) for value in fields[3:6]))  # type: ignore[arg-type]
    if len(points) != 5:
        raise RuntimeError(f"expected five route waypoints in {path}")
    return points


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    fixture = Path(__file__).resolve().parent
    template = (
        root
        / ".codex_tmp"
        / "DTW_Path_Classification"
        / "T4L_benzene"
        / "opes_explore"
        / "template"
    )
    parser.add_argument("--pdb", type=Path, default=template / "prot_lig.pdb")
    parser.add_argument("--ligand-itp", type=Path, default=template / "BNZ_GMX.itp")
    parser.add_argument("--route", type=Path, default=fixture / "hj_literature_route_a.tpath")
    parser.add_argument("--manifest", type=Path, default=fixture / "source_manifest_hj.json")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--tolerance-A", type=float, default=1.0e-9)
    args = parser.parse_args(argv)

    try:
        pdb = args.pdb.resolve()
        ligand_itp = args.ligand_itp.resolve()
        route = args.route.resolve()
        manifest = json.loads(args.manifest.resolve().read_text(encoding="utf-8"))
        if sha256(pdb) != manifest["reference_structure_sha256"]:
            raise RuntimeError("reference PDB hash does not match H/J manifest")
        if sha256(route) != manifest["route_sha256"]:
            raise RuntimeError("route hash does not match H/J manifest")
        ca, ligand = parse_pdb(pdb)
        masses = parse_ligand_masses(ligand_itp)
        total_mass = sum(masses[name] for name, _ in ligand)
        start = tuple(
            sum(masses[name] * position[axis] for name, position in ligand)
            / total_mass
            for axis in range(3)
        )
        inner = mean([ca[residue] for residue in INNER_RESIDUES])
        outer = mean([ca[residue] for residue in OUTER_RESIDUES])
        middle = scale(add(inner, outer), 0.5)
        outward = normalize(subtract(outer, start))
        derived = [
            start,
            middle,
            outer,
            add(outer, scale(outward, EXTENSIONS_A[0])),
            add(outer, scale(outward, EXTENSIONS_A[1])),
        ]
        accepted = route_points(route)
        errors = [
            math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))
            for a, b in zip(derived, accepted)
        ]
        max_error = max(errors)
        passed = max_error <= args.tolerance_A
        report = {
            "schema": "titan-path-4w52-hj-route-derivation-report-v1",
            "status": "PASS" if passed else "FAIL",
            "reference_pdb_sha256": sha256(pdb),
            "ligand_itp_sha256": sha256(ligand_itp),
            "route_sha256": sha256(route),
            "inner_gate_residues": list(INNER_RESIDUES),
            "outer_gate_residues": list(OUTER_RESIDUES),
            "terminal_extensions_A": list(EXTENSIONS_A),
            "derived_waypoints_A": derived,
            "accepted_waypoints_A": accepted,
            "waypoint_error_A": errors,
            "max_waypoint_error_A": max_error,
            "tolerance_A": args.tolerance_A,
        }
        payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if args.report is not None:
            output = args.report.resolve()
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(payload, encoding="utf-8", newline="\n")
        print(payload, end="")
        return 0 if passed else 1
    except (OSError, RuntimeError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"TITAN_PATH_HJ_ROUTE_DERIVATION_ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
