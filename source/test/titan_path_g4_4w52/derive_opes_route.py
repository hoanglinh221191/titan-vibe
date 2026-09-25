#!/usr/bin/env python3
"""Reproduce the 4W52 pilot centerline from pinned public OPES structures."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


SELECTED_STRUCTURES = (0, 117, 128, 131, 129, 130)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_index(path: Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    active: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and "]" in line:
            name = line[1 : line.index("]")].strip()
            # The public index contains a duplicate BNZ header.  Bind the first
            # definition and reject neither silently nor by concatenating it.
            active = name if name not in groups else None
            groups.setdefault(name, [])
        elif active is not None:
            groups[active].extend(int(token) - 1 for token in line.split())
    for required in ("C-alpha", "BNZ"):
        if required not in groups or not groups[required]:
            raise RuntimeError(f"index group {required!r} is missing from {path}")
    return groups


def read_gro_coordinates(path: Path) -> np.ndarray:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:
        raise RuntimeError(f"truncated GRO file: {path}")
    atom_count = int(lines[1].strip())
    atom_lines = lines[2 : 2 + atom_count]
    if len(atom_lines) != atom_count:
        raise RuntimeError(f"GRO atom count mismatch in {path}")
    # GRO positions are in nm; TITAN route coordinates are in angstrom.
    return np.asarray(
        [
            (float(line[20:28]), float(line[28:36]), float(line[36:44]))
            for line in atom_lines
        ],
        dtype=np.float64,
    ) * 10.0


def make_molecules_whole(
    *, gmx: str, tpr: Path, source: Path, destination: Path
) -> None:
    completed = subprocess.run(
        [
            gmx,
            "trjconv",
            "-s",
            str(tpr),
            "-f",
            str(source),
            "-o",
            str(destination),
            "-pbc",
            "mol",
        ],
        input="0\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode != 0 or not destination.is_file():
        raise RuntimeError(
            f"gmx trjconv failed for {source} with exit code "
            f"{completed.returncode}:\n{completed.stdout}"
        )


def align_to_reference(
    mobile: np.ndarray, reference: np.ndarray, fit_indices: list[int]
) -> tuple[np.ndarray, float]:
    mobile_fit = mobile[fit_indices]
    reference_fit = reference[fit_indices]
    mobile_center = mobile_fit.mean(axis=0)
    reference_center = reference_fit.mean(axis=0)
    covariance = (mobile_fit - mobile_center).T @ (
        reference_fit - reference_center
    )
    left, _, right_t = np.linalg.svd(covariance)
    rotation = right_t.T @ left.T
    if np.linalg.det(rotation) < 0.0:
        right_t[-1, :] *= -1.0
        rotation = right_t.T @ left.T
    aligned = (mobile - mobile_center) @ rotation.T + reference_center
    residual = aligned[fit_indices] - reference_fit
    rmsd = float(np.sqrt(np.mean(np.sum(residual * residual, axis=1))))
    return aligned, rmsd


def route_waypoints(path: Path) -> np.ndarray:
    points: list[tuple[float, float, float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        fields = raw.strip().split()
        if fields and fields[0] == "waypoint" and len(fields) == 7:
            points.append((float(fields[3]), float(fields[4]), float(fields[5])))
    if len(points) != len(SELECTED_STRUCTURES):
        raise RuntimeError(
            f"route has {len(points)} waypoints; expected {len(SELECTED_STRUCTURES)}"
        )
    return np.asarray(points, dtype=np.float64)


def validate_manifest_hashes(
    manifest: dict[str, object], template: Path, structures: Path
) -> dict[str, str]:
    observed: dict[str, str] = {}
    reference_files = manifest["reference_files"]
    route_inputs = manifest["route_input_files"]
    if not isinstance(reference_files, dict) or not isinstance(route_inputs, dict):
        raise RuntimeError("manifest hash sections must be JSON objects")
    files = {
        "prot_lig.pdb_sha256": template / "prot_lig.pdb",
        "lysozyme.gro_sha256": template / "lysozyme.gro",
        "lysozyme.top_sha256": template / "lysozyme.top",
        **{
            f"eq_structure{index}.gro_sha256": structures
            / f"eq_structure{index}.gro"
            for index in SELECTED_STRUCTURES
        },
    }
    for key, path in files.items():
        expected_section = reference_files if key in reference_files else route_inputs
        expected = str(expected_section.get(key, "")).lower()
        actual = sha256(path)
        observed[key] = actual
        if actual != expected:
            raise RuntimeError(
                f"source hash mismatch for {path}: expected {expected}, observed {actual}"
            )
    return observed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    fixture = Path(__file__).resolve().parent
    parser.add_argument(
        "--repository",
        type=Path,
        default=root / ".codex_tmp" / "DTW_Path_Classification",
    )
    parser.add_argument("--route", type=Path, default=fixture / "opes_route_a.tpath")
    parser.add_argument(
        "--manifest", type=Path, default=fixture / "source_manifest.json"
    )
    parser.add_argument("--report", type=Path)
    parser.add_argument("--tolerance-A", type=float, default=0.02)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument(
        "--tpr",
        type=Path,
        default=root
        / "docs"
        / "titan_path"
        / "gates"
        / "G4_4w52_pilot_20260825"
        / "run_001"
        / "initial.tpr",
    )
    args = parser.parse_args(argv)

    try:
        repository = args.repository.resolve()
        template = repository / "T4L_benzene" / "opes_explore" / "template"
        structures = repository / "T4L_benzene" / "opes_flooding" / "structures"
        manifest = json.loads(args.manifest.resolve().read_text(encoding="utf-8"))
        observed_hashes = validate_manifest_hashes(manifest, template, structures)
        groups = read_index(template / "index.ndx")
        tpr = args.tpr.resolve()
        if not tpr.is_file():
            raise RuntimeError(f"whole-molecule conversion TPR is absent: {tpr}")
        with tempfile.TemporaryDirectory(prefix="titan_path_4w52_") as temporary:
            whole_dir = Path(temporary)
            whole_files: dict[int, Path] = {}
            for index in SELECTED_STRUCTURES:
                whole = whole_dir / f"eq_structure{index}_whole.gro"
                make_molecules_whole(
                    gmx=args.gmx,
                    tpr=tpr,
                    source=structures / f"eq_structure{index}.gro",
                    destination=whole,
                )
                whole_files[index] = whole
            reference = read_gro_coordinates(whole_files[0])
            waypoints = []
            fit_rmsd = []
            for index in SELECTED_STRUCTURES:
                coordinates = read_gro_coordinates(whole_files[index])
                if coordinates.shape != reference.shape:
                    raise RuntimeError(f"atom-count mismatch in structure {index}")
                aligned, rmsd = align_to_reference(
                    coordinates, reference, groups["C-alpha"]
                )
                # Benzene is symmetric; the arithmetic all-atom center is the
                # published route center used for this pilot.
                waypoints.append(aligned[groups["BNZ"]].mean(axis=0))
                fit_rmsd.append(rmsd)
        derived = np.asarray(waypoints)
        accepted = route_waypoints(args.route.resolve())
        errors = np.linalg.norm(derived - accepted, axis=1)
        max_error = float(errors.max())
        passed = max_error <= args.tolerance_A
        report = {
            "schema": "titan-path-4w52-route-derivation-report-v1",
            "status": "PASS" if passed else "FAIL",
            "selected_structures": list(SELECTED_STRUCTURES),
            "alignment_group": "C-alpha",
            "ligand_group": "BNZ",
            "source_sha256": observed_hashes,
            "whole_molecule_tpr_sha256": sha256(tpr),
            "route_sha256": sha256(args.route.resolve()),
            "derived_waypoints_A": derived.tolist(),
            "accepted_waypoints_A": accepted.tolist(),
            "waypoint_error_A": errors.tolist(),
            "max_waypoint_error_A": max_error,
            "tolerance_A": args.tolerance_A,
            "alignment_rmsd_A": fit_rmsd,
        }
        payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if args.report is not None:
            output = args.report.resolve()
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(payload, encoding="utf-8", newline="\n")
        print(payload, end="")
        return 0 if passed else 1
    except (OSError, RuntimeError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"TITAN_PATH_ROUTE_DERIVATION_ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
