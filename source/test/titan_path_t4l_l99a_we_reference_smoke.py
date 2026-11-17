#!/usr/bin/env python3
"""Focused regression gate for the frozen T4L L99A WE reference dossier."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
MODULE_DIR = ROOT / "tools" / "titan_path_reference"
if not MODULE_DIR.is_dir():
    MODULE_DIR = ROOT / "tools"
sys.path.insert(0, str(MODULE_DIR))

import build_t4l_l99a_we_reference as builder  # noqa: E402
from evaluate_reference_path import evaluate  # noqa: E402
from reference_metrics import classify_reference, validate_polyline  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_alignment() -> None:
    source = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 3.0),
    ]
    target = [(5.0 - y, -2.0 + x, 7.0 + z) for x, y, z in source]
    transform = builder.fit_rigid_transform(source, target)
    observed = [builder.transform_point(transform, point) for point in source]
    require(transform.rmsd_angstrom < 1.0e-12, "synthetic rigid RMSD drift")
    require(
        max(math.dist(left, right) for left, right in zip(observed, target))
        < 1.0e-12,
        "synthetic rigid transform drift",
    )


def test_dossier(dossier_path: Path) -> None:
    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    require(
        dossier.get("schema") == "titan-path-t4l-l99a-we-reference-dossier-v1",
        "dossier schema mismatch",
    )
    require(
        dossier.get("status") == "REFERENCE_DOSSIER_BUILT_NO_TITAN_OUTCOMES",
        "dossier status mismatch",
    )
    require(dossier["system"]["reference_pdb_id"] == "1L83", "wrong PDB frame")
    require(
        dossier["we_source"]["article_doi"] == "10.1016/j.bpj.2018.01.014",
        "correct WE DOI is absent",
    )
    legacy = dossier["legacy_provenance_correction"]
    require(legacy["incorrect_doi_present"] is True, "legacy DOI mismatch not retained")
    require(legacy["correct_doi_present"] is False, "legacy DOI unexpectedly corrected")
    require(
        legacy["g4c3_construct_differences_from_we_reference"]
        == [{"residue_number": 153, "reference": "PHE", "observed": "ALA"}],
        "G4C3 F153A mismatch changed",
    )
    require(
        [item["residue_number"] for item in legacy["rcsb_4w52_construct_differences_from_we_reference"]]
        == [54, 97, 163, 164],
        "RCSB 4W52 construct audit changed",
    )

    routes = dossier["routes"]
    require(set(routes) == {"blue", "orange", "pink", "cyan"}, "route class set mismatch")
    require(
        {label: routes[label]["published_model_count"] for label in routes}
        == {"blue": 25, "orange": 23, "pink": 24, "cyan": 26},
        "published Voronoi model counts changed",
    )
    require(
        sum(route["dense_holdout_count"] for route in routes.values()) == 78,
        "dense holdout count mismatch",
    )
    require(
        max(route["alignment_ca_rmsd_max_angstrom"] for route in routes.values()) < 3.1,
        "1L83 CA alignment RMSD exceeds frozen pre-outcome gate",
    )
    require(
        max(route["reference_bound_start_offset_angstrom"] for route in routes.values()) < 1.1,
        "reference bound-center offset exceeds frozen pre-outcome gate",
    )

    references = {
        label: validate_polyline(route["dense_centerline_angstrom"], label)
        for label, route in routes.items()
    }
    for label, route in routes.items():
        result = classify_reference(references[label], references)
        require(result["assigned_class"] == label, f"self-classification failed: {label}")
        require(result["dtw_margin_angstrom"] > 1.0, f"reference margin collapsed: {label}")
        sparse = [
            references[label][index]
            for index in route["sparse_landmark_indices_0based"]
        ]
        sparse_result = classify_reference(sparse, references)
        require(
            sparse_result["assigned_class"] == label,
            f"sparse landmark classification failed: {label}",
        )
        landmark_path = dossier_path.parent / route["sparse_landmark_file"]
        require(landmark_path.is_file(), f"sparse landmark file absent: {label}")
        require(sha256(landmark_path) == route["sparse_landmark_sha256"], "landmark hash drift")
        landmark = json.loads(landmark_path.read_text(encoding="utf-8"))
        require(landmark["executable_titan_route"] is False, "unbound landmarks became executable")

    dense_csv = dossier_path.parent / dossier["dense_centerline_csv"]
    require(dense_csv.is_file(), "dense centerline CSV absent")
    require(sha256(dense_csv) == dossier["dense_centerline_csv_sha256"], "dense CSV hash drift")

    with tempfile.TemporaryDirectory(prefix="titan_path_reference_smoke_") as temp:
        candidate = Path(temp) / "cyan_candidate.csv"
        with candidate.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("x_angstrom", "y_angstrom", "z_angstrom"),
                lineterminator="\n",
            )
            writer.writeheader()
            for x, y, z in references["cyan"]:
                writer.writerow(
                    {
                        "x_angstrom": f"{x:.9f}",
                        "y_angstrom": f"{y:.9f}",
                        "z_angstrom": f"{z:.9f}",
                    }
                )
        accepted = evaluate(dossier_path, candidate, "cyan")
        rejected = evaluate(dossier_path, candidate, "orange")
        require(accepted["assigned_class"] == "cyan", "CSV evaluator class drift")
        require(accepted["reference_consistent"] is True, "matching prescription rejected")
        require(rejected["reference_consistent"] is False, "wrong prescription accepted")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dossier", type=Path, required=True)
    args = parser.parse_args()
    try:
        test_alignment()
        test_dossier(args.dossier.resolve())
        print(
            "T4L_L99A_WE_REFERENCE_SMOKE_PASS "
            "classes=4 dense_points=98 holdout_points=78 executable_routes=0"
        )
        return 0
    except (OSError, ValueError, KeyError, RuntimeError, json.JSONDecodeError) as error:
        print(f"T4L_L99A_WE_REFERENCE_SMOKE_ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
