#!/usr/bin/env python3
"""Focused tests for the metadata-only G5 preflight contract."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
TOOL = REPO / "tools" / "g5_readonly_preflight.py"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(role: str, path: Path, policy: str = "always") -> dict[str, object]:
    return {
        "artifact_role": role,
        "path": str(path),
        "expected_bytes": path.stat().st_size,
        "expected_sha256": digest(path),
        "hash_policy": policy,
    }


def run_tool(spec: Path, output: Path, table: Path, cap: int = 8) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(TOOL),
            "--spec",
            str(spec),
            "--output",
            str(output),
            "--table",
            str(table),
            "--max-hash-bytes",
            str(cap),
        ],
        check=False,
        capture_output=True,
        text=True,
    )


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="g5_readonly_preflight_") as tmp:
        root = Path(tmp)
        protocol = root / "protocol.md"
        protocol.write_text("current protocol\n", encoding="utf-8")
        metadata = root / "status.json"
        metadata.write_text('{"frames": 7, "outcomes": false}\n', encoding="utf-8")
        large = root / "trajectory.xtc"
        large.write_bytes(b"0123456789abcdef")
        legacy = root / "legacy_sha256.txt"
        legacy.write_text(f"{'0' * 64}  {protocol}\n", encoding="utf-8")

        spec_payload = {
            "schema_version": 1,
            "protocol_id": "g5-smoke-v1",
            "mode": "metadata_only",
            "outcome_execution_authorized": False,
            "forbidden_path_fragments": ["outcome_table.csv"],
            "datasets": [
                {
                    "dataset_id": "reference",
                    "data_role": "development_reference",
                    "known_outcome_exposure": True,
                    "confirmatory_eligible": False,
                    "development_replay_candidate": False,
                    "trajectory_replay_authorized": False,
                    "expected_preflight_status": "REFERENCE_ONLY",
                    "bindings": [binding("metadata", metadata)],
                    "json_assertions": [
                        {"path": str(metadata), "pointer": "/frames", "expected": 7}
                    ],
                },
                {
                    "dataset_id": "lineage_blocked",
                    "data_role": "development_descriptive",
                    "known_outcome_exposure": True,
                    "confirmatory_eligible": False,
                    "development_replay_candidate": False,
                    "trajectory_replay_authorized": False,
                    "expected_preflight_status": "BLOCKED_PROTOCOL_LINEAGE",
                    "bindings": [
                        binding("protocol", protocol),
                        binding("legacy_manifest", legacy),
                    ],
                    "lineage_checks": [
                        {"manifest_path": str(legacy), "target_path": str(protocol)}
                    ],
                },
                {
                    "dataset_id": "large_deferred",
                    "data_role": "development_trajectory",
                    "known_outcome_exposure": True,
                    "confirmatory_eligible": False,
                    "development_replay_candidate": True,
                    "trajectory_replay_authorized": False,
                    "expected_preflight_status": "NEEDS_FULL_HASH_FOR_REPLAY",
                    "bindings": [binding("trajectory", large, "defer_if_large")],
                },
            ],
        }
        spec = root / "spec.json"
        spec.write_text(json.dumps(spec_payload), encoding="utf-8")
        before = {
            path: (path.stat().st_size, path.stat().st_mtime_ns, digest(path))
            for path in (protocol, metadata, large, legacy, spec)
        }
        result = run_tool(spec, root / "report.json", root / "report.tsv")
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        report = json.loads((root / "report.json").read_text(encoding="utf-8"))
        if report["status"] != "PASS" or report["outcome_execution_authorized"]:
            raise RuntimeError("valid metadata-only fixture did not pass closed")
        after = {
            path: (path.stat().st_size, path.stat().st_mtime_ns, digest(path))
            for path in before
        }
        if before != after:
            raise RuntimeError("preflight modified an allowlisted input")

        bad_spec = root / "bad_spec.json"
        bad_payload = json.loads(json.dumps(spec_payload))
        bad_payload["datasets"][0]["bindings"][0]["path"] = str(
            root / "outcome_table.csv"
        )
        bad_spec.write_text(json.dumps(bad_payload), encoding="utf-8")
        bad = run_tool(bad_spec, root / "bad.json", root / "bad.tsv")
        if bad.returncode == 0 or "outcome-like path is forbidden" not in bad.stderr:
            raise RuntimeError("preflight accepted an outcome-like binding")

        mismatch_spec = root / "mismatch_spec.json"
        mismatch_payload = json.loads(json.dumps(spec_payload))
        mismatch_payload["datasets"][2]["expected_preflight_status"] = (
            "READY_DEVELOPMENT"
        )
        mismatch_spec.write_text(json.dumps(mismatch_payload), encoding="utf-8")
        mismatch = run_tool(
            mismatch_spec, root / "mismatch.json", root / "mismatch.tsv"
        )
        if mismatch.returncode == 0:
            raise RuntimeError("preflight accepted an unexpected readiness status")

    print(
        "PASS g5_readonly_preflight_smoke "
        "read_only=PASS lineage=BLOCKED large_hash=DEFERRED outcomes=REJECT"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
