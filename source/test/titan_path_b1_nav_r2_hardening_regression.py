#!/usr/bin/env python3
"""Optional arithmetic/controller regression against the exact real B1 R2 bundle.

This consumes archived decisions. It does not re-select candidates, evaluate
clearance/energies on molecular coordinates, or execute MD. No synthetic bundle
can substitute for the frozen SHA256 identity below.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
from typing import Any

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
EXPECTED_SHA256 = "a06b63e903fb53b86734e857b9b248143104f4b1a55051c5f84f0409c052bc8d"
P2_SHA256 = "36a17877db91381af9d4ecaaf4ce8d4f837553338c8438c7dcaa965cafaa2d03"
CELLS = tuple(f"branch_{branch}__{route}" for branch in ("a", "b")
              for route in ("blue", "cyan", "orange", "pink"))
CONTINUE_CELLS = {"branch_a__blue", "branch_a__cyan", "branch_b__cyan",
                  "branch_b__orange", "branch_b__pink"}
STATUS_VECTOR_TOLERANCE_A = 1.0e-7
FILE_VECTOR_TOLERANCE_A = 5.01e-5
ANGLE_TOLERANCE_DEG = 1.0e-6


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class Archive:
    def __init__(self, handle: tarfile.TarFile):
        self.handle = handle
        self.files: dict[str, tarfile.TarInfo] = {}
        seen: set[str] = set()
        for member in handle.getmembers():
            name = member.name
            parts = PurePosixPath(name).parts
            require(not name.startswith("/") and "\\" not in name and
                    all(part not in ("..", ".") for part in parts),
                    f"unsafe archive member {name!r}")
            require(name not in seen, f"duplicate archive member {name}")
            seen.add(name)
            if member.isdir():
                continue
            require(member.isfile(), f"non-regular archive member {name}")
            self.files[name] = member

    def read(self, name: str) -> bytes:
        require(name in self.files, f"archive member missing: {name}")
        stream = self.handle.extractfile(self.files[name])
        require(stream is not None, f"cannot read archive member: {name}")
        assert stream is not None
        return stream.read()

    def json(self, name: str) -> Any:
        return json.loads(self.read(name).decode("utf-8"))


def status_values(payload: bytes) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in payload.decode("utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        key, value = line.split("=", 1)
        require(key not in result, f"duplicate status field {key}")
        result[key] = value
    return result


def finite(status: dict[str, str], key: str) -> float:
    value = float(status[key])
    require(math.isfinite(value), f"non-finite status field {key}")
    return value


def vector(status: dict[str, str], prefix: str) -> list[float]:
    return [finite(status, f"{prefix}_{axis}_A") for axis in "XYZ"]


def parse_vector_file(payload: bytes) -> list[float]:
    first = payload.decode("utf-8").splitlines()[0].split()
    require(len(first) >= 3, "vector file has fewer than three components")
    value = [float(token) for token in first[:3]]
    require(all(math.isfinite(item) for item in value), "non-finite vector file")
    require(math.hypot(*value) > 0.0, "zero vector file")
    return value


def source_hashes() -> dict[str, str]:
    names = ["tools/titan_path_controller.py", "tools/titan_path_handoff.py",
             "tools/titan_path_native_contract.py", "src/titan_path/control.cpp",
             "src/titan_path/control.h", "test/titan_path_b1_nav_r2_turn_probe.cpp",
             "test/titan_path_b1_nav_r2_hardening_regression.py"]
    return {name: sha256_file(ROOT / name) for name in names}


def run(args: argparse.Namespace) -> dict[str, Any]:
    bundle = args.bundle.resolve(strict=True)
    require(bundle.is_file(), "bundle must be a regular file")
    observed_hash = sha256_file(bundle)
    require(observed_hash == EXPECTED_SHA256,
            f"R2 bundle SHA256 mismatch: expected {EXPECTED_SHA256}, observed {observed_hash}")

    before_sources = source_hashes()
    sys.path.insert(0, str(ROOT / "tools"))
    spec = importlib.util.spec_from_file_location(
        "b1_r2_current_controller", ROOT / "tools" / "titan_path_controller.py")
    require(spec is not None and spec.loader is not None, "cannot load current controller")
    assert spec is not None and spec.loader is not None
    controller = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(controller)

    output = args.output_dir.resolve()
    require(not output.exists(), f"output directory already exists: {output}")
    output.mkdir(parents=True, exist_ok=False)
    records: dict[str, Any] = {}
    input_rows: list[str] = []
    config = {
        "reference_displacement_contract": "match_titan_probe",
        "allowed_probe_lengths_A": [1.0, 0.5, 0.25],
        "rate_nm_per_ps": 0.005, "dt_ps": 0.002, "segment_steps": 10000,
    }
    with tarfile.open(bundle, "r:gz") as top_handle:
        top = Archive(top_handle)
        replay = top.json("generated/NAV_P2_REPLAY_REPORT.json")
        binding = top.json("generated/NAV_P2_INPUT_BINDING.json")
        require(replay.get("md_executed") is False and binding.get("md_executed") is False,
                "R2 evidence must state md_executed=false")
        require(set(replay["records"]) == set(CELLS), "unexpected R2 record set")
        bound = {entry["cell_id"]: entry for entry in binding["cells"]}
        require(len(binding["cells"]) == len(bound) == 8 and set(bound) == set(CELLS),
                "unexpected or duplicated binding cell set")
        p2_payload = top.read("baselines/P2_R2_RESULT.tar.gz")
        require(sha256_bytes(p2_payload) == P2_SHA256,
                "embedded P2 parameter/vector baseline hash mismatch")
        with tarfile.open(fileobj=io.BytesIO(p2_payload), mode="r:gz") as p2_handle:
            p2 = Archive(p2_handle)
            for cell in CELLS:
                prefix = f"outputs/nav_p2_fixed_snapshot_replays/{cell}/"
                status_payload = top.read(prefix + "pull_direction_status.txt")
                status_hash = sha256_bytes(status_payload)
                require(status_hash == replay["records"][cell]["status_sha256"],
                        f"{cell}: status hash/report mismatch")
                status = status_values(status_payload)
                require(status["CLEARANCE_MODE"] == "swept-linear-center",
                        f"{cell}: unexpected clearance mode")
                require(finite(status, "MIN_CLEARANCE_THRESHOLD_A") == 1.5,
                        f"{cell}: unexpected hard threshold")
                old_payload = p2.read(
                    f"outputs/topology_complete_terminal_replays/{cell}/DE_out_vector.dat")
                require(sha256_bytes(old_payload) == bound[cell]["old_vector_sha256"],
                        f"{cell}: old vector does not match frozen binding")
                old_vector = parse_vector_file(old_payload)
                record: dict[str, Any] = {
                    "archived_status": status["STATUS"], "archived_mode": status["MODE"],
                    "archived_reason_detail": status["REASON_DETAIL"],
                    "status_sha256": status_hash,
                    "old_vector_sha256": sha256_bytes(old_payload),
                }
                if cell not in CONTINUE_CELLS:
                    require(status["STATUS"] == "STOP", f"{cell}: expected frozen STOP")
                    stop_kind = controller.validate_route_stop_diagnostics(status)
                    require(stop_kind == "FORWARD_PROGRESS_EXHAUSTED",
                            f"{cell}: unexpected STOP diagnostics result {stop_kind}")
                    record["controller_stop_diagnostics_pass"] = True
                    records[cell] = record
                    continue

                require(status["STATUS"] == "CONTINUE", f"{cell}: expected frozen CONTINUE")
                controller.validate_continue_route_status(
                    status, expected_clearance_mode="swept-linear-center")
                probe_contract = controller.validate_route_probe_contract(status, config)
                requested = vector(status, "ROUTE_TARGET_VECTOR")
                require(abs(math.hypot(*requested) - finite(status, "REQUESTED_PULL_DISTANCE_A"))
                        <= STATUS_VECTOR_TOLERANCE_A,
                        f"{cell}: ROUTE_TARGET_VECTOR is not the requested displacement")
                max_turn = finite(status, "MAX_TURN_ANGLE_DEG")
                probe = finite(status, "NEXT_REFERENCE_DISPLACEMENT_A")
                input_rows.append(cell + "\t" + "\t".join(
                    format(value, ".17g") for value in requested + old_vector + [max_turn, probe]))
                output_payload = top.read(prefix + "DE_out_vector.dat")
                require(sha256_bytes(output_payload) == replay["records"][cell]["vector_sha256"],
                        f"{cell}: output vector hash/report mismatch")
                record.update({
                    "controller_continue_validation_pass": True,
                    "probe_contract": probe_contract,
                    "requested_vector_A": requested,
                    "old_vector": old_vector,
                    "archived_next_vector_A": vector(status, "NEXT_VECTOR"),
                    "archived_vector_file_A": parse_vector_file(output_payload),
                    "ligand_com_A": vector(status, "LIGAND_COM"),
                    "archived_requested_angle_deg": finite(status, "ROUTE_TARGET_ANGLE_DEG"),
                    "archived_applied_angle_deg": finite(status, "NEXT_VECTOR_ANGLE_DEG"),
                    "archived_turn_capped": status["TURN_CAPPED"] == "1",
                    "selected_probe_A": probe,
                })
                records[cell] = record

    compiler = shutil.which(args.cxx)
    require(compiler is not None, f"C++ compiler not found: {args.cxx}")
    helper = output / ("turn_probe.exe" if os.name == "nt" else "turn_probe")
    command = [compiler, "-std=c++20", "-O2", "-I" + str(ROOT / "src"),
               str(ROOT / "test" / "titan_path_b1_nav_r2_turn_probe.cpp"),
               str(ROOT / "src" / "titan_path" / "control.cpp"), "-o", str(helper)]
    built = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                           encoding="utf-8", errors="replace", check=False)
    (output / "compile.log").write_text(built.stdout + built.stderr, encoding="utf-8")
    require(built.returncode == 0, "C++ arithmetic helper compilation failed; see compile.log")
    inputs = "\n".join(input_rows) + "\n"
    (output / "turn_probe_input.tsv").write_text(inputs, encoding="utf-8")
    completed = subprocess.run([str(helper)], input=inputs, capture_output=True,
                               text=True, encoding="utf-8", errors="strict", check=False)
    (output / "turn_probe_output.tsv").write_text(completed.stdout, encoding="utf-8")
    require(completed.returncode == 0, "C++ turn/probe evaluation failed: " + completed.stderr)
    processed: set[str] = set()
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        require(len(fields) == 8, "unexpected C++ helper output schema")
        cell = fields[0]
        require(cell in CONTINUE_CELLS and cell not in processed,
                "unexpected or duplicated C++ helper output cell")
        processed.add(cell)
        record = records[cell]
        requested_angle, applied_angle = map(float, fields[2:4])
        next_vector = list(map(float, fields[4:7]))
        length = float(fields[7])
        require(all(math.isfinite(value) for value in
                    [requested_angle, applied_angle, length] + next_vector),
                f"{cell}: non-finite C++ output")
        require((fields[1] == "1") == record["archived_turn_capped"],
                f"{cell}: turn-cap decision changed")
        status_error = max(abs(a - b) for a, b in zip(next_vector, record["archived_next_vector_A"]))
        file_error = max(abs(a - b) for a, b in zip(next_vector, record["archived_vector_file_A"]))
        require(status_error <= STATUS_VECTOR_TOLERANCE_A,
                f"{cell}: current turn/probe vector differs from archived status: {status_error}")
        require(file_error <= FILE_VECTOR_TOLERANCE_A,
                f"{cell}: current turn/probe vector differs from archived 4dp file: {file_error}")
        require(abs(length - record["selected_probe_A"]) <= 1.0e-12,
                f"{cell}: current probe scaling length mismatch")
        require(abs(requested_angle - record["archived_requested_angle_deg"]) <= ANGLE_TOLERANCE_DEG
                and abs(applied_angle - record["archived_applied_angle_deg"]) <= ANGLE_TOLERANCE_DEG,
                f"{cell}: current turn angles differ from archived telemetry")
        record.update({"current_turn_probe_pass": True, "current_next_vector_A": next_vector,
                       "status_vector_max_abs_error_A": status_error,
                       "file_vector_max_abs_error_A": file_error,
                       "current_endpoint_A": [a + b for a, b in zip(record["ligand_com_A"], next_vector)]})
    require(processed == CONTINUE_CELLS, "not all five selected vectors were recomputed")
    require(source_hashes() == before_sources, "source changed during regression; use a new output directory")
    report = {
        "status": "PASS_REAL_B1_R2_ARCHIVED_CONTRACT_AND_TURN_PROBE_REGRESSION",
        "evidence_kind": "HASH_BOUND_REAL_FIXED_SNAPSHOT_REPLAY_OUTPUTS",
        "bundle_sha256": observed_hash, "p2_baseline_sha256": P2_SHA256,
        "controller_status_pass_count": 8, "stop_diagnostics_pass_count": 3,
        "probe_contract_pass_count": 5, "turn_probe_numeric_match_count": 5,
        "candidate_reselection_executed": False, "molecular_clearance_recomputed": False,
        "md_executed": False, "synthetic_input_allowed": False,
        "status_vector_tolerance_A": STATUS_VECTOR_TOLERANCE_A,
        "file_vector_tolerance_A": FILE_VECTOR_TOLERANCE_A,
        "angle_tolerance_deg": ANGLE_TOLERANCE_DEG,
        "source_sha256": before_sources, "helper_sha256": sha256_file(helper),
        "compiler_command": command, "probe_contract_test_config": config,
        "records": records,
    }
    destination = output / "B1_NAV_R2_HARDENING_REGRESSION.json"
    destination.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    destination.with_suffix(".json.sha256").write_text(
        sha256_file(destination) + "  " + destination.name + "\n", encoding="ascii")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path,
                        help="Exact downloaded B1 NAV-P2 R2 analysis bundle; SHA256 is fixed in this test")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="New local regression output directory; never modifies downloaded evidence")
    parser.add_argument("--cxx", default="clang++")
    args = parser.parse_args()
    try:
        result = run(args)
        print("PASS B1_NAV_R2_HARDENING_REGRESSION statuses=8/8 probes=5/5 "
              "turn_probe_matches=5/5 candidate_reselection=false md_executed=false "
              f"bundle_sha256={result['bundle_sha256']}")
        print(f"REPORT={args.output_dir.resolve() / 'B1_NAV_R2_HARDENING_REGRESSION.json'}")
        return 0
    except (OSError, RuntimeError, ValueError, KeyError, tarfile.TarError) as error:
        print(f"FAIL B1_NAV_R2_HARDENING_REGRESSION: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
