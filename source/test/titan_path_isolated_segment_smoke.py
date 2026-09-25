#!/usr/bin/env python3
"""Fast synthetic tests; never invoke native MD."""
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import titan_path_isolated_segment as isolated
import titan_path_controller as c
import titan_path_controller_native_smoke as native_fixture


class IsolatedSmoke(unittest.TestCase):
    def setUp(self):
        self.fixture = native_fixture.NativeControllerSmoke()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root

    def test_annotated_archived_vector_is_preserved_and_checked(self):
        path = self.root / "vector.dat"
        path.write_text("0.5000 0.0000 0.0000 Angle:0.0000 -- Modul:0.5000 Angstrom\n# retained\n")
        self.fixture.rebind()
        c.validate_initial_nav(self.fixture.load())
        for text in (".5 0 0 Angle:0 -- Modul:nan Angstrom", ".5 0 0 Angle:20 -- Modul:.5 Angstrom", ".5 0 0 ignored suffix"):
            path.write_text(text)
            with self.assertRaises(c.ControllerError):
                c.parse_navigation_vector_file(path, self.fixture.status)

    def test_native_time_coverage_not_just_exit_zero(self):
        good = "Reading frame 0 time 20.000\nLast frame 50 time 30.000\n"
        self.assertEqual(isolated.check_time_coverage(good, 20, 30, complete=True)["frames"], 51)
        for bad in (good.replace("30.000", "29.000"), good.replace("20.000", "21.000"), "exit code zero"):
            with self.assertRaises(c.ControllerError):
                isolated.check_time_coverage(bad, 20, 30, complete=True)

    def test_actual_native_part_names_are_discovered(self):
        directory = self.root / "invocation"
        directory.mkdir()
        names = ["run.part0007.log", "pullx.part0007.xvg", "pullf.part0007.xvg", "run.part0007.edr", "run.part0007.xtc"]
        for name in [*names, "state.cpt", "final.part0007.gro"]:
            (directory / name).write_text("SYNTHETIC native output")
        observed = {"simulation_part": 7, "output_files": [str(directory / name) for name in names]}
        output = isolated.discover_output_map(directory, observed, require_gro=True)
        self.assertEqual(Path(output["gro"]).name, "final.part0007.gro")
        with self.assertRaises(c.ControllerError):
            isolated.discover_output_map(directory, {**observed, "simulation_part": 8}, require_gro=True)
        (directory / "final.gro").write_text("ambiguous")
        with self.assertRaises(c.ControllerError):
            isolated.discover_output_map(directory, observed, require_gro=True)

    def test_r1_resume_refuses_before_any_command(self):
        with patch.object(isolated, "command", side_effect=AssertionError("must not execute")), self.assertRaisesRegex(c.ControllerError, "does not support resume"):
            isolated.run_segment({}, resume=True)

    def test_existing_complete_or_partial_attempt_never_repeats(self):
        work = self.root / "isolated"
        work.mkdir()
        state = work / "isolated_state.json"
        for status in ("COMPLETE", "INTERRUPTED"):
            state.write_text(json.dumps({"status": status}))
            before = state.read_bytes()
            with patch.object(isolated, "preflight", side_effect=AssertionError("must not re-execute")), self.assertRaisesRegex(c.ControllerError, "no automatic retry"):
                isolated.run_segment({"work_dir": work})
            self.assertEqual(state.read_bytes(), before)

    def test_partial_native_endpoint_never_becomes_success(self):
        work = self.root / "isolated"
        directory = work / "invocations/invocation_0000"
        directory.mkdir(parents=True)
        expected = {"natoms": 3, "init_step": 10000, "nsteps": 5000, "dt_ps": .002}
        partial = {"natoms": 3, "step": 12500, "time_ps": 25, "simulation_part": 3}
        with patch.object(isolated, "dump_checkpoint", return_value=partial), self.assertRaisesRegex(c.ControllerError, "separate recovery"):
            isolated.inspect_invocation({"work_dir": work}, directory, expected, 10000, -.1, complete=True, seed_part=2)

    def test_source_hash_drift_rejected(self):
        source = self.root / "original.cpt"
        source.write_text("original")
        binding = self.root / "source_binding.json"
        binding.write_text(json.dumps({"source_files_sha256": {str(source): c.sha256_file(source)}}))
        config = {"initial_nav": {"binding_file": binding}, "work_dir": self.root / "output"}
        isolated.source_hashes(config)
        source.write_text("changed")
        with self.assertRaisesRegex(c.ControllerError, "source hash changed"):
            isolated.source_hashes(config)


if __name__ == "__main__":
    unittest.main(verbosity=2)
