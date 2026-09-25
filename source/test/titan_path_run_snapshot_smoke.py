#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from titan_path_run_snapshot import (  # noqa: E402
    SnapshotError,
    _normalize_old_direction,
    _parser,
    _write_old_direction,
    _build_titan_command,
    _nav_context,
    NAV_P2_OPTIONS,
    main,
)


class TitanPathRunSnapshotSmoke(unittest.TestCase):
    def _fixture(self, base: Path):
        route=base/"route.tpath"
        route.write_text("synthetic route\n")
        arg=base/"arg.in"
        arg.write_text("-pbcn\n-ff=CHARMM36\n")
        status=base/"previous.txt"
        status.write_text("SCHEMA=TITAN_PATH_STATUS_V1\nSTATUS=CONTINUE\nMODE=Recovery\n"
            "ROUTE_ID=route-A\nROUTE_FILE_SHA256="+hashlib.sha256(route.read_bytes()).hexdigest()+"\n"
            "SNAPSHOT_TIME_PS=20\nPATH_PROGRESS_A=5\nPATH_TOTAL_A=25\n"
            "NEXT_VECTOR_X_A=1\nNEXT_VECTOR_Y_A=0\nNEXT_VECTOR_Z_A=0\n")
        args=argparse.Namespace(profile="nav-p2",route_id="route-A",time_ps=30,
            previous_status=str(status),previous_status_sha256=hashlib.sha256(status.read_bytes()).hexdigest(),progress_window_A=1.)
        return args,arg,route,status

    def test_nav_profile_matches_replay_and_binds_continuity(self):
        with tempfile.TemporaryDirectory() as temporary:
            base=Path(temporary)
            args,arg,route,status=self._fixture(base)
            with patch.dict(os.environ,{"GMXLIB":str(base)}):
                context=_nav_context(args,arg,route,(1.,0.,0.))
            command=_build_titan_command(args,titan=base/"titan",arg_file=arg,
                snapshot=base/"snapshot.pdb",topology=base/"topol.top",index=base/"index.ndx",
                route=route,fflib=base,nav_context=context)
            self.assertIn("-common_topo",command)
            self.assertNotIn("-topo",command)
            start=command.index("-pathway-probe-lengths")
            self.assertEqual(command[start:start+len(NAV_P2_OPTIONS)],list(NAV_P2_OPTIONS))
            self.assertEqual(command[-4:],["-pathway-previous-status",str(status.resolve()),"-pathway-continuity-window","1"])
            self.assertEqual(context["previous_progress_A"],5.)

    def test_nav_rejects_wrong_hash_missing_state_and_unbound_direction(self):
        with tempfile.TemporaryDirectory() as temporary:
            base=Path(temporary)
            args,arg,route,status=self._fixture(base)
            with patch.dict(os.environ,{"GMXLIB":str(base)}):
                with self.assertRaisesRegex(SnapshotError,"old direction"):
                    _nav_context(args,arg,route,(0.,1.,0.))
                args.previous_status_sha256="0"*64
                with self.assertRaisesRegex(SnapshotError,"SHA256 mismatch"):
                    _nav_context(args,arg,route,(1.,0.,0.))
                args.previous_status=None
                with self.assertRaisesRegex(SnapshotError,"hash-bound"):
                    _nav_context(args,arg,route,(1.,0.,0.))

    def test_nav_rejects_future_terminal_duplicate_status_and_bad_environment(self):
        with tempfile.TemporaryDirectory() as temporary:
            base=Path(temporary)
            args,arg,route,status=self._fixture(base)
            original=status.read_text()
            for old,new in [("STATUS=CONTINUE","STATUS=STOP"),
                            ("SNAPSHOT_TIME_PS=20","SNAPSHOT_TIME_PS=40"),
                            ("PATH_PROGRESS_A=5","PATH_PROGRESS_A=nan"),
                            ("MODE=Recovery","MODE=Recovery\nMODE=Advance")]:
                status.write_text(original.replace(old,new))
                args.previous_status_sha256=hashlib.sha256(status.read_bytes()).hexdigest()
                with self.assertRaises(SnapshotError),patch.dict(os.environ,{"GMXLIB":str(base)}):
                    _nav_context(args,arg,route,(1.,0.,0.))
            status.write_text(original)
            args.previous_status_sha256=hashlib.sha256(status.read_bytes()).hexdigest()
            with self.assertRaisesRegex(SnapshotError,"GMXLIB"),patch.dict(os.environ,{"GMXLIB":""}):
                _nav_context(args,arg,route,(1.,0.,0.))
            arg.write_text("-pbctopo\n")
            with self.assertRaisesRegex(SnapshotError,"pbcn"),patch.dict(os.environ,{"GMXLIB":str(base)}):
                _nav_context(args,arg,route,(1.,0.,0.))

    def test_legacy_command_remains_explicitly_legacy(self):
        args=argparse.Namespace(profile="legacy",route_id="r")
        p=Path("dummy")
        command=_build_titan_command(args,titan=p,arg_file=p,snapshot=p,topology=p,
                                    index=p,route=p,fflib=p)
        self.assertIn("-topo",command)
        self.assertNotIn("-common_topo",command)
        self.assertNotIn("-pathway-navigation",command)

    def test_large_finite_old_vector_is_normalized_without_overflow(self):
        self.assertEqual(_normalize_old_direction([1e200,0,0]),(1.,0.,0.))
        with self.assertRaises(SnapshotError):
            _normalize_old_direction([1.7e308,1.7e308,0])

    def test_mocked_full_adapter_records_hashes_and_rejects_threshold_drift(self):
        for threshold in ("1.5", "1.4"):
            with self.subTest(threshold=threshold), tempfile.TemporaryDirectory() as temporary:
                base=Path(temporary)
                args,arg,route,status=self._fixture(base)
                for name in ("titan", "topol.top", "index.ndx"):
                    (base/name).write_text("SYNTHETIC fixture, not a native input\n")
                (base/"state.gro").write_text("SYNTHETIC\n2\nunused\nunused\n1 1 1\n")
                work=base/"new_segment"
                work.mkdir()
                snapshot=work/"snapshot.pdb"
                actual={"SCHEMA":"TITAN_PATH_STATUS_V1","STATUS":"CONTINUE",
                    "NAVIGATION_ENABLED":"1","GEOMETRY_ONLY":"0","FIXED_PROBE_ENABLED":"1",
                    "CLEARANCE_MODE":"swept-linear-center","CLEARANCE_PBC_PREPROCESSOR":"pbcn",
                    "PROJECTION_CONTINUITY":"BOUNDED_PREVIOUS_PROGRESS","ROUTE_ID":"route-A",
                    "ROUTE_FILE_SHA256":hashlib.sha256(route.read_bytes()).hexdigest(),
                    "MIN_CLEARANCE_THRESHOLD_A":threshold,"ENERGY_KEY_EPSILON":".05",
                    "ADVANCE_MAX_TURN_ANGLE_DEG":"20","RECOVERY_MAX_TURN_ANGLE_DEG":"60",
                    "NAVIGATION_TUBE_WEIGHT":"1","NAVIGATION_BACKTRACK_WEIGHT":"1",
                    "NAVIGATION_COST_EPSILON":".02","NAVIGATION_SEED_RADIAL_SAMPLES":"8",
                    "NAVIGATION_SEED_AMPLITUDE":".35","CONFIGURED_PROBE_LENGTHS_A":"1,.5,.25",
                    "SNAPSHOT_TIME_PS":"30","PREVIOUS_PATH_PROGRESS_A":"5"}
                def convert(*unused,**kw):
                    snapshot.write_text("TITLE fake\nATOM      1\nATOM      2\nEND\n")
                    return 0.001
                def fake_titan(command,**kw):
                    (work/"pull_direction_status.txt").write_text("".join(k+"="+v+"\n" for k,v in actual.items()))
                    (work/"DE_out_vector.dat").write_text("0.5 0 0\n")
                    return subprocess.CompletedProcess(command,0)
                argv=["--titan",str(base/"titan"),"--gro",str(base/"state.gro"),
                    "--snapshot",str(snapshot),"--time-ps","30","--old-direction","1","0","0",
                    "--arg-file",str(arg),"--topology",str(base/"topol.top"),"--index",str(base/"index.ndx"),
                    "--route",str(route),"--route-id","route-A","--fflib",str(base),"--profile","nav-p2",
                    "--previous-status",str(status),"--previous-status-sha256",args.previous_status_sha256,
                    "--progress-window-A","1"]
                prior=Path.cwd()
                try:
                    os.chdir(work)
                    with patch.dict(os.environ,{"GMXLIB":str(base)}), patch("titan_path_run_snapshot._run_checked",side_effect=convert), patch("titan_path_run_snapshot.subprocess.run",side_effect=fake_titan) as mocked:
                        code=main(argv)
                        self.assertEqual(mocked.call_count,1)
                finally:
                    os.chdir(prior)
                if threshold=="1.5":
                    self.assertEqual(code,0)
                    runtime=json.loads((work/"direction_provider_runtime.json").read_text())
                    self.assertEqual(runtime["input_sha256"]["snapshot"],hashlib.sha256(snapshot.read_bytes()).hexdigest())
                    self.assertEqual(runtime["output_sha256"]["status"],hashlib.sha256((work/"pull_direction_status.txt").read_bytes()).hexdigest())
                    self.assertIn("-common_topo",runtime["command"])
                else:
                    self.assertEqual(code,1)
                    self.assertFalse((work/"direction_provider_runtime.json").exists())
    def test_old_direction_is_normalized_and_atomically_replaces_stale_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "DE_out_vector.dat"
            path.write_text("stale direction\n", encoding="utf-8")
            direction = _normalize_old_direction([0.0, 3.0, 4.0])
            self.assertEqual(direction, (0.0, 0.6, 0.8))
            _write_old_direction(path, direction)
            self.assertEqual(
                [float(value) for value in path.read_text(encoding="utf-8").split()],
                [0.0, 0.6, 0.8],
            )
            self.assertEqual(list(path.parent.glob("DE_out_vector.dat.tmp.*")), [])

    def test_zero_old_direction_is_rejected(self) -> None:
        with self.assertRaisesRegex(SnapshotError, "must be non-zero"):
            _normalize_old_direction([0.0, 0.0, 0.0])

    def test_cli_requires_three_old_direction_components(self) -> None:
        required = [
            "--titan", "titan", "--gro", "state.gro", "--snapshot", "state.pdb",
            "--time-ps", "20", "--old-direction", "1", "0", "0",
            "--arg-file", "arg.in", "--topology", "topol.top", "--index", "index.ndx",
            "--route", "route.tpath", "--route-id", "route-A", "--fflib", "fflib",
        ]
        parsed = _parser().parse_args(required)
        self.assertEqual(parsed.old_direction, [1.0, 0.0, 0.0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
