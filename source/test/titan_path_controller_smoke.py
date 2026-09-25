#!/usr/bin/env python3

from __future__ import annotations

import json
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from titan_path_controller import (  # noqa: E402
    CONFIG_SCHEMA,
    ControllerError,
    _initial_state,
    _normalize_resume_controller_contract,
    expand_command,
    load_and_validate_config,
    parse_key_value_file,
    read_xvg_boundary,
    reference_displacement_steps,
    run_controller,
    status_vector,
    validate_continue_route_status,
    validate_route_probe_contract,
    validate_route_stop_diagnostics,
)


class TitanPathControllerSmoke(unittest.TestCase):
    def test_reference_displacement_contract(self) -> None:
        for displacement, steps, duration in (
            (1.0, 10000, 20.0),
            (0.5, 5000, 10.0),
            (0.25, 2500, 5.0),
        ):
            represented = reference_displacement_steps(
                displacement, rate_nm_per_ps=0.005, dt_ps=0.002
            )
            self.assertEqual(represented["steps"], steps)
            self.assertAlmostEqual(represented["duration_ps"], duration)
            self.assertAlmostEqual(represented["residual_A"], 0.0)
        with self.assertRaisesRegex(ControllerError, "integer MD step"):
            reference_displacement_steps(
                0.25005, rate_nm_per_ps=0.005, dt_ps=0.002
            )

        config = {
            "reference_displacement_contract": "match_titan_probe",
            "allowed_probe_lengths_A": [1.0, 0.5, 0.25],
            "rate_nm_per_ps": 0.005,
            "dt_ps": 0.002,
            "segment_steps": 10000,
        }
        status = {
            "FIXED_PROBE_ENABLED": "1",
            "CONFIGURED_PROBE_LENGTH_A": "1.0",
            "CONFIGURED_PROBE_LENGTHS_A": "1.0,0.5,0.25",
            "NEXT_REFERENCE_DISPLACEMENT_A": "0.5",
            "NEXT_VECTOR_X_A": "0.3",
            "NEXT_VECTOR_Y_A": "0.4",
            "NEXT_VECTOR_Z_A": "0",
            "APPLIED_PROGRESS_A": "1.0",
            "APPLIED_DELTA_PROGRESS_A": "0.5",
            "APPLIED_PATH_DISTANCE_A": "0.2",
        }
        checked = validate_route_probe_contract(status, config)
        self.assertEqual(checked["next_segment_steps"], 5000)
        self.assertAlmostEqual(checked["next_segment_duration_ps"], 10.0)
        for changes, message in (
            ({"NEXT_REFERENCE_DISPLACEMENT_A": "0"}, "positive"),
            ({"NEXT_REFERENCE_DISPLACEMENT_A": "0.25"}, "NEXT_VECTOR norm"),
            (
                {"CONFIGURED_PROBE_LENGTHS_A": "1.0,0.25"},
                "configured probe ladder",
            ),
            (
                {"CONFIGURED_PROBE_LENGTH_A": "0.5"},
                "first allowed_probe_lengths_A",
            ),
            (
                {
                    "NEXT_REFERENCE_DISPLACEMENT_A": "0.75",
                    "CONFIGURED_PROBE_LENGTH_A": "0.75",
                    "NEXT_VECTOR_X_A": "0.45",
                    "NEXT_VECTOR_Y_A": "0.6",
                },
                "allowed_probe_lengths",
            ),
        ):
            with self.assertRaisesRegex(ControllerError, message):
                validate_route_probe_contract({**status, **changes}, config)
        with self.assertRaisesRegex(ControllerError, "requires TITAN fixed probe"):
            validate_route_probe_contract(
                {**status, "FIXED_PROBE_ENABLED": "0"}, config
            )

    def test_xvg_boundary_uses_explicit_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "pullx.xvg"
            path.write_text(
                "# comment\n@ legend\n"
                "0.0 1.0 0.8 1 0 1 0 0 0 1 0 1\n"
                "20.0 1.2 0.9 1 0 1.2 0 0 0 1 0 1.2\n",
                encoding="utf-8",
            )
            boundary = read_xvg_boundary(
                path,
                time_column=0,
                coordinate_column=1,
                reference_column=2,
                coordinate_component_columns=[3, 4, 5],
                reference_com_columns=[6, 7, 8],
                ligand_com_columns=[9, 10, 11],
            )
            self.assertEqual(boundary["time_ps"], 20.0)
            self.assertAlmostEqual(boundary["coordinate_nm"], 1.2)
            self.assertAlmostEqual(boundary["reference_nm"], 0.9)
            self.assertEqual(boundary["coordinate_components_nm"], (1.0, 0.0, 1.2))
            self.assertEqual(boundary["reference_com_nm"], (0.0, 0.0, 0.0))
            self.assertEqual(boundary["ligand_com_nm"], (1.0, 0.0, 1.2))

    def test_status_duplicate_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "status.txt"
            path.write_text("STATUS=CONTINUE\nSTATUS=STOP\n", encoding="utf-8")
            with self.assertRaises(ControllerError):
                parse_key_value_file(path)

    def test_status_vector_and_command_expansion(self) -> None:
        status = {
            "NEXT_VECTOR_X_A": "1",
            "NEXT_VECTOR_Y_A": "2",
            "NEXT_VECTOR_Z_A": "3",
        }
        self.assertEqual(status_vector(status, "NEXT_VECTOR"), (1.0, 2.0, 3.0))
        self.assertEqual(
            expand_command(["titan", "-f", "{gro}"], {"gro": "state.gro"}),
            ["titan", "-f", "state.gro"],
        )

    def test_route_stop_diagnostics_contract(self) -> None:
        scan = {
            "STATUS": "STOP",
            "REASON": "NO_VALID_FORWARD_DIRECTION",
            "REASON_DETAIL": "SCAN_SHELL_EXHAUSTED",
            "CANDIDATE_PROPOSALS": "12",
            "CANDIDATE_EVALUATIONS": "0",
            "SCAN_SHELL_EXCLUSIONS": "12",
            "CLEARANCE_REJECTIONS": "0",
            "ENERGY_EVALUATIONS": "0",
        }
        self.assertEqual(
            validate_route_stop_diagnostics(scan), "SCAN_SHELL_EXHAUSTED"
        )
        clearance = {
            **scan,
            "REASON_DETAIL": "CLEARANCE_EXHAUSTED",
            "CANDIDATE_PROPOSALS": "10",
            "CANDIDATE_EVALUATIONS": "6",
            "SCAN_SHELL_EXCLUSIONS": "4",
            "CLEARANCE_REJECTIONS": "6",
        }
        self.assertEqual(
            validate_route_stop_diagnostics(clearance), "CLEARANCE_EXHAUSTED"
        )
        self.assertIsNone(validate_route_stop_diagnostics({"STATUS": "STOP"}))
        invalid = {**clearance, "CLEARANCE_REJECTIONS": "5"}
        with self.assertRaisesRegex(ControllerError, "exhaustion contract"):
            validate_route_stop_diagnostics(invalid)
        recovery = {
            "STATUS": "STOP",
            "REASON": "NO_RECOVERY_DESCENT_DIRECTION",
            "REASON_DETAIL": "RECOVERY_DESCENT_EXHAUSTED",
            "MODE": "Recovery",
            "FIXED_PROBE_ENABLED": "1",
            "RECOVERY_DIRECT_PROBE_USED": "1",
        }
        self.assertEqual(
            validate_route_stop_diagnostics(recovery),
            "RECOVERY_DESCENT_EXHAUSTED",
        )
        ordinary_recovery = {**recovery, "RECOVERY_DIRECT_PROBE_USED": "0"}
        self.assertEqual(
            validate_route_stop_diagnostics(ordinary_recovery),
            "RECOVERY_DESCENT_EXHAUSTED",
        )
        with self.assertRaisesRegex(ControllerError, "fixed-probe Recovery STOP"):
            validate_route_stop_diagnostics(
                {**recovery, "FIXED_PROBE_ENABLED": "0"}
            )
        progress = {
            "STATUS": "STOP", "MODE": "Advance",
            "REASON": "NO_FORWARD_PROGRESS_CANDIDATE",
            "REASON_DETAIL": "FORWARD_PROGRESS_EXHAUSTED",
            "VALID_ARCHIVE_SIZE": "236", "FORWARD_ARCHIVE_SIZE": "0",
            "FORWARD_PROGRESS_REJECTIONS": "236", "EPSILON_ELIGIBLE_SIZE": "0",
            "NAVIGATION_ELIGIBLE_SIZE": "0",
        }
        self.assertEqual(validate_route_stop_diagnostics(progress), "FORWARD_PROGRESS_EXHAUSTED")
        for changed in ({"MODE": "Recovery"}, {"STATUS": "CONTINUE"}, {"VALID_ARCHIVE_SIZE": "0"}, {"FORWARD_ARCHIVE_SIZE": "1"}, {"FORWARD_PROGRESS_REJECTIONS": "235"}, {"EPSILON_ELIGIBLE_SIZE": "1"}, {"NAVIGATION_ELIGIBLE_SIZE": "1"}):
            with self.assertRaises(ControllerError):
                validate_route_stop_diagnostics({**progress, **changed})

    def test_advance_status_requires_strict_forward_selection(self) -> None:
        valid = {
            "MODE": "Advance",
            "CLEARANCE_MODE": "endpoint-center",
            "PATH_PROGRESS_A": "0.0",
            "TARGET_PROGRESS_A": "0.25",
            "VALID_ARCHIVE_SIZE": "25",
            "FORWARD_ARCHIVE_SIZE": "17",
            "EPSILON_ELIGIBLE_SIZE": "1",
        }
        checked = validate_continue_route_status(valid)
        self.assertEqual(checked["reported_forward_archive_size"], 17)
        for changes in (
            {"TARGET_PROGRESS_A": "0.0"},
            {"FORWARD_ARCHIVE_SIZE": "0"},
            {"EPSILON_ELIGIBLE_SIZE": "18"},
        ):
            invalid = {**valid, **changes}
            with self.assertRaises(ControllerError):
                validate_continue_route_status(invalid)
        navigation = {
            **valid,
            "NAVIGATION_ENABLED": "1",
            "NAVIGATION_COST_KIND": "ALGORITHMIC_DIMENSIONLESS",
            "NAVIGATION_COST_EPSILON": "0.02",
            "TARGET_NAVIGATION_TUBE_COST": "0.01",
            "TARGET_NAVIGATION_BACKTRACK_COST": "0.0",
            "TARGET_NAVIGATION_COST": "0.01",
            "NAVIGATION_ELIGIBLE_SIZE": "1",
            "ENERGY_ELIGIBLE_MIN_NAVIGATION_COST": "0.0",
        }
        checked_navigation = validate_continue_route_status(navigation)
        self.assertTrue(checked_navigation["navigation_enabled"])
        with self.assertRaisesRegex(ControllerError, "decomposition"):
            validate_continue_route_status(
                {**navigation, "TARGET_NAVIGATION_COST": "0.02"}
            )
        with self.assertRaisesRegex(ControllerError, "energy pool"):
            validate_continue_route_status(
                {**navigation, "NAVIGATION_ELIGIBLE_SIZE": "20"}
            )
        swept = {
            **navigation,
            "CLEARANCE_MODE": "swept-linear-center",
            "CLEARANCE_SCOPE": "CONTINUOUS_LINEAR_TRANSLATION",
            "CLEARANCE_COORDINATE_FRAME": "UPSTREAM_PBC_RECONSTRUCTED_CARTESIAN",
            "CLEARANCE_PERIODIC_DISTANCE_POLICY": "DIRECT_IN_RECONSTRUCTED_FRAME",
            "CLEARANCE_PBC_PREPROCESSOR": "pbcn",
            "SWEPT_CLEARANCE_PROBE_EVALUATIONS": "12",
            "SWEPT_CLEARANCE_PAIR_CHECKS": "1200",
            "SWEPT_CLEARANCE_REJECTIONS": "2",
            "SWEPT_INTERIOR_ONLY_REJECTIONS": "1",
            "TARGET_SWEPT_MIN_CLEARANCE_A": "1.6",
            "TARGET_ENDPOINT_MIN_CLEARANCE_A": "2.0",
            "TARGET_MIN_CLEARANCE_A": "1.6",
            "MIN_CLEARANCE_THRESHOLD_A": "1.5",
            "TARGET_SWEPT_MIN_FRACTION": "0.4",
            "TARGET_SWEPT_MIN_COMPLETE": "1",
        }
        checked_swept = validate_continue_route_status(swept)
        self.assertEqual(checked_swept["clearance_mode"], "swept-linear-center")
        with self.assertRaisesRegex(ControllerError, "witness"):
            validate_continue_route_status(
                {**swept, "TARGET_ENDPOINT_MIN_CLEARANCE_A": "1.4"}
            )
        with self.assertRaisesRegex(ControllerError, "controller config"):
            validate_continue_route_status(
                swept, expected_clearance_mode="endpoint-center"
            )
        with self.assertRaisesRegex(ControllerError, "CLEARANCE_MODE"):
            validate_continue_route_status(
                {key: value for key, value in valid.items() if key != "CLEARANCE_MODE"}
            )
        recovery = {
            "MODE": "Recovery",
            "CLEARANCE_MODE": "endpoint-center",
            "PATH_PROGRESS_A": "5.0",
            "TARGET_PROGRESS_A": "5.5",
            "REQUESTED_TARGET_PROGRESS_A": "5.5",
            "APPLIED_PROGRESS_A": "5.1",
            "APPLIED_DELTA_PROGRESS_A": "0.1",
            "PATH_DISTANCE_A": "4.5",
            "APPLIED_PATH_DISTANCE_A": "3.7",
            "FIXED_PROBE_ENABLED": "1",
            "RECOVERY_DIRECT_PROBE_USED": "1",
            "VALID_ARCHIVE_SIZE": "1",
            "FORWARD_ARCHIVE_SIZE": "1",
            "EPSILON_ELIGIBLE_SIZE": "1",
        }
        checked_recovery = validate_continue_route_status(recovery)
        self.assertEqual(checked_recovery["mode"], "Recovery")
        with self.assertRaisesRegex(ControllerError, "does not approach"):
            validate_continue_route_status(
                {**recovery, "APPLIED_PATH_DISTANCE_A": "4.5"}
            )
        with self.assertRaisesRegex(ControllerError, "inconsistent"):
            validate_continue_route_status(
                {**recovery, "REQUESTED_TARGET_PROGRESS_A": "5.4"}
            )

    def test_pre_nav_resume_contract_normalization(self) -> None:
        old = {
            "max_turn_angle_deg": 20.0,
            "direction_kind": "explicit_route",
        }
        normalized = _normalize_resume_controller_contract(old)
        self.assertEqual(normalized["recovery_max_turn_angle_deg"], 20.0)
        self.assertEqual(
            normalized["reference_displacement_contract"],
            "legacy_direction_only",
        )
        self.assertEqual(normalized["allowed_probe_lengths_A"], [])
        self.assertEqual(normalized["expected_clearance_mode"], "endpoint-center")
        self.assertNotIn("reference_displacement_contract", old)

    def test_config_preflight_resolves_and_hashes_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in (
                "template.mdp",
                "topol.top",
                "route.tpath",
                "initial.tpr",
                "index.ndx",
            ):
                (root / name).write_text(
                    "integrator = md\n" if name.endswith(".mdp") else name,
                    encoding="utf-8",
                )
            config = {
                "schema": CONFIG_SCHEMA,
                "work_dir": "run",
                "gmx_command": ["gmx"],
                "titan_command": ["titan", "-f", "{gro}"],
                "template_mdp": "template.mdp",
                "topology": "topol.top",
                "route_file": "route.tpath",
                "expected_route_id": "route-A",
                "expected_reference_sha256": "0" * 64,
                "index": "index.ndx",
                "initial_tpr": "initial.tpr",
                "initial_checkpoint": None,
                "initial_direction": [0, 0, 2],
                "dt_ps": 0.002,
                "segment_steps": 100,
                "max_segments": 2,
                "rate_nm_per_ps": 0.001,
                "spring_k_kj_mol_nm2": 1000,
                "pullx_columns": {
                    "time": 0,
                    "coordinate": 1,
                    "reference": 2,
                    "coordinate_components": [3, 4, 5],
                    "reference_com": [6, 7, 8],
                    "ligand_com": [9, 10, 11],
                },
            }
            path = root / "controller.json"
            path.write_text(json.dumps(config), encoding="utf-8")
            loaded = load_and_validate_config(path)
            self.assertEqual(loaded["initial_direction"], (0.0, 0.0, 1.0))
            self.assertEqual(
                loaded["recovery_max_turn_angle_deg"],
                loaded["max_turn_angle_deg"],
            )
            self.assertEqual(loaded["expected_clearance_mode"], "endpoint-center")
            self.assertEqual(loaded["template_mdp"], (root / "template.mdp").resolve())

            nonroute_matched = {
                **config,
                "direction_kind": "fixed",
                "expected_provider_id": "fixed-v1",
                "reference_displacement_contract": "match_titan_probe",
                "allowed_probe_lengths_A": [0.002],
                "recovery_max_turn_angle_deg": 60,
            }
            nonroute_path = root / "controller_nonroute_matched.json"
            nonroute_path.write_text(json.dumps(nonroute_matched), encoding="utf-8")
            with self.assertRaisesRegex(ControllerError, "explicit_route"):
                load_and_validate_config(nonroute_path)

            duplicate_template = root / "template_duplicate.mdp"
            duplicate_template.write_text(
                "integrator = md\ndt = 0.002\ndt = 0.004\n",
                encoding="utf-8",
            )
            duplicate_path = root / "controller_duplicate_mdp.json"
            duplicate_path.write_text(
                json.dumps({**config, "template_mdp": duplicate_template.name}),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ControllerError, "duplicate MDP key 'dt'"):
                load_and_validate_config(duplicate_path)

    def test_fake_two_segment_controller_cycle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            route = root / "route.tpath"
            route.write_text("fake-route\n", encoding="utf-8")
            for name in ("topol.top", "initial.tpr", "index.ndx"):
                (root / name).write_text(name, encoding="utf-8")
            (root / "template.mdp").write_text(
                "integrator = md\ndt = 0.004\npull = yes\n"
                "pull-coord1-rate = 0.9\n",
                encoding="utf-8",
            )
            fake_gmx = root / "fake_gmx.py"
            fake_gmx.write_text(
                """from pathlib import Path
import sys

def value(flag):
    return Path(sys.argv[sys.argv.index(flag) + 1])

mode = sys.argv[1]
if mode == 'mdrun':
    prefix = value('-deffnm')
    gro = prefix.with_suffix('.gro')
    cpt = value('-cpo')
    pullx = value('-px')
    gro.parent.mkdir(parents=True, exist_ok=True)
    gro.write_text('fake gro\\n', encoding='utf-8')
    cpt.write_text('fake checkpoint\\n', encoding='utf-8')
    previous = 0.0
    if pullx.exists():
        rows = [line for line in pullx.read_text(encoding='utf-8').splitlines()
                if line and line[0] not in '#@']
        if rows:
            previous = float(rows[-1].split()[0])
    with pullx.open('a', encoding='utf-8') as handle:
        handle.write(
            f'{previous + 0.2:.3f} 1.0 0.9 1 0 1.2 0 0 0 1 0 1.2\\n'
        )
elif mode == 'grompp':
    output = value('-o')
    output.write_text('fake tpr\\n', encoding='utf-8')
else:
    raise SystemExit(9)
""",
                encoding="utf-8",
            )
            fake_titan = root / "fake_titan.py"
            fake_titan.write_text(
                """from pathlib import Path
import hashlib
import sys

route = Path(sys.argv[1])
time_ps = float(sys.argv[2])
behavior = sys.argv[3] if len(sys.argv) > 3 else 'complete'
used_old_vector = '0' if behavior == 'bad-old-vector' else '1'
probe_length = 0.001 if behavior in ('probe-half', 'recovery-turn') else 1.0
segment = int(Path.cwd().name.rsplit('_', 1)[1])
if segment >= 1 and behavior in ('stop', 'stop-clearance', 'stop-scan'):
    status = 'STOP'
    reason = 'NO_VALID_FORWARD_DIRECTION'
elif segment >= 1:
    status = 'COMPLETE'
    reason = 'ROUTE_TERMINAL_REACHED'
else:
    status = 'CONTINUE'
    reason = 'OK'
lines = [
    'SCHEMA=TITAN_PATH_STATUS_V1', f'STATUS={status}', f'REASON={reason}',
    'ROUTE_KIND=EXPLICIT', 'ROUTE_ID=route-A', 'REFERENCE_SHA256=' + '0' * 64,
    'ROUTE_FILE_SHA256=' + hashlib.sha256(route.read_bytes()).hexdigest(),
    'CLEARANCE_MODE=endpoint-center',
    'MAX_TURN_ANGLE_DEG=' + ('60' if behavior == 'recovery-turn' and segment == 0 else '20'),
    f'USED_OLD_VECTOR={used_old_vector}',
    f'SNAPSHOT_TIME_PS={time_ps}',
]
if status == 'STOP' and behavior == 'stop-clearance':
    lines += [
        'REASON_DETAIL=CLEARANCE_EXHAUSTED',
        'CANDIDATE_PROPOSALS=10', 'CANDIDATE_EVALUATIONS=6',
        'SCAN_SHELL_EXCLUSIONS=4',
        'CLEARANCE_REJECTIONS=6', 'ENERGY_EVALUATIONS=0',
    ]
elif status == 'STOP' and behavior == 'stop-scan':
    lines += [
        'REASON_DETAIL=SCAN_SHELL_EXHAUSTED',
        'CANDIDATE_PROPOSALS=12', 'CANDIDATE_EVALUATIONS=0',
        'SCAN_SHELL_EXCLUSIONS=12',
        'CLEARANCE_REJECTIONS=0', 'ENERGY_EVALUATIONS=0',
    ]
if status == 'CONTINUE':
    mode = 'Recovery' if behavior == 'recovery-turn' else 'Advance'
    angle = 40.0 if behavior == 'recovery-turn' else 0.0
    x = probe_length * (0.766044443118978 if behavior == 'recovery-turn' else 1.0)
    y = probe_length * (0.6427876096865394 if behavior == 'recovery-turn' else 0.0)
    lines += [
        f'MODE={mode}', 'PATH_PROGRESS_A=0', 'TARGET_PROGRESS_A=0.25',
        'VALID_ARCHIVE_SIZE=2', 'FORWARD_ARCHIVE_SIZE=2',
        'EPSILON_ELIGIBLE_SIZE=1',
        f'NEXT_VECTOR_X_A={x}', f'NEXT_VECTOR_Y_A={y}',
        'NEXT_VECTOR_Z_A=0',
        f'NEXT_VECTOR_ANGLE_DEG={angle}',
        'ROUTE_TARGET_VECTOR_X_A=1', 'ROUTE_TARGET_VECTOR_Y_A=0',
        'ROUTE_TARGET_VECTOR_Z_A=0', 'ROUTE_TARGET_ANGLE_DEG=0',
        'TURN_CAPPED=0', 'TARGET_CLEARANCE_VECTOR=NEXT_VECTOR',
        'REFERENCE_COM_X_A=0', 'REFERENCE_COM_Y_A=0', 'REFERENCE_COM_Z_A=0',
        'LIGAND_COM_X_A=20', 'LIGAND_COM_Y_A=0', 'LIGAND_COM_Z_A=0',
    ]
    if behavior in ('probe-half', 'recovery-turn'):
        lines += [
            'FIXED_PROBE_ENABLED=1', 'CONFIGURED_PROBE_LENGTH_A=0.002',
            'CONFIGURED_PROBE_LENGTHS_A=0.002,0.001,0.0005',
            'NEXT_REFERENCE_DISPLACEMENT_A=0.001',
            'APPLIED_PROGRESS_A=0.25', 'APPLIED_DELTA_PROGRESS_A=0.25',
            'APPLIED_PATH_DISTANCE_A=0.1',
        ]
    if behavior == 'recovery-turn':
        lines += ['PATH_DISTANCE_A=0.2', 'RECOVERY_DIRECT_PROBE_USED=0']
Path('pull_direction_status.txt').write_text('\\n'.join(lines) + '\\n', encoding='utf-8')
""",
                encoding="utf-8",
            )
            config = {
                "schema": CONFIG_SCHEMA,
                "work_dir": "run",
                "gmx_command": [sys.executable, str(fake_gmx)],
                "titan_command": [
                    sys.executable,
                    str(fake_titan),
                    str(route),
                    "{time_ps}",
                ],
                "template_mdp": "template.mdp",
                "topology": "topol.top",
                "route_file": "route.tpath",
                "expected_route_id": "route-A",
                "expected_reference_sha256": "0" * 64,
                "index": "index.ndx",
                "initial_tpr": "initial.tpr",
                "initial_checkpoint": None,
                "initial_direction": [1, 0, 0],
                "dt_ps": 0.002,
                "segment_steps": 100,
                "max_segments": 3,
                "rate_nm_per_ps": 0.001,
                "spring_k_kj_mol_nm2": 1000,
                "max_turn_angle_deg": 20,
                "pullx_columns": {
                    "time": 0,
                    "coordinate": 1,
                    "reference": 2,
                    "coordinate_components": [3, 4, 5],
                    "reference_com": [6, 7, 8],
                    "ligand_com": [9, 10, 11],
                },
            }
            config_path = root / "controller.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            loaded = load_and_validate_config(config_path)
            state = run_controller(loaded, resume=False)
            self.assertEqual(state["status"], "COMPLETE")
            state_bytes = (loaded["work_dir"] / "controller_state.json").read_bytes()
            with patch("titan_path_controller._run_command", side_effect=AssertionError("terminal resume attempted an external command")):
                self.assertEqual(run_controller(loaded, resume=True), json.loads(json.dumps(state)))
            self.assertEqual((loaded["work_dir"] / "controller_state.json").read_bytes(), state_bytes)
            self.assertTrue(Path(state["checkpoint"]).is_file())
            terminal_record = json.loads(
                Path(state["last_terminal_record"]).read_text(encoding="utf-8")
            )
            self.assertEqual(
                terminal_record["schema"], "titan-path-controller-terminal-v1"
            )
            self.assertEqual(terminal_record["outcome"], "COMPLETE")
            record_path = root / "run" / "segments" / "segment_0000" / "segment_record.json"
            record = json.loads(record_path.read_text(encoding="utf-8"))
            self.assertAlmostEqual(record["handoff"]["extension_residual_nm"], 0.0)
            self.assertAlmostEqual(record["handoff"]["new_coordinate_nm"], 1.0)
            self.assertEqual(
                record["handoff_coordinate_source"],
                "GROMACS_PULLX_COORDINATE_COMPONENTS",
            )
            self.assertEqual(
                record["timing_wall_seconds"]["schema"],
                "titan-path-segment-timing-v1",
            )
            self.assertEqual(record["remaining_direction_increments"], 0)
            next_mdp = record_path.parent / "next.mdp"
            self.assertIn("pull-coord1-start = no", next_mdp.read_text(encoding="utf-8"))
            self.assertIn("dt = 0.002", next_mdp.read_text(encoding="utf-8"))
            self.assertIn(
                "pull-coord1-rate = 0.001", next_mdp.read_text(encoding="utf-8")
            )
            second_log = root / "run" / "segments" / "segment_0001" / "mdrun.log"
            self.assertIn("-cpi", second_log.read_text(encoding="utf-8"))

            probe_config = {
                **config,
                "work_dir": "run_probe_half",
                "reference_displacement_contract": "match_titan_probe",
                "allowed_probe_lengths_A": [0.002, 0.001, 0.0005],
                "recovery_max_turn_angle_deg": 60,
                "titan_command": [
                    sys.executable,
                    str(fake_titan),
                    str(route),
                    "{time_ps}",
                    "probe-half",
                ],
            }
            probe_path = root / "controller_probe_half.json"
            probe_path.write_text(json.dumps(probe_config), encoding="utf-8")
            # Old NAV fixtures lacked proof for their first TPR. Strict NAV
            # now rejects them; the native-controller smoke covers the fully
            # bound initial transfer and the shorter first segment.
            with self.assertRaisesRegex(ControllerError, "route_continuity_window_A"):
                load_and_validate_config(probe_path)

            recovery_turn_config = {
                **probe_config,
                "work_dir": "run_recovery_turn",
                "titan_command": [
                    sys.executable,
                    str(fake_titan),
                    str(route),
                    "{time_ps}",
                    "recovery-turn",
                ],
            }
            recovery_turn_path = root / "controller_recovery_turn.json"
            recovery_turn_path.write_text(
                json.dumps(recovery_turn_config), encoding="utf-8"
            )
            with self.assertRaisesRegex(ControllerError, "route_continuity_window_A"):
                load_and_validate_config(recovery_turn_path)

            legacy_config = {
                **config,
                "work_dir": "run_legacy",
                "spring_handoff_mode": "legacy_reset",
            }
            legacy_path = root / "controller_legacy.json"
            legacy_path.write_text(json.dumps(legacy_config), encoding="utf-8")
            legacy_state = run_controller(
                load_and_validate_config(legacy_path), resume=False
            )
            self.assertEqual(legacy_state["status"], "COMPLETE")
            legacy_record_path = (
                root
                / "run_legacy"
                / "segments"
                / "segment_0000"
                / "segment_record.json"
            )
            legacy_record = json.loads(
                legacy_record_path.read_text(encoding="utf-8")
            )
            self.assertEqual(legacy_record["spring_handoff_mode"], "legacy_reset")
            self.assertEqual(
                legacy_record["handoff"]["schema_version"],
                "titan-path-legacy-spring-reset-v1",
            )
            self.assertIsNotNone(legacy_record["force_continuous_counterfactual"])
            legacy_mdp = (legacy_record_path.parent / "next.mdp").read_text(
                encoding="utf-8"
            )
            self.assertIn("pull-coord1-start = yes", legacy_mdp)
            self.assertIn("pull-coord1-init = -0.00020000000000000001", legacy_mdp)

            stop_config = {
                **config,
                "work_dir": "run_stop",
                "titan_command": [
                    sys.executable,
                    str(fake_titan),
                    str(route),
                    "{time_ps}",
                    "stop",
                ],
            }
            stop_path = root / "controller_stop.json"
            stop_path.write_text(json.dumps(stop_config), encoding="utf-8")
            stop_state = run_controller(
                load_and_validate_config(stop_path), resume=False
            )
            self.assertEqual(stop_state["status"], "STOPPED")
            with patch("titan_path_controller._run_command", side_effect=AssertionError("STOPPED resume attempted an external command")):
                self.assertEqual(run_controller(load_and_validate_config(stop_path), resume=True), json.loads(json.dumps(stop_state)))
            self.assertEqual(
                stop_state["stop_reason"], "NO_VALID_FORWARD_DIRECTION"
            )
            stop_terminal_path = Path(stop_state["last_terminal_record"])
            stop_terminal = json.loads(
                stop_terminal_path.read_text(encoding="utf-8")
            )
            self.assertEqual(stop_terminal["outcome"], "STOP")
            self.assertEqual(
                stop_terminal["reason"], "NO_VALID_FORWARD_DIRECTION"
            )
            self.assertEqual(stop_terminal["segment"], 1)
            self.assertTrue(stop_terminal_path.with_name("terminal.complete").is_file())
            self.assertNotIn("stop_reason_detail", stop_state)
            self.assertNotIn("reason_detail", stop_terminal)

            for behavior, detail in (
                ("stop-clearance", "CLEARANCE_EXHAUSTED"),
                ("stop-scan", "SCAN_SHELL_EXHAUSTED"),
            ):
                detail_config = {
                    **config,
                    "work_dir": f"run_{behavior}",
                    "titan_command": [
                        sys.executable,
                        str(fake_titan),
                        str(route),
                        "{time_ps}",
                        behavior,
                    ],
                }
                detail_path = root / f"controller_{behavior}.json"
                detail_path.write_text(json.dumps(detail_config), encoding="utf-8")
                detail_state = run_controller(
                    load_and_validate_config(detail_path), resume=False
                )
                self.assertEqual(detail_state["status"], "STOPPED")
                self.assertEqual(detail_state["stop_reason_detail"], detail)
                detail_terminal = json.loads(
                    Path(detail_state["last_terminal_record"]).read_text(
                        encoding="utf-8"
                    )
                )
                self.assertEqual(detail_terminal["reason_detail"], detail)
                self.assertEqual(
                    detail_terminal["titan_status"]["REASON_DETAIL"], detail
                )

            bad_old_config = {
                **config,
                "work_dir": "run_bad_old_vector",
                "titan_command": [
                    sys.executable,
                    str(fake_titan),
                    str(route),
                    "{time_ps}",
                    "bad-old-vector",
                ],
            }
            bad_old_path = root / "controller_bad_old_vector.json"
            bad_old_path.write_text(json.dumps(bad_old_config), encoding="utf-8")
            with self.assertRaisesRegex(
                ControllerError, "did not evaluate against the active old direction"
            ):
                run_controller(
                    load_and_validate_config(bad_old_path), resume=False
                )

            resume_config = {**config, "work_dir": "run_interrupted_segment0"}
            resume_path = root / "controller_interrupted_segment0.json"
            resume_path.write_text(json.dumps(resume_config), encoding="utf-8")
            resume_loaded = load_and_validate_config(resume_path)
            resume_work = Path(resume_loaded["work_dir"])
            resume_work.mkdir(parents=True, exist_ok=True)
            (resume_work / "controller_state.json").write_text(
                json.dumps(_initial_state(resume_loaded)), encoding="utf-8"
            )
            orphan_checkpoint = resume_work / "run.cpt"
            orphan_checkpoint.write_text(
                "interrupted segment-0 checkpoint\n", encoding="utf-8"
            )
            resumed_state = run_controller(resume_loaded, resume=True)
            self.assertEqual(resumed_state["status"], "COMPLETE")
            adoption = resumed_state["checkpoint_adoptions"][0]
            self.assertEqual(
                adoption["reason"],
                "INTERRUPTED_SEGMENT0_CHECKPOINT_RECOVERY",
            )
            resume_log = (
                resume_work / "segments" / "segment_0000" / "mdrun.log"
            ).read_text(encoding="utf-8")
            self.assertIn("-cpi", resume_log)


if __name__ == "__main__":
    unittest.main(verbosity=2)
