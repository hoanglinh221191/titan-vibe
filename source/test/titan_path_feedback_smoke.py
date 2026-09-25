"""Synthetic boundary-consumer tests; no native TITAN or MD execution."""
from pathlib import Path
import copy
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"))
from titan_path_observation import Polyline, Tracker, digest
from titan_path_feedback import (FeedbackError, LoopPolicy, boundary_observation,
    consume_boundary, new_supervisor, paired_contract, seal, unseal)

IDENTITY = {"route_id": "synthetic", **{key: "a"*64 for key in
            ("route_sha256", "reference_sha256", "policy_sha256", "trajectory_sha256", "seed_sha256")}}
SELECTION = {"mode": "Advance", "source_kind": "NATIVE_TITAN_SELECTION", "sha256": "d"*64}


def qualification():
    # A test fixture, not a qualification of any real route/frame.
    return seal({"schema": "titan-path-frame-qualification-v1", "qualified": True,
                 "selection_basis": "PRE_SMD_REFERENCE_AND_EQUILIBRATION",
                 **{key: IDENTITY[key] for key in ("route_sha256", "reference_sha256", "policy_sha256")},
                 "input_sha256": {"SYNTHETIC": "e"*64},
                 "frame_definition_sha256": "f"*64, "validation_report_sha256": "f"*64})


def samples(start=0., progress=1., delta=.5, invalid_index=None):
    tracker = Tracker(IDENTITY, progress)
    route = Polyline([[0, 0, 0], [20, 0, 0]])
    return [tracker.observe(route, [progress+i*delta/25, 0, 0], start+i*.2,
            digest([start, i]), quality_reasons=["SYNTHETIC_INVALID"] if i == invalid_index else [])
            for i in range(26)]


def sensitivity(rows):
    return [{"observation_sha256": row["observation_sha256"], "time_ps": row["time_ps"],
             "reference_sensitive": False} for row in rows]


def boundary(rows=None, **kwargs):
    rows = samples() if rows is None else rows
    args = dict(identity=IDENTITY, start_ps=rows[0]["time_ps"], end_ps=rows[-1]["time_ps"],
                total_length_A=20., last_selection=SELECTION,
                endpoint_snapshot_sha256="b"*64, qualification=qualification())
    flags = kwargs.pop("sensitivity", sensitivity(rows))
    args.update(kwargs)
    return boundary_observation(rows, flags, **args)


def supervisor(arm="adaptive", **kwargs):
    args = dict(identity={k: v for k, v in IDENTITY.items() if k not in ("trajectory_sha256", "seed_sha256")},
                checkpoint_sha256="1"*64, execution_contract_sha256="3"*64,
                start_time_ps=0., initial_progress_A=1., arm=arm)
    args.update(kwargs)
    return new_supervisor(**args)


def consume(state, result, **kwargs):
    args = dict(expected_state_sha256=state["record_sha256"],
                seed_checkpoint_sha256=state["checkpoint_sha256"],
                endpoint_checkpoint_sha256=digest([result, "synthetic_checkpoint"]))
    args.update(kwargs)
    return consume_boundary(state, result, **args)


class FeedbackSmoke(unittest.TestCase):
    def test_valid_boundary_is_not_selection_or_md(self):
        b = boundary()
        self.assertEqual(b["state"], "SELECTOR_INPUT_READY")
        self.assertEqual(b["observed_progress_A"], 1.5)
        self.assertFalse(b["selection_decision_emitted"])
        self.assertFalse(b["md_authorized"])
        self.assertNotIn("STATUS", b)

    def test_probe_is_not_observed_movement(self):
        result = boundary(samples(delta=2.))
        self.assertAlmostEqual(result["observed_progress_A"], 3.)
        self.assertEqual(LoopPolicy().probe_A, .25)

    def test_backward_motion_stays_observed(self):
        result = boundary(samples(progress=3., delta=-1.))
        self.assertEqual(result["observed_progress_A"], 2.)

    def test_missing_frame_qualification_holds(self):
        self.assertIn("PHYSICAL_FRAME_NOT_QUALIFIED", boundary(qualification=None)["reasons"])

    def test_mismatched_qualification_rejected(self):
        q = unseal(qualification()); q["policy_sha256"] = "0"*64
        with self.assertRaisesRegex(FeedbackError, "IDENTITY"):
            boundary(qualification=seal(q))

    def test_outcome_selected_frame_rejected(self):
        q = unseal(qualification()); q["selection_basis"] = "BEST_NAV5_PROGRESS"
        with self.assertRaisesRegex(FeedbackError, "INDEPENDENT"):
            boundary(qualification=seal(q))

    def test_bare_qualified_boolean_rejected(self):
        with self.assertRaises((FeedbackError, TypeError)):
            boundary(qualification={"qualified": True})

    def test_midsegment_sensitivity_holds_even_valid_endpoint(self):
        rows = samples(); flags = sensitivity(rows); flags[3]["reference_sensitive"] = True
        self.assertEqual(boundary(rows, sensitivity=flags)["reasons"], ["REFERENCE_FRAME_SENSITIVE"])

    def test_broken_chain_never_promoted(self):
        result = boundary(samples(invalid_index=3))
        self.assertIn("OBSERVATION_CHAIN_UNRESOLVED", result["reasons"])
        self.assertIsNone(result["observed_progress_A"])

    def test_endpoint_only_sensitivity_rejected(self):
        rows = samples()
        with self.assertRaisesRegex(FeedbackError, "COVERAGE"):
            boundary(rows, sensitivity=sensitivity(rows)[-1:])

    def test_sensitivity_hash_mismatch(self):
        rows = samples(); flags = sensitivity(rows); flags[2]["observation_sha256"] = "b"*64
        with self.assertRaisesRegex(FeedbackError, "SENSITIVITY_IDENTITY"):
            boundary(rows, sensitivity=flags)

    def test_nonstrict_boolean_rejected(self):
        rows = samples(); flags = sensitivity(rows); flags[2]["reference_sensitive"] = "false"
        with self.assertRaisesRegex(FeedbackError, "SENSITIVITY_IDENTITY"):
            boundary(rows, sensitivity=flags)

    def test_truncated_span_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "TIME_MISMATCH"):
            boundary(start_ps=.2)

    def test_unmatched_probe_cadence_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "PROBE_CADENCE"):
            LoopPolicy(probe_A=1.)

    def test_integer_steps_required(self):
        with self.assertRaisesRegex(FeedbackError, "NONINTEGER"):
            LoopPolicy(dt_ps=.003)

    def test_model_selection_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "MODEL_SELECTION"):
            boundary(last_selection={**SELECTION, "source_kind": "RIGID_MODEL_ONLY"})

    def test_transit_not_mixed_into_c1(self):
        with self.assertRaisesRegex(FeedbackError, "NOT_ADMITTED"):
            boundary(last_selection={**SELECTION, "mode": "CORNER_TRANSIT"})

    def test_nan_rejected_even_before_promotion(self):
        with self.assertRaises((ValueError, FeedbackError)):
            boundary(total_length_A=float("nan"))

    def test_supplied_geometry_length_checked(self):
        with self.assertRaisesRegex(FeedbackError, "OUTSIDE_ROUTE"):
            boundary(total_length_A=.5)

    def test_tampered_observation_rejected(self):
        rows = samples(); rows[-1]["validated_progress_A"] = 10.
        with self.assertRaises(ValueError):
            boundary(rows)

    def test_six_consecutive_boundaries_and_budget(self):
        state = supervisor()
        for i in range(6):
            state = consume(state, boundary(samples(i*5., 1+i*.5)))
            self.assertEqual(state["completed_segments"], i+1)
        self.assertEqual(state["status"], "BUDGET_COMPLETE")
        self.assertEqual(state["observed_progress_A"], 4.)
        with self.assertRaisesRegex(FeedbackError, "CANNOT_RESUME"):
            consume(state, boundary(samples(30., 4.)))

    def test_stale_state_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "STALE_SUPERVISOR"):
            consume(supervisor(), boundary(), expected_state_sha256="0"*64)

    def test_wrong_checkpoint_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "CHECKPOINT_CHAIN"):
            consume(supervisor(), boundary(), seed_checkpoint_sha256="0"*64)

    def test_time_gap_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "TIME_CHAIN"):
            consume(supervisor(), boundary(samples(5.)))

    def test_progress_reset_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "RESET_FORBIDDEN"):
            consume(supervisor(), boundary(samples(progress=2.)))

    def test_hold_is_terminal(self):
        state = consume(supervisor(), boundary(qualification=None))
        self.assertEqual(state["status"], "HOLD")
        with self.assertRaisesRegex(FeedbackError, "CANNOT_RESUME"):
            consume(state, boundary(samples(5., 1.5)))

    def test_invalid_start_is_retained_hold(self):
        state = consume(supervisor(), boundary(samples(invalid_index=0)))
        self.assertEqual(state["status"], "HOLD")

    def test_cannot_promote_hold_by_changing_state_only(self):
        b = unseal(boundary(qualification=None)); b.update(state="SELECTOR_INPUT_READY", reasons=[])
        with self.assertRaises(FeedbackError):
            consume(supervisor(), seal(b))

    def test_input_objects_are_not_mutated(self):
        state = supervisor(); old = copy.deepcopy(state)
        consume(state, boundary())
        self.assertEqual(state, old)

    def test_matched_pair(self):
        p = paired_contract(supervisor(), supervisor("fixed"))
        self.assertTrue(p["same_initial_checkpoint_and_protocol"])
        self.assertFalse(p["native_handoff_verified"])
        self.assertFalse(p["scientific_pass"])

    def test_unequal_exposure_is_not_complete_pair(self):
        p = paired_contract(consume(supervisor(), boundary()), supervisor("fixed"))
        self.assertFalse(p["equal_observed_exposure"])
        self.assertEqual(p["common_end_time_ps"], 0.)

    def test_different_checkpoint_pair_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "PAIR_START"):
            paired_contract(supervisor(), supervisor("fixed", checkpoint_sha256="2"*64))

    def test_different_native_execution_contract_pair_rejected(self):
        with self.assertRaisesRegex(FeedbackError, "PAIR_START"):
            paired_contract(supervisor(), supervisor("fixed", execution_contract_sha256="2"*64))


if __name__ == "__main__":
    unittest.main()
