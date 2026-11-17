"""Offline toy invariants, not physical validation or MD permission."""
import math
from pathlib import Path
import sys
import unittest
from dataclasses import replace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_corner_reward_prototype import CornerPolicy,Polyline,cone_candidates,lookahead_goal,select_geometry,swept_clearance,virtual_episode


class PrototypeTests(unittest.TestCase):
    def setUp(self):
        self.p=CornerPolicy();self.line=Polyline([[-2,0,0],[0,0,0],[0,3,0]])
        self.routes=[self.line]*3;self.lig=np.array([[.2,-.7,0.]])
        self.rec=np.array([[0,0,10.]]);self.s=[2.]*3;self.old=[1.,0,0]
        self.goals=[lookahead_goal(r,2.,self.p) for r in self.routes]

    def choose(self,transition=True,rec=None):
        return select_geometry(self.rec if rec is None else rec,self.lig,[1.],self.routes,self.s,self.old,self.goals,self.p,transition)

    def test_reward_does_not_bypass_strict_gate(self):self.assertIsNone(self.choose(False)['selected'])
    def test_corner_gate_has_geometric_candidate(self):self.assertIsNotNone(self.choose()['selected'])
    def test_no_md_authority(self):self.assertFalse(self.choose()['md_authorized'])
    def test_energy_explicitly_unevaluated(self):self.assertFalse(self.choose()['energy_filter_evaluated'])
    def test_no_collision_bypass(self):self.assertIsNone(self.choose(rec=np.array([[.2,-.7,0.]]))['selected'])
    def test_probe_and_cone(self):
        for d in cone_candidates(self.old,self.goals,self.lig[0],[[0,1,0]],self.p):
            self.assertAlmostEqual(np.linalg.norm(d),.25,places=12)
            self.assertLessEqual(math.degrees(math.acos(np.clip(d[0]/.25,-1,1))),20+1e-8)
    def test_fixed_goal_at_outgoing_lookahead(self):np.testing.assert_allclose(self.goals[0],[0,1,0],atol=1e-14)
    def test_nonvertex_entry_rejected(self):
        with self.assertRaises(ValueError):lookahead_goal(self.line,2.1,self.p)
    def test_terminal_entry_rejected(self):
        with self.assertRaises(ValueError):lookahead_goal(self.line,5.,self.p)
    def test_collision_on_interior_sweep(self):
        self.assertAlmostEqual(swept_clearance([[.125,0,0]],[[0,0,0]],[.25,0,0]),0)
    def test_budget_bounded(self):
        result=virtual_episode(self.rec,self.lig,[1.],self.routes,self.s,self.old,self.p)
        self.assertLessEqual(len(result['steps']),8);self.assertFalse(result['md_executed'])
        for step in result['steps']:
            if step['selected']:
                self.assertGreater(min(step['selected']['goal_approach_A']),1e-6)
                self.assertGreaterEqual(min(step['selected']['delta_progress_A']),-1e-6)
    def test_input_arrays_unchanged(self):
        before=self.lig.copy();virtual_episode(self.rec,self.lig,[1.],self.routes,self.s,self.old,self.p)
        np.testing.assert_array_equal(before,self.lig)
    def test_bad_masses_rejected(self):
        with self.assertRaises(ValueError):select_geometry(self.rec,self.lig,[0],self.routes,self.s,self.old,self.goals,self.p,True)
    def test_three_frames_required(self):
        with self.assertRaises(ValueError):select_geometry(self.rec,self.lig,[1.],self.routes[:2],self.s,self.old,self.goals,self.p,True)
    def test_safety_change_rejected(self):
        with self.assertRaises(ValueError):replace(self.p,clearance_A=1.)
    def test_budget_extension_rejected(self):
        with self.assertRaises(ValueError):replace(self.p,max_virtual_steps=100)


if __name__=='__main__':unittest.main()
