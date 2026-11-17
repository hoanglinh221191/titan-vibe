"""Pure coordinator tests: no energy executable or MD invocation."""
import sys,unittest
from pathlib import Path
from dataclasses import replace
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_short_path_planner import SearchPolicy,plan_short_path


class Tests(unittest.TestCase):
    def rank(self,records,width,terminal):return sorted(range(len(records)),key=lambda i:records[i]['energy_key'])[:width]
    def expand(self,state):
        depth=len(state.steps)+1
        return [{'delta_A':[.25,0,0],'progress_A':[0 if depth<3 else .25]*3,'route_distance_A':[1]*3,'clearance_A':2.}]
    def runplan(self,expand=None,score=None,rank=None,policy=None):
        return plan_short_path([0]*3,[1]*3,[1,0,0],expand or self.expand,score or (lambda nodes:[1.]*len(nodes)),rank or self.rank,policy or SearchPolicy())
    def test_flat_prefix_can_lead_to_terminal(self):self.assertEqual(self.runplan()['status'],'OFFLINE_PATH_PROPOSED')
    def test_first_step_stays_small(self):self.assertEqual(self.runplan()['first_displacement_A'],[.25,0,0])
    def test_no_md_authority(self):self.assertFalse(self.runplan()['md_authorized'])
    def test_no_completed_path_no_step(self):
        r=self.runplan(policy=SearchPolicy(horizon_steps=2));self.assertIsNone(r['first_displacement_A'])
    def test_budget_no_partial_batch(self):
        r=self.runplan(policy=SearchPolicy(energy_evaluation_budget=1));self.assertEqual(r['energy_evaluations'],1)
        self.assertEqual(r['search_stop_reason'],'ENERGY_BUDGET_EXHAUSTED')
    def test_empty_frontier(self):self.assertEqual(self.runplan(expand=lambda s:[])['search_stop_reason'],'FRONTIER_EXHAUSTED')
    def test_large_probe_rejected(self):
        def expand(s):r=self.expand(s);r[0]['delta_A']=[1,0,0];return r
        self.assertIsNone(self.runplan(expand=expand)['selected_path'])
    def test_bad_clearance_rejected(self):
        def expand(s):r=self.expand(s);r[0]['clearance_A']=1.49;return r
        self.assertIsNone(self.runplan(expand=expand)['selected_path'])
    def test_turn_rejected(self):
        def expand(s):r=self.expand(s);r[0]['delta_A']=[0,.25,0];return r
        self.assertIsNone(self.runplan(expand=expand)['selected_path'])
    def test_backtrack_rejected(self):
        def expand(s):r=self.expand(s);r[0]['progress_A']=[-.01]*3;return r
        self.assertIsNone(self.runplan(expand=expand)['selected_path'])
    def test_infinite_energy_rejected(self):
        with self.assertRaises(ValueError):self.runplan(score=lambda nodes:[float('inf')]*len(nodes))
    def test_wrong_energy_count_rejected(self):
        with self.assertRaises(ValueError):self.runplan(score=lambda nodes:[])
    def test_bad_rank_rejected(self):
        with self.assertRaises(ValueError):self.runplan(rank=lambda *args:[100])
    def test_beam_bound(self):
        with self.assertRaises(ValueError):SearchPolicy(beam_width=100)
    def test_safety_fixed(self):
        with self.assertRaises(ValueError):SearchPolicy(minimum_clearance_A=1.)
    def test_stage_expansion_bound(self):
        with self.assertRaises(ValueError):self.runplan(expand=lambda s:self.expand(s)*264)
    def test_beam_preserves_alternative_to_dead_end(self):
        import math
        def expand(s):
            if not s.steps:
                a=self.expand(s)[0];b=dict(a);b['delta_A']=[.25*math.cos(.1),.25*math.sin(.1),0]
                return [a,b]
            if s.steps[0]['delta_A'][1]==0:return []
            depth=len(s.steps)+1
            return [{'delta_A':list(s.steps[0]['delta_A']),'progress_A':[.25 if depth>=3 else 0]*3,'route_distance_A':[1]*3,'clearance_A':2.}]
        self.assertIsNone(self.runplan(expand=expand,policy=SearchPolicy(beam_width=1))['selected_path'])
        self.assertIsNotNone(self.runplan(expand=expand,policy=SearchPolicy(beam_width=2))['selected_path'])
    def test_initial_recovery_state_not_accepted_as_advance(self):
        with self.assertRaises(ValueError):plan_short_path([0]*3,[5]*3,[1,0,0],self.expand,lambda x:[1]*len(x),self.rank)
    def test_provider_cannot_mutate_retained_steps(self):
        def expand(s):
            if s.steps:s.steps[0]['delta_A'][0]=100
            return self.expand(s)
        self.assertEqual(self.runplan(expand=expand)['first_displacement_A'],[.25,0,0])


if __name__=='__main__':unittest.main()
