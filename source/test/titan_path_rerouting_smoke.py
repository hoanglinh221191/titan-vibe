import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_rerouting import TrackedRoute,shared_mode,interval_consensus,local_reroute_reasons,sampler_phase
from titan_path_observation import ObservationError


class ReroutingTests(unittest.TestCase):
    def test_closer_branch_cannot_steal_state(self):
        r=TrackedRoute([[0,0,0],[4,0,0],[4,4,0]])
        p=r.project_tracked([2,1,0],5,1)
        self.assertEqual((p.s,p.segment),(5,1));self.assertIn('OFF_BRANCH_NEARER_PROJECTION',r.tracking['warnings'])
    def test_zero_step_idempotent(self):
        r=TrackedRoute([[0,0,0],[4,0,0],[4,4,0]])
        p=r.project_tracked([2,1,0],5,1);q=r.project_tracked([2,1,0],p.s,p.segment)
        self.assertEqual(p,q)
    def test_adjacent_forward_transition(self):
        r=TrackedRoute([[0,0,0],[4,0,0],[4,4,0]])
        p=r.project_tracked([4.2,.3,0],3.9,0)
        self.assertEqual(p.segment,1);self.assertAlmostEqual(p.s,4.3)
        self.assertEqual(p,r.project_tracked([4.2,.3,0],p.s,p.segment))
    def test_adjacent_reverse_transition(self):
        r=TrackedRoute([[0,0,0],[4,0,0],[4,4,0]])
        p=r.project_tracked([3.5,-.1,0],4.1,1);self.assertEqual(p.segment,0)
        self.assertEqual(p,r.project_tracked([3.5,-.1,0],p.s,p.segment))
    def test_bad_branch(self):
        with self.assertRaises(ObservationError):TrackedRoute([[0,0,0],[4,0,0]]).project_tracked([1,0,0],1,2)
    def test_branch_progress_mismatch(self):
        with self.assertRaises(ObservationError):TrackedRoute([[0,0,0],[4,0,0],[4,4,0]]).project_tracked([1,0,0],1,1)
    def test_nonlocal_jump(self):
        with self.assertRaises(ObservationError):TrackedRoute([[0,0,0],[10,0,0]]).project_tracked([8,0,0],1,0)
    def test_recovery_threshold_straddle(self):
        self.assertEqual(shared_mode([2.999,2.996,3.001],'Recovery'),'Recovery')
    def test_all_frames_exit(self):self.assertEqual(shared_mode([2.9]*3,'Recovery'),'Advance')
    def test_advance_threshold_straddle(self):self.assertEqual(shared_mode([3.99,4.01,4.02],'Advance'),'Advance')
    def test_all_frames_enter(self):self.assertEqual(shared_mode([4.1]*3,'Advance'),'Recovery')
    def test_nan_mode(self):
        with self.assertRaises(ValueError):shared_mode([np.nan]*3,'Advance')
    def test_interval_overlap(self):
        self.assertTrue(interval_consensus([[9.75,10.46],[9.77,10.43],[9.73,10.50]])['compatible'])
    def test_point_float_tolerance(self):
        self.assertTrue(interval_consensus([[6.,6.],[6.+1e-14,6.+1e-14],[6.-1e-14,6.-1e-14]])['compatible'])
    def test_separated_intervals(self):self.assertFalse(interval_consensus([[1,1.1],[2,2.1],[3,3.1]])['compatible'])
    def test_invalid_intervals(self):
        with self.assertRaises(ValueError):interval_consensus([[2,1]]*3)
    def test_bounded_reroute_domain(self):
        self.assertFalse(local_reroute_reasons(1.8,40,[-.1]*3,[.1]*3,[4.3]*3,[4.2]*3))
    def test_reroute_cannot_ignore_clearance(self):
        self.assertIn('IMMEDIATE_CLEARANCE',local_reroute_reasons(1.4,40,[0]*3,[.1]*3,[4.3]*3,[4.2]*3))
    def test_reroute_not_unbounded(self):
        self.assertIn('REROUTE_EXCURSION_LIMIT',local_reroute_reasons(2,40,[0]*3,[.1]*3,[5]*3,[4.2]*3))
    def test_reroute_goal_required(self):
        self.assertIn('REROUTE_NO_GOAL_APPROACH',local_reroute_reasons(2,40,[0]*3,[0]*3,[4.3]*3,[4.2]*3))
    def test_random_idempotence(self):
        r=TrackedRoute([[0,0,0],[4,0,0],[4,4,0]]);rng=np.random.default_rng(81);checked=0
        for _ in range(1000):
            point=rng.uniform([-1,-1,-.2],[5,5,.2]);branch=int(rng.integers(2));s=2 if branch==0 else 6
            try:p=r.project_tracked(point,s,branch)
            except ObservationError:continue
            self.assertEqual(p,r.project_tracked(point,p.s,p.segment));checked+=1
        self.assertGreater(checked,500)
    def test_sampler_reproducible(self):
        self.assertEqual(sampler_phase(23,7,'reroute'),sampler_phase(23,7,'reroute'))
        self.assertEqual(sampler_phase(0,7,'reroute'),0.)
    def test_sampler_seed_domains_differ(self):
        self.assertEqual(len({sampler_phase(s,0,'reroute') for s in (11,23,37,53,71)}),5)
        self.assertNotEqual(sampler_phase(23,0,'reroute'),sampler_phase(23,1,'reroute'))
    def test_sampler_bounds(self):
        with self.assertRaises(ValueError):sampler_phase(-1,0,'reroute')


if __name__=='__main__':unittest.main()
