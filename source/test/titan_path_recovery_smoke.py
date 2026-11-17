import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_recovery import RecoveryPolicy, immediate_reasons, start_excursion, observe_excursion
from titan_path_multiframe_feedback import frame_mode_reason


class RecoveryTests(unittest.TestCase):
    def test_absolute_bound_not_clearance(self):
        self.assertEqual(immediate_reasons(1.7,[0]*3,[4.3]*3),['IMMEDIATE_DISTANCE_BOUND'])
    def test_recovery_can_descend_outside_corridor(self):
        self.assertEqual(immediate_reasons(1.7,[0]*3,[4.1]*3,recovery=True,start_distance=[4.2]*3),[])
    def test_recovery_cannot_move_outward(self):
        self.assertIn('RECOVERY_NO_DISTANCE_DESCENT',immediate_reasons(1.7,[0]*3,[4.3]*3,recovery=True,start_distance=[4.2]*3))
    def test_clearance_still_hard(self):
        self.assertIn('IMMEDIATE_CLEARANCE',immediate_reasons(1.4,[0]*3,[4.1]*3,recovery=True,start_distance=[4.2]*3))
    def test_bounded_progress_loss(self):
        self.assertIn('RECOVERY_PROGRESS_LOSS',immediate_reasons(2.,[-.6,0,0],[4.1]*3,recovery=True,start_distance=[4.2]*3))
    def test_missing_projection(self):
        self.assertEqual(immediate_reasons(2.,None,None,projection_error=True),['IMMEDIATE_PROJECTION_UNAVAILABLE'])
    def test_separate_reasons(self):
        self.assertEqual(set(immediate_reasons(1.,[-1]*3,[5]*3,ambiguous=True)),
                         {'IMMEDIATE_CLEARANCE','IMMEDIATE_BACKTRACK','IMMEDIATE_DISTANCE_BOUND','IMMEDIATE_PROJECTION_AMBIGUOUS'})
    def test_finite(self):
        with self.assertRaises(ValueError): immediate_reasons(2.,[0]*3,[np.nan]*3)
    def test_advance_boolean_parity(self):
        rng=np.random.default_rng(112)
        for _ in range(1000):
            clearance=float(rng.uniform(1,3)); ds=rng.uniform(-.2,.2,3); d=rng.uniform(0,5,3); ambiguous=bool(rng.integers(2))
            old=clearance>=1.5 and min(ds)>=-1e-6 and max(d)<=4 and not ambiguous
            self.assertEqual(not immediate_reasons(clearance,ds,d,ambiguous=ambiguous),old)
    def test_outward_budget_not_reset(self):
        s=start_excursion([4.2]*3,[6.]*3)
        for x in (4.25,4.3,4.35): s=observe_excursion(s,[x]*3,[6.]*3)
        self.assertEqual(s['status'],'ACTIVE')
        self.assertEqual(observe_excursion(s,[4.4]*3,[6.]*3)['status'],'STOP_OUTWARD_BUDGET')
        self.assertEqual(s['anchor_distance_A'],[4.2]*3)
    def test_distance_cap(self):
        self.assertEqual(observe_excursion(start_excursion([4.2]*3,[6.]*3),[4.71]*3,[6.]*3)['status'],'STOP_EXCURSION_DISTANCE')
    def test_return(self):
        self.assertEqual(observe_excursion(start_excursion([4.2]*3,[6.]*3),[3.99]*3,[6.]*3)['status'],'RETURN_OBSERVED')
    def test_deadline(self):
        s=start_excursion([4.2]*3,[6.]*3)
        for _ in range(10): s=observe_excursion(s,[4.2]*3,[6.]*3)
        self.assertEqual(s['status'],'STOP_RETURN_DEADLINE')
        with self.assertRaises(ValueError): observe_excursion(s,[4.]*3,[6.]*3)
    def test_progress_budget(self):
        self.assertEqual(observe_excursion(start_excursion([4.2]*3,[6.]*3),[4.1]*3,[5.4]*3)['status'],'STOP_EXCURSION_PROGRESS_LOSS')
    def test_modes_missing_not_disagreement(self):
        self.assertEqual(frame_mode_reason([{'mode':None}]*3),'FRAME_PROJECTION_UNAVAILABLE')
        self.assertEqual(frame_mode_reason([{'mode':'Advance'},{'mode':'Recovery'},{'mode':'Advance'}]),'FRAME_MODE_DISAGREEMENT')
        self.assertIsNone(frame_mode_reason([{'mode':'Advance'}]*3))


if __name__=='__main__':unittest.main()
