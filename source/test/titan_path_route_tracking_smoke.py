import sys,unittest
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_route_tracking import BoundedRoute
from titan_path_observation import Polyline,Policy,ObservationError


class TrackingTests(unittest.TestCase):
    def test_adjacent_uncertainty_retained(self):
        r=BoundedRoute([[-2,0,0],[0,0,0],[0,2,0]])
        p=r.continuity([-.2,.2,0],2,1)
        self.assertFalse(r.rivals([-.2,.2,0],p,2,Policy()))
        self.assertTrue(r.tracking['accepted_adjacent_ambiguity'])
        self.assertAlmostEqual(r.tracking['progress_interval_A'][1]-r.tracking['progress_interval_A'][0],.4)
    def test_strict_still_flags(self):
        r=Polyline([[-2,0,0],[0,0,0],[0,2,0]]);p=r.continuity([-.2,.2,0],2,1)
        self.assertTrue(r.rivals([-.2,.2,0],p,2,Policy()))
    def test_nonadjacent_rival_not_waived(self):
        r=BoundedRoute([[0,0,0],[2,0,0],[2,.2,0],[0,.2,0]])
        p=r.continuity([1.7,.1,0],2.1,1)
        self.assertTrue(r.rivals([1.7,.1,0],p,2.1,Policy()))
    def test_window_recovery(self):
        r=BoundedRoute([[0,0,0],[10,0,0]]);p=r.continuity([3,0,0],1.5,1)
        self.assertEqual(p.s,3);self.assertTrue(r.tracking['expanded_window'])
    def test_large_jump_stays_stopped(self):
        r=BoundedRoute([[0,0,0],[10,0,0]])
        with self.assertRaises(ObservationError):r.continuity([5,0,0],1,1)
    def test_strict_success_unchanged(self):
        r=BoundedRoute([[0,0,0],[10,0,0]])
        self.assertEqual(r.continuity([2,1,0],1.5,1),Polyline([[0,0,0],[10,0,0]]).continuity([2,1,0],1.5,1))
        self.assertFalse(r.tracking['expanded_window'])
    def test_observer_opt_in_retains_interval_uncertainty(self):
        from titan_path_multiframe_feedback_smoke import fixture,CORE
        from titan_path_multiframe_feedback import observe_segment
        from titan_path_feedback import seal,unseal
        state,raw,whole,pull=fixture();body=unseal(state);body['progress_A']=[1.8]*3;body['distance_A']=[.2]*3;state=seal(body)
        for a,b in zip(raw,whole):a['xyz'][2603:]=[-.2,.2,0];b['xyz'][2603:]=[-.2,.2,0]
        pull[:,9:12]=[-.02,.02,0]
        args=dict(core=list(range(81)),reference=CORE,waypoints=[[-2,0,0],[0,0,0],[0,2,0]],ligand_masses=[12.]*12,
                  source_checkpoint_sha256='b'*64,endpoint_checkpoint_sha256='c'*64,trajectory_sha256='d'*64,endpoint_snapshot_sha256='e'*64)
        strict=observe_segment(state,raw,whole,pull,**args)
        bounded=observe_segment(state,raw,whole,pull,**args,route_tracking='bounded-adjacent-route-tracking-v1')
        self.assertEqual(strict[0]['status'],'HOLD');self.assertEqual(bounded[0]['status'],'READY')
        self.assertFalse(bounded[0]['global_frame_qualified'])
        self.assertTrue(all(f['projection_tracking']['accepted_adjacent_ambiguity'] for r in bounded[2] for f in r['frames']))


if __name__=='__main__':unittest.main()
