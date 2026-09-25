"""NAV6-A focused synthetic regression; no MD or original-result writes."""
from pathlib import Path
import copy
import json
import math
import subprocess
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_observation import (ObservationError, Polyline, Policy, Tracker,
    anchor_quality, register, digest, validate_observation_chain, image_discontinuities)

ORACLE=None
if '--oracle' in sys.argv:
    i=sys.argv.index('--oracle'); ORACLE=Path(sys.argv[i+1]).resolve(); del sys.argv[i:i+2]
IDENTITY={'route_id':'synthetic', **{k:'a'*64 for k in ('route_sha256','reference_sha256','policy_sha256','trajectory_sha256','seed_sha256')}}
ANCHORS=np.array([[0.,0.,0.],[10.,0.,0.],[0.,5.,0.]])


class ObservationSmoke(unittest.TestCase):
    def tracker(self,start=0.,policy=Policy()): return Tracker(IDENTITY,start,policy)
    def sample(self,tracker,route,x,t,**kwargs): return tracker.observe(route,x,t,'b'*64,**kwargs)

    def test_straight_projection(self):
        p=Polyline([[0,0,0],[10,0,0]]).project([3,2,0]); self.assertEqual((p.s,p.d),(3,2))

    def test_large_total_small_steps(self):
        route=Polyline([[0,0,0],[10,0,0]]); tracker=self.tracker()
        rows=[self.sample(tracker,route,[x,0,0],i*.02) for i,x in enumerate(np.linspace(0,2.14,101))]
        self.assertTrue(rows[-1]['chain_valid']); self.assertAlmostEqual(rows[-1]['validated_progress_A'],2.14)
        validate_observation_chain(rows,IDENTITY)

    def test_single_large_step_is_clipped(self):
        route=Polyline([[0,0,0],[10,0,0]]); t=self.tracker()
        self.sample(t,route,[0,0,0],0); row=self.sample(t,route,[2.14,0,0],.2)
        self.assertIn('PROJECTION_CONTINUITY_WINDOW_EXHAUSTED',row['reasons'])
        self.assertIsNone(row['pairwise_local_progress_A']); self.assertIsNone(row['validated_progress_A'])

    def test_invalid_chain_never_recovers(self):
        route=Polyline([[0,0,0],[10,0,0]]); t=self.tracker()
        rows=[self.sample(t,route,[0,0,0],0),self.sample(t,route,[2,0,0],.2),self.sample(t,route,[2.1,0,0],.4)]
        self.assertTrue(rows[-1]['frame_valid']); self.assertFalse(rows[-1]['chain_valid'])
        self.assertIsNone(rows[-1]['validated_progress_A']); validate_observation_chain(rows,IDENTITY)

    def test_backward_is_observation_not_mode_error(self):
        route=Polyline([[0,0,0],[10,0,0]]); t=self.tracker(3)
        row=self.sample(t,route,[2.5,0,0],0)
        self.assertTrue(row['frame_valid']); self.assertEqual(row['validated_progress_A'],2.5)

    def test_real_end_is_not_artificial_clipping(self):
        p=Polyline([[0,0,0],[1,0,0]]).continuity([2,0,0],.5,1)
        self.assertEqual(p.s,1)

    def test_low_window_boundary_clipped(self):
        with self.assertRaisesRegex(ObservationError,'WINDOW_EXHAUSTED'):
            Polyline([[0,0,0],[10,0,0]]).continuity([1,0,0],3,1)

    def test_corner_vertex_is_one_minimum(self):
        route=Polyline([[0,0,0],[1,0,0],[1,1,0]])
        self.assertEqual(len(route.stationary_minima([1.2,-.2,0])),1)
        self.assertTrue(self.sample(self.tracker(.9),route,[1.2,-.2,0],0)['frame_valid'])

    def test_bend_has_no_false_endpoint_rival(self):
        route=Polyline([[0,0,0],[1,0,0],[1,1,0]])
        minima=route.stationary_minima([.8,-.1,0])
        self.assertEqual(len(minima),1); self.assertAlmostEqual(minima[0].s,.8)

    def test_near_self_approach_ambiguous(self):
        route=Polyline([[0,0,0],[.5,0,0],[.5,.02,0],[0,.02,0]])
        row=self.sample(self.tracker(.5),route,[.2,.01,0],0)
        self.assertIn('PROJECTION_AMBIGUOUS',row['reasons']); self.assertFalse(row['chain_valid'])

    def test_temporal_gap_separate_from_geometry(self):
        route=Polyline([[0,0,0],[10,0,0]]); t=self.tracker()
        self.sample(t,route,[0,0,0],0); row=self.sample(t,route,[.1,0,0],.4)
        self.assertIn('INSUFFICIENT_TEMPORAL_RESOLUTION',row['reasons'])

    def test_nonmonotone_time_rejected(self):
        route=Polyline([[0,0,0],[10,0,0]]); t=self.tracker()
        self.sample(t,route,[0,0,0],1)
        with self.assertRaisesRegex(ObservationError,'NONMONOTONIC'): self.sample(t,route,[0,0,0],1)

    def test_identity_mismatch(self):
        with self.assertRaisesRegex(ObservationError,'IDENTITY_MISMATCH'):
            self.sample(self.tracker(),Polyline([[0,0,0],[10,0,0]]),[0,0,0],0,identity={**IDENTITY,'route_id':'other'})

    def test_bad_hash_rejected(self):
        with self.assertRaises(ObservationError): Tracker({**IDENTITY,'seed_sha256':'bad'},0)

    def test_native_quality_error_poison_chain(self):
        row=self.sample(self.tracker(),Polyline([[0,0,0],[10,0,0]]),[0,0,0],0,quality_reasons=['NATIVE_COM_IMAGE_MISMATCH'])
        self.assertFalse(row['chain_valid'])

    def test_planar_three_anchors_valid(self):
        self.assertGreater(anchor_quality(ANCHORS)['sigma2_over_sigma1'],.1)
        r,t,rmsd,*_=register(ANCHORS,ANCHORS+1)
        np.testing.assert_allclose(r,np.eye(3),atol=1e-14); self.assertLess(rmsd,1e-13)

    def test_collinear_current_rejected(self):
        with self.assertRaisesRegex(ObservationError,'ILL_CONDITIONED'):
            register(ANCHORS,[[0,0,0],[1,0,0],[2,0,0]])

    def test_skinny_reference_rejected(self):
        ref=[[0,0,0],[1,0,0],[.5,1e-8,0]]
        with self.assertRaisesRegex(ObservationError,'ILL_CONDITIONED'): register(ref,ref)

    def test_zero_extent_rejected(self):
        with self.assertRaises(ObservationError): register(np.zeros((3,3)),ANCHORS)

    def test_invalid_coordinates(self):
        for invalid in ([[0,0,0]],[[0,0,0],[0,0,0]],[[0,0,0],[math.nan,0,0]]):
            with self.assertRaises(ObservationError): Polyline(invalid)

    def test_invalid_policy(self):
        with self.assertRaises(ObservationError): Policy(half_window=0)

    def test_unrelated_image_event_does_not_invalidate_primary(self):
        old=np.zeros((4,3)); new=old.copy(); new[3,0]=10
        check=image_discontinuities(new,old,np.eye(3)*10,[0,1,2])
        self.assertEqual(check['changed_atom_ids'],[3]); self.assertFalse(check['primary_changed'])

    def test_primary_image_event_invalidates_primary(self):
        old=np.zeros((4,3)); new=old.copy(); new[1,0]=10
        self.assertTrue(image_discontinuities(new,old,np.eye(3)*10,[0,1,2])['primary_changed'])

    def test_image_check_rejects_singular_box(self):
        with self.assertRaises(ObservationError): image_discontinuities(np.zeros((4,3)),np.zeros((4,3)),np.zeros((3,3)),[0])

    def test_rigid_motion_invariance(self):
        angle=.61; rot=np.array([[math.cos(angle),-math.sin(angle),0],[math.sin(angle),math.cos(angle),0],[0,0,1]])
        points=np.array([[0,0,0],[5,0,0],[5,5,0]],dtype=float); point=np.array([3.,1.,0.])
        before=Polyline(points).continuity(point,3,1)
        fit,trans,*_=register(ANCHORS,ANCHORS@rot.T+[10,-7,4])
        after=Polyline(points@fit.T+trans).continuity(point@rot.T+[10,-7,4],3,1)
        self.assertAlmostEqual(before.s,after.s); self.assertAlmostEqual(before.d,after.d)

    def test_tamper_is_rejected(self):
        row=self.sample(self.tracker(),Polyline([[0,0,0],[10,0,0]]),[0,0,0],0)
        row['global_progress_A']=.2
        with self.assertRaisesRegex(ObservationError,'INTEGRITY'): validate_observation_chain([row],IDENTITY)

    def test_no_continue_or_mode_emitted(self):
        row=self.sample(self.tracker(),Polyline([[0,0,0],[10,0,0]]),[0,0,0],0)
        self.assertNotIn('mode',row); self.assertNotIn('CONTINUE',json.dumps(row))

    def test_emitted_identity_does_not_alias_tracker(self):
        tracker=self.tracker(); route=Polyline([[0,0,0],[10,0,0]])
        row=self.sample(tracker,route,[0,0,0],0); row['identity']['route_id']='tampered'
        second=self.sample(tracker,route,[.1,0,0],.2)
        self.assertEqual(second['identity']['route_id'],'synthetic')

    def test_string_boolean_rejected_after_rehash(self):
        row=self.sample(self.tracker(),Polyline([[0,0,0],[10,0,0]]),[0,0,0],0)
        row['frame_valid']='true'; row.pop('observation_sha256'); row['observation_sha256']=digest(row)
        with self.assertRaisesRegex(ObservationError,'SCHEMA'): validate_observation_chain([row],IDENTITY)

    def test_validated_progress_cannot_be_replaced_after_rehash(self):
        row=self.sample(self.tracker(),Polyline([[0,0,0],[10,0,0]]),[0,0,0],0)
        row['validated_progress_A']=3; row.pop('observation_sha256'); row['observation_sha256']=digest(row)
        with self.assertRaisesRegex(ObservationError,'PROMOTION'): validate_observation_chain([row],IDENTITY)

    def test_rehashed_promotion_still_rejected(self):
        row=self.sample(self.tracker(),Polyline([[0,0,0],[10,0,0]]),[3,0,0],0)
        row['chain_valid']=True; row['validated_progress_A']=3
        row.pop('observation_sha256'); row['observation_sha256']=digest(row)
        with self.assertRaisesRegex(ObservationError,'PROMOTION'): validate_observation_chain([row],IDENTITY)

    @unittest.skipUnless(ORACLE,'C++ oracle path not provided')
    def test_cpp_geometry_and_clipping(self):
        payload=[]; expected=[]
        for points in ([[0,0,0],[10,0,0]],[[0,0,0],[1,0,0],[1,1,0]],[[0,0,0],[.5,0,0],[.5,.02,0],[0,.02,0]]):
            route=Polyline(points)
            for s0 in np.linspace(0,route.total,7):
                for point in ([0,0,0],[.2,.01,0],[1.2,-.2,0],[2.14,.3,0]):
                    values=[*ANCHORS.ravel(),*ANCHORS.ravel(),*np.ravel(points),*point,s0,1.0]
                    payload.append(f'3 {len(points)}\n'+' '.join(format(float(v),'.17g') for v in values)+'\n')
                    g=route.project(point,float(s0))
                    try: local=route.continuity(point,float(s0),1); kind='VALID'
                    except ObservationError: local=None; kind='CLIPPED'
                    expected.append((g,local,kind))
        run=subprocess.run([str(ORACLE)],input=''.join(payload),text=True,capture_output=True,check=True,timeout=30)
        self.assertEqual(len(run.stdout.splitlines()),len(expected))
        for line,(g,p,kind) in zip(run.stdout.splitlines(),expected):
            s,d,rmsd,seg,state,local_s,local_d,local_seg=line.split()
            self.assertAlmostEqual(float(s),g.s,places=8); self.assertAlmostEqual(float(d),g.d,places=8)
            self.assertEqual(state,kind)
            if p:
                self.assertAlmostEqual(float(local_s),p.s,places=8); self.assertAlmostEqual(float(local_d),p.d,places=8)


if __name__=='__main__': unittest.main()
