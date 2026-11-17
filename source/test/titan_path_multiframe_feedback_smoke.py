"""Synthetic observer/state/handoff tests; never calls MD."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_multiframe_feedback import new_state, observe_segment, observed_native_input, admit_direction, make_handoff, mode_at
from titan_path_feedback import FeedbackError, seal, unseal

IDENTITY={'route_id':'SYNTHETIC',**{k:'a'*64 for k in ('route_sha256','reference_sha256','frame_definition_sha256','protocol_sha256')}}
CORE=np.array([[x,y,z] for x in range(3) for y in range(3) for z in range(9)],dtype=float)

def fixture():
    raw=[];whole=[];pull=[]
    for i in range(26):
        xyz=np.zeros((2615,3));xyz[:81]=CORE;xyz[2603:]=[3+.01*i,1,0]
        frame={'time_ps':i*.2,'xyz':xyz,'box':np.eye(3)*100,'ids':list(range(2615))}
        raw.append(frame);whole.append(copy.deepcopy(frame))
        com=xyz[2603]/10
        pull.append([i*.2,com[0],0,*com,0,0,0,*com])
    state=new_state(IDENTITY,'b'*64,0.,[3.]*3,[1.]*3,'Advance','adaptive')
    return state,raw,whole,np.array(pull)


class Smoke(unittest.TestCase):
    def test_rerouting_does_not_waive_bad_com(self):
        s,r,w,p=fixture();p[10,9]+=1
        state,boundary,_=self.observe((s,r,w,p),route_tracking='branch-aware-rerouting-v1')
        self.assertEqual(state['status'],'HOLD');self.assertIn('LIGAND_COM_NATIVE_MISMATCH',boundary['reasons'])
    def test_rerouting_does_not_waive_image_break(self):
        s,r,w,p=fixture();r[10]['xyz'][0,0]+=100
        state,boundary,_=self.observe((s,r,w,p),route_tracking='branch-aware-rerouting-v1')
        self.assertEqual(state['status'],'HOLD');self.assertIn('PRIMARY_IMAGE_DISCONTINUITY',boundary['reasons'])
    def test_rerouting_missing_projection_stays_hold(self):
        from titan_path_rerouting import TrackedRoute
        from titan_path_observation import ObservationError
        with patch.object(TrackedRoute,'project_tracked',side_effect=ObservationError('TRACKED_PROGRESS_JUMP_UNRESOLVED')):
            state,boundary,_=self.observe(route_tracking='branch-aware-rerouting-v1')
        self.assertEqual(state['status'],'HOLD');self.assertIn('FRAME_PROJECTION_UNAVAILABLE',boundary['reasons'])
    def test_missing_projection_is_not_mode_disagreement(self):
        from titan_path_observation import Polyline, ObservationError
        with patch.object(Polyline,'continuity',side_effect=ObservationError('PROJECTION_CONTINUITY_WINDOW_EXHAUSTED')):
            state,boundary,rows=self.observe()
        self.assertEqual(state['status'],'HOLD')
        self.assertIn('FRAME_PROJECTION_UNAVAILABLE',boundary['reasons'])
        self.assertNotIn('FRAME_MODE_DISAGREEMENT',boundary['reasons'])
        self.assertTrue(all(f['projection_error']=='PROJECTION_CONTINUITY_WINDOW_EXHAUSTED' for r in rows for f in r['frames']))
    def observe(self,f=None,**kw):
        s,r,w,p=fixture() if f is None else f
        return observe_segment(s,r,w,p,core=list(range(81)),reference=CORE,waypoints=[[0,0,0],[10,0,0]],
            ligand_masses=[12.]*12,source_checkpoint_sha256='b'*64,endpoint_checkpoint_sha256='c'*64,
            trajectory_sha256='d'*64,endpoint_snapshot_sha256='e'*64,**kw)
    def test_valid_finite_frame_contract_keeps_global_hold(self):
        s,b,rows=self.observe();self.assertEqual(s['status'],'READY');self.assertEqual(b['state'],'OBSERVATION_VALID')
        self.assertFalse(s['global_frame_qualified']);self.assertFalse(s['md_authorized']);self.assertEqual(len(rows),26)
    def test_measured_motion_not_probe(self):
        s,r,w,p=fixture()
        for i in range(26):
            r[i]['xyz'][2603:,0]=3+.02*i;w[i]['xyz'][2603:,0]=3+.02*i;p[i,9]=.3+.002*i
        a,_,_=self.observe((s,r,w,p));self.assertAlmostEqual(a['progress_A'][0],3.5)
    def test_missing_frame_is_error(self):
        s,r,w,p=fixture();del r[10];del w[10]
        with self.assertRaises(FeedbackError):self.observe((s,r,w,p))
    def test_native_xtc_time_quantization(self):
        s,r,w,p=fixture();r[1]['time_ps']-=3e-6;w[1]['time_ps']-=3e-6
        self.assertEqual(self.observe((s,r,w,p))[0]['status'],'READY')
    def test_wrong_pull_sample_time_is_error(self):
        s,r,w,p=fixture();p[1,0]+=.0009
        with self.assertRaises(FeedbackError):self.observe((s,r,w,p))
    def test_intermediate_image_break_never_reacquires(self):
        s,r,w,p=fixture();r[10]['xyz'][0,0]+=100
        a,b,rows=self.observe((s,r,w,p));self.assertEqual(a['status'],'HOLD');self.assertFalse(rows[-1]['chain_valid'])
        with self.assertRaises(FeedbackError):observed_native_input(a,10.)
    def test_nonprimary_image_change_is_diagnostic(self):
        s,r,w,p=fixture();r[10]['xyz'][100,0]+=100
        a,_,_=self.observe((s,r,w,p));self.assertEqual(a['status'],'READY')
    def test_bad_com_holds(self):
        s,r,w,p=fixture();p[10,9]+=1
        self.assertEqual(self.observe((s,r,w,p))[0]['status'],'HOLD')
    def test_seed_progress_reset_holds(self):
        s,r,w,p=fixture();b=unseal(s);b['progress_A']=[2.5]*3
        self.assertEqual(self.observe((seal(b),r,w,p))[0]['status'],'HOLD')
    def test_state_checksum_tampering_rejected(self):
        s,r,w,p=fixture();s['time_ps']=1
        with self.assertRaises(FeedbackError):self.observe((s,r,w,p))
    def test_native_input_is_not_forged_continue(self):
        s,_,_=self.observe();text=observed_native_input(s,10.)
        self.assertIn('TITAN_OBSERVED_BOUNDARY_V1',text);self.assertNotIn('CONTINUE',text)
    def test_handoff_reports_vector_jump(self):
        h,j=make_handoff([5,1,1.1,1,0,0,0,0,0,1,0,0],[1,0,0],[0,1,0])
        self.assertAlmostEqual(h.scalar_force_before_kj_mol_nm,h.scalar_force_after_kj_mol_nm)
        self.assertGreater(j['vector_force_jump_kj_mol_nm'],0);self.assertFalse(j['vector_force_continuity_claimed'])
    def test_budget_complete_is_terminal(self):
        s,r,w,p=fixture();b=unseal(s);b['segments']=5
        a,_,_=self.observe((seal(b),r,w,p));self.assertEqual(a['status'],'BUDGET_COMPLETE')
    def test_explicit_new_campaign_budget_does_not_inherit_six_rounds(self):
        from titan_path_feedback import LoopPolicy
        s,r,w,p=fixture();b=unseal(s);b['segments']=5
        a,_,_=self.observe((seal(b),r,w,p),policy=LoopPolicy(max_segments=100))
        self.assertEqual(a['status'],'READY')
        b['segments']=99
        a,_,_=self.observe((seal(b),r,w,p),policy=LoopPolicy(max_segments=100))
        self.assertEqual(a['status'],'BUDGET_COMPLETE')
        b['segments']=100
        with self.assertRaises(FeedbackError):self.observe((seal(b),r,w,p),policy=LoopPolicy(max_segments=100))
    def test_tangent_arm_is_recorded_without_global_qualification(self):
        s=new_state(IDENTITY,'b'*64,0.,[3.]*3,[1.]*3,'Advance','tangent')
        self.assertEqual(s['arm'],'tangent');self.assertFalse(s['global_frame_qualified'])
    def test_mode_hysteresis(self):
        self.assertEqual(mode_at(3.5,'Recovery'),'Recovery');self.assertEqual(mode_at(3.5,'Advance'),'Advance')
    def test_veto_stop_cannot_admit(self):
        s,_,_,_=fixture()
        with self.assertRaises(FeedbackError):admit_direction(s,{'md_authorized':False,'frames':[{}]*3,'accepted':False},[.25,0,0],'c'*64,'NATIVE_TITAN_SELECTION')

if __name__=='__main__':unittest.main(verbosity=2)
