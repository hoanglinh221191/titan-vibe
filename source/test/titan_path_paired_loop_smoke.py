"""Orchestration safety checks with mocked workers; never MD."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_paired_loop as p
from titan_path_feedback import FeedbackError

class Smoke(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(prefix='PAIR-SYNTHETIC-');self.addCleanup(t.cleanup);self.root=Path(t.name)
    def test_no_local_md(self):
        with patch.dict(os.environ,{},clear=True):
            with self.assertRaisesRegex(FeedbackError,'MD_REQUIRES_BATCH'):p.run_one(self.root,{}, {},self.root,{})
    def test_unrecorded_batch_job_rejected(self):
        (self.root/'generated').mkdir();(self.root/'generated/SUBMISSION.json').write_text('{"array_job":"100"}')
        with patch.dict(os.environ,{'SLURM_JOB_PARTITION':'batch','SLURM_ARRAY_JOB_ID':'200'},clear=True):
            with self.assertRaisesRegex(FeedbackError,'UNRECORDED'):p.run_one(self.root,{}, {},self.root,{})
    def test_afterany_dependencies_and_single_submission(self):
        calls=[]
        def submit(command,**kw):
            calls.append(command)
            class R: returncode=0;stdout=str(100+len(calls))+'\n';stderr=''
            return R()
        with patch.object(p,'verify',return_value={}),patch.object(p.subprocess,'run',side_effect=submit):
            p.submit(self.root)
            with self.assertRaises(FileExistsError):p.submit(self.root)
        self.assertEqual(len(calls),3);self.assertIn('--dependency=afterany:101',calls[1]);self.assertIn('--dependency=afterany:102',calls[2])
        self.assertIn('--array=0-4%2',calls[1])
    def test_partial_submit_preserved_no_retry(self):
        class R:returncode=1;stdout='';stderr='SYNTHETIC'
        with patch.object(p,'verify',return_value={}),patch.object(p.subprocess,'run',return_value=R()) as call:
            with self.assertRaises(FeedbackError):p.submit(self.root)
            self.assertTrue((self.root/'generated/SUBMIT_build.json').exists());self.assertEqual(call.call_count,1)
    def test_missing_workers_create_diagnostic_bundle(self):
        cfg={'cells':[{'cell':'synthetic_'+str(i)} for i in range(5)]}
        with patch.object(p,'verify',return_value=cfg):self.assertEqual(p.collect(self.root),1)
        self.assertTrue((self.root/'downloads/FAILED.json').exists());self.assertFalse((self.root/'downloads/READY.json').exists())
    def test_all_method_stops_are_not_scientific_pass(self):
        cfg={'cells':[{'cell':'synthetic_'+str(i)} for i in range(5)]}
        for c in cfg['cells']:p.write(self.root/'outputs/pairs'/c['cell']/'PAIR_RESULT.json',{'cell':c['cell'],'status':'PAIR_STOP','scientific_pass':False})
        with patch.object(p,'verify',return_value=cfg):self.assertEqual(p.collect(self.root),0)
        self.assertFalse(p.read(self.root/'generated/COLLECTION.json')['scientific_pass'])
    def exercise_pair(self,behavior):
        cell={'cell':'synthetic','time_ps':0.}
        p.write(self.root/'generated/BUILD_IDENTITY.json',{'binary_sha256':{}})
        p.write(self.root/'generated/PREPARE_READY.json',{'status':'PASS_PREPARE_NO_MD','build_identity_sha256':p.sha(self.root/'generated/BUILD_IDENTITY.json')})
        p.write(self.root/'baselines/r2/synthetic/VETO.json',{'frames':[{'start_progress_A':3.}]*3})
        calls=[]
        def initial(root,c,arm):return {'state':{'arm':arm},'segments':[],'veto':{},'delta':[.25,0,0],'selection_sha':'a'*64}
        def prepared(root,c,context,path):path.mkdir(parents=True);return {}
        def proposal(root,c,context,path,adaptive):return {'accepted':not (behavior=='stop_second_round' and adaptive),'reason':'SYNTHETIC_STOP'}
        def run(root,c,context,path,receipt):
            arm=context['state']['arm'];calls.append(arm)
            if behavior=='fixed_failure' and arm=='fixed':raise RuntimeError('SYNTHETIC_FAILURE')
            status='HOLD' if behavior=='adaptive_hold' and arm=='adaptive' else 'READY'
            context['segments'].append({'observer_state':status,'observed_frame_progress_A':[3.1]*3})
            return status
        with patch.object(p,'verify',return_value={'cells':[cell]}),patch.object(p,'initial_context',side_effect=initial),\
             patch.object(p,'admit_direction',side_effect=lambda state,*rest:(state,{'SYNTHETIC':True})),\
             patch.object(p,'prepare_segment',side_effect=prepared),patch.object(p,'native_proposal',side_effect=proposal),patch.object(p,'run_one',side_effect=run),\
             patch.object(p,'validate_pending_pair',side_effect=RuntimeError('SYNTHETIC_TPR_GATE') if behavior=='bad_tpr_gate' else None):
            p.run_pair(self.root,0)
        return calls,p.read(self.root/'outputs/pairs/synthetic/PAIR_RESULT.json')
    def test_pair_stop_prevents_next_round_md(self):
        calls,result=self.exercise_pair('stop_second_round')
        self.assertEqual(calls,['adaptive','fixed']);self.assertEqual(result['common_validated_rounds'],1)
        self.assertEqual(result['status'],'PAIR_STOP')
    def test_hold_round_is_not_validated_exposure(self):
        calls,result=self.exercise_pair('adaptive_hold')
        self.assertEqual(calls,['adaptive','fixed']);self.assertEqual(result['common_completed_rounds'],1)
        self.assertEqual(result['common_validated_rounds'],0);self.assertEqual(result['status'],'PAIR_HOLD')
    def test_partial_pair_keeps_unpaired_exposure(self):
        calls,result=self.exercise_pair('fixed_failure')
        self.assertEqual(result['status'],'TECHNICAL_FAILURE');self.assertEqual(result['common_validated_rounds'],0)
        self.assertEqual(result['paired_metrics']['adaptive']['actual_completed_duration_ps'],5)
    def test_six_round_cap_is_hard(self):
        calls,result=self.exercise_pair('complete')
        self.assertEqual(len(calls),12);self.assertEqual(result['common_validated_rounds'],6)

    def test_pair_tpr_gate_blocks_both_md_arms(self):
        calls,result=self.exercise_pair('bad_tpr_gate')
        self.assertEqual(calls,[]);self.assertEqual(result['status'],'TECHNICAL_FAILURE')
        self.assertIn('SYNTHETIC_TPR_GATE',result['error'])

if __name__=='__main__':unittest.main(verbosity=2)
