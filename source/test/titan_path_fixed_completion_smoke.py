"""Fixed first-round orchestration tests; all MD/submission calls mocked."""
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_fixed_completion as f
import titan_path_paired_loop as p
from titan_path_feedback import FeedbackError, seal


class Fixed(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(prefix='FIXED-SYNTHETIC-');self.addCleanup(t.cleanup);self.root=Path(t.name)

    def env(self,index=0):
        return {'SLURM_JOB_PARTITION':'batch','SLURM_ARRAY_JOB_ID':'102','SLURM_ARRAY_TASK_ID':str(index),
                'SLURM_JOB_NUM_NODES':'1','SLURM_NTASKS':'24','SLURM_CPUS_PER_TASK':'2'}

    def submission(self):f.write(self.root/'generated/SUBMISSION.json',{'array_job':'102'})

    def test_no_local_md(self):
        with patch.dict(os.environ,{},clear=True),self.assertRaisesRegex(ValueError,'MD_REQUIRES_BATCH'):f.allocation_guard(self.root,0)

    def test_wrong_indices(self):
        for i in (-1,5,None,True,'0'):
            with self.assertRaisesRegex(ValueError,'INDEX'):f.allocation_guard(self.root,i)

    def test_wrong_job_task_or_resources(self):
        self.submission()
        for key,value in [('SLURM_ARRAY_JOB_ID','103'),('SLURM_ARRAY_TASK_ID','1'),('SLURM_NTASKS','48'),('SLURM_CPUS_PER_TASK','1'),('SLURM_JOB_NUM_NODES','2')]:
            e=self.env();e[key]=value
            with patch.dict(os.environ,e,clear=True),self.assertRaises(ValueError):f.allocation_guard(self.root,0)

    def test_recorded_allocation(self):
        self.submission()
        with patch.dict(os.environ,self.env(),clear=True):f.allocation_guard(self.root,0)

    def test_submission_chain_and_duplicate_guard(self):
        calls=[]
        def sbatch(command,**kw):
            calls.append(command)
            return type('R',(),{'returncode':0,'stdout':str(100+len(calls))+'\n','stderr':''})()
        with patch.object(f,'verify'),patch.object(f,'sha',side_effect=lambda p:f.INPUT_SHA if str(p)=='old' else f.REAUDIT_SHA),patch.object(f.subprocess,'run',side_effect=sbatch):
            f.submit(self.root,Path('old'),Path('audit'))
            with self.assertRaises(FileExistsError):f.submit(self.root,Path('old'),Path('audit'))
        self.assertEqual(len(calls),3);self.assertIn('--array=0-4%2',calls[1])
        self.assertIn('--dependency=afterany:101',calls[1]);self.assertIn('--dependency=afterany:102',calls[2])
        self.assertNotIn('adaptive',' '.join(sum(calls,[])))

    def test_partial_submission_no_retry(self):
        result=type('R',(),{'returncode':1,'stdout':'','stderr':'SYNTHETIC_FAILURE'})()
        with patch.object(f,'verify'),patch.object(f,'sha',side_effect=lambda p:f.INPUT_SHA if str(p)=='old' else f.REAUDIT_SHA),patch.object(f.subprocess,'run',return_value=result) as run:
            with self.assertRaisesRegex(ValueError,'NO_RETRY'):f.submit(self.root,Path('old'),Path('audit'))
        self.assertEqual(run.call_count,1);self.assertTrue((self.root/'generated/SUBMIT_prepare.json').exists())

    def test_malformed_submission_receipt(self):
        result=type('R',(),{'returncode':0,'stdout':'not a job','stderr':''})()
        with patch.object(f,'verify'),patch.object(f,'sha',side_effect=lambda p:f.INPUT_SHA if str(p)=='old' else f.REAUDIT_SHA),patch.object(f.subprocess,'run',return_value=result):
            with self.assertRaisesRegex(ValueError,'UNCERTAIN'):f.submit(self.root,Path('old'),Path('audit'))

    def test_missing_prepare_creates_failure_without_md(self):
        self.submission()
        with patch.object(f,'verify'),patch.dict(os.environ,self.env(),clear=True),patch.object(p,'run_one',side_effect=AssertionError('must not execute')):
            self.assertEqual(f.run_fixed(self.root,0),1)
        r=f.read(self.root/'outputs/fixed'/f.CELLS[0]/'FIXED_RESULT.json')
        self.assertFalse(r['md_dispatch_attempted']);self.assertIsNone(r['actual_fixed_ps'])

    def test_local_compatibility_cannot_authorize_md(self):
        self.submission();f.write(self.root/'generated/PREPARE_READY.json',{'status':'PASS_LOCAL_PREPARE_COMPATIBILITY'})
        with patch.object(f,'verify'),patch.dict(os.environ,self.env(),clear=True),patch.object(p,'run_one',side_effect=AssertionError('must not execute')):
            self.assertEqual(f.run_fixed(self.root,0),1)

    def test_worker_cannot_repeat(self):
        self.submission();(self.root/'generated'/('worker_'+f.CELLS[0]+'.guard')).mkdir()
        with patch.object(f,'verify'),patch.dict(os.environ,self.env(),clear=True),self.assertRaises(FileExistsError):f.run_fixed(self.root,0)

    def boundary(self):
        return {'state':'OBSERVATION_VALID','start_time_ps':20,'end_time_ps':25,'source_checkpoint_sha256':'a'*64,
                'identity':{'route':'SYNTHETIC'},'frame_progress_A':[4,4.1,4.2],'reasons':[]}

    def test_valid_pair_is_not_adaptive_superiority(self):
        b=self.boundary();r=f.compare_pair(b,copy.deepcopy(b),[3,3,3])
        self.assertEqual(r['validated_common_ps'],5);self.assertFalse(r['adaptive_benefit_tested']);self.assertFalse(r['scientific_pass'])
        self.assertEqual(r['adaptive_minus_fixed_progress_A_by_frame'],[0,0,0])

    def test_hold_keeps_exposure_but_no_validated_pair(self):
        a=self.boundary();b=self.boundary();a.update(state='HOLD',frame_progress_A=None,reasons=['PROJECTION_AMBIGUOUS'])
        r=f.compare_pair(a,b,[3,3,3])
        self.assertEqual(r['completed_common_ps'],5);self.assertEqual(r['validated_common_ps'],0)
        self.assertIsNone(r['adaptive_minus_fixed_progress_A_by_frame']);self.assertFalse(r['further_md_authorized'])

    def test_pair_seed_and_time_mismatch(self):
        for key,value in [('source_checkpoint_sha256','b'*64),('start_time_ps',21),('identity',{'other':'route'})]:
            a=self.boundary();b=self.boundary();b[key]=value
            with self.assertRaisesRegex(ValueError,'PAIR_START_IDENTITY'):f.compare_pair(a,b,[3,3,3])

    def test_nonfinite_progress_rejected(self):
        a=self.boundary();a['frame_progress_A'][1]=float('nan')
        with self.assertRaisesRegex(ValueError,'FINITE'):f.compare_pair(a,self.boundary(),[3,3,3])

    def test_failure_inclusive_collection(self):
        (self.root/'generated').mkdir()
        with patch.object(f,'verify'):self.assertEqual(f.collect(self.root),1)
        self.assertTrue((self.root/'downloads/FAILED.json').exists())
        self.assertFalse((self.root/'downloads/READY.json').exists())
        self.assertIsNone(f.read(self.root/'generated/AGGREGATE.json')['retained_adaptive_ps'])

    def test_all_held_pairs_still_terminal_not_scientific(self):
        records=[{'cell':c,'status':'PAIR_HOLD','scientific_pass':False} for c in f.CELLS]
        self.assertEqual(f.package_results(self.root,records),0)
        r=f.read(self.root/'generated/AGGREGATE.json')
        self.assertEqual(r['technical_pairs'],5);self.assertEqual(r['validated_pairs'],0)
        self.assertEqual(r['validated_new_fixed_ps'],25);self.assertFalse(r['automatic_next_round'])

    def test_runtime_root_dispatch_uses_new_identity(self):
        runtime=self.root/'new';retained=self.root/'retained';dst=self.root/'new/outputs/prepared';dst.mkdir(parents=True)
        f.write(runtime/'generated/SUBMISSION.json',{'array_job':'102'});(dst/'segment.tpr').write_text('SYNTHETIC')
        context={'state':{'record_sha256':'state'},'checkpoint':retained/'seed.cpt'}
        receipt={'tpr_sha256':f.sha(dst/'segment.tpr'),'state_sha256':'state','source_sha256':{}}
        with patch.dict(os.environ,self.env(),clear=True),patch.object(p,'validate_receipt'),patch.object(p,'logged',side_effect=RuntimeError('STOP_BEFORE_MD')) as logged:
            with self.assertRaisesRegex(RuntimeError,'STOP_BEFORE_MD'):p.run_one(retained,{},context,dst,receipt,execution_root=runtime)
        self.assertIn(str(runtime/'scripts/gmx_dispatch.sh'),logged.call_args.args[0])
        self.assertIn('mdrun',logged.call_args.args[0]);self.assertIn(retained/'seed.cpt',logged.call_args.args[0])
        self.assertEqual(logged.call_count,1)

    def worker_fixture(self,outcome):
        self.submission();name=f.CELLS[0];retained=self.root/'generated/retained';reference=self.root/'generated/adaptive_reaudit'
        f.write(retained/'MANIFEST.json',{});f.write(reference/'RESULT_MANIFEST.json',{})
        f.write(retained/'inputs/REPLAY_PLAN.json',{'cells':[{'cell':name}]})
        dst=self.root/'outputs/fixed'/name/'segment_00/prepared';dst.mkdir(parents=True)
        (dst/'segment.tpr').write_text('SYNTHETIC TPR');receipt={'synthetic':True};f.write(dst/'PREPARED.json',receipt)
        f.write(self.root/'generated/PREPARE_READY.json',{'status':'PASS_FIXED_PREPARE_NO_MD','local_compatibility_only':False,
            'failed_bundle_sha256':f.INPUT_SHA,'reaudit_bundle_sha256':f.REAUDIT_SHA,'engine_sha256':f.ENGINE_SHA,
            'cells':{name:{'receipt_sha256':f.sha(dst/'PREPARED.json'),'tpr_sha256':f.sha(dst/'segment.tpr')}}})
        context={'state':{'arm':'fixed'},'segments':[]};calls=[]
        def once(root,cell,ctx,segment,proof,**kw):
            calls.append((ctx['state']['arm'],segment,kw))
            if outcome=='failure':raise RuntimeError('SYNTHETIC_AFTER_DISPATCH')
            ctx['segments'].append({'SYNTHETIC':True});f.write(segment/'RESULT.json',{'SYNTHETIC':True});return outcome
        with patch.object(f,'verify'),patch.dict(os.environ,self.env(),clear=True),patch.object(f,'check_engine'),\
             patch.object(f,'bound_receipt',return_value=(context,None,receipt)),patch.object(p,'run_one',side_effect=once),\
             patch.object(p,'native_proposal',side_effect=AssertionError('no direction search')):
            code=f.run_fixed(self.root,0)
        return code,calls,f.read(self.root/'outputs/fixed'/name/'FIXED_RESULT.json')

    def test_fixed_runs_exactly_one_segment(self):
        code,calls,r=self.worker_fixture('READY')
        self.assertEqual(code,0);self.assertEqual(len(calls),1);self.assertEqual(calls[0][0],'fixed')
        self.assertEqual(calls[0][2]['execution_root'],self.root);self.assertEqual(r['fixed_validated_ps'],5)

    def test_fixed_hold_does_not_extend(self):
        code,calls,r=self.worker_fixture('HOLD')
        self.assertEqual(code,0);self.assertEqual(len(calls),1);self.assertEqual(r['observer_state'],'HOLD')
        self.assertFalse(r['further_md_authorized'])

    def test_failed_dispatch_never_retries_or_claims_zero_execution(self):
        code,calls,r=self.worker_fixture('failure')
        self.assertEqual(code,1);self.assertEqual(len(calls),1);self.assertTrue(r['md_dispatch_attempted'])
        self.assertIsNone(r['actual_fixed_ps']);self.assertFalse(r['automatic_retry'])

    def collector_fixture(self,tamper=False):
        f.write(self.root/'generated/PREPARE_READY.json',{'status':'PASS_FIXED_PREPARE_NO_MD'})
        f.write(self.root/'generated/retained/MANIFEST.json',{})
        f.write(self.root/'generated/adaptive_reaudit/RESULT_MANIFEST.json',{})
        for cell in f.CELLS:
            a=self.boundary();b=self.boundary()
            if cell=='branch_b__cyan':a.update(state='HOLD',frame_progress_A=None,reasons=['PROJECTION_AMBIGUOUS'])
            a=seal(a);b=seal(b)
            f.write(self.root/'generated/adaptive_reaudit/generated/cells'/cell/'BOUNDARY.json',a)
            f.write(self.root/'generated/retained/baselines/r2'/cell/'VETO.json',{'frames':[{'start_progress_A':3}]*3})
            old=self.root/'generated/retained/outputs/pairs'/cell/'adaptive/segment_00/prepared/md/invocations/invocation_0000'
            old.mkdir(parents=True);(old/'pullf.part0002.xvg').write_text('20 1\n25 2\n')
            dst=self.root/'outputs/fixed'/cell/'segment_00/prepared';dst.mkdir(parents=True)
            (dst/'native.dat').write_text('SYNTHETIC_NATIVE')
            native={'start_time_ps':20,'end_time_ps':25,'md_executed':True,'output_map':{'SYNTHETIC':str(dst/'native.dat')},
                    'output_hashes':{'SYNTHETIC':f.sha(dst/'native.dat')},'observation_boundary_sha256':b['record_sha256'],
                    'peak_absolute_scalar_force_kj_mol_nm':2}
            f.write(dst/'RESULT.json',native);f.write(dst/'observed/BOUNDARY.json',b)
            f.write(self.root/'outputs/fixed'/cell/'FIXED_RESULT.json',{'status':'PASS_FIXED_FIRST_SEGMENT','result_sha256':f.sha(dst/'RESULT.json')})
        if tamper:(self.root/'outputs/fixed'/f.CELLS[0]/'segment_00/prepared/native.dat').write_text('TAMPERED')
        with patch.object(f,'verify'):code=f.collect(self.root)
        return code,f.read(self.root/'generated/AGGREGATE.json')

    def test_collector_matches_five_pairs_and_preserves_one_hold(self):
        code,r=self.collector_fixture()
        self.assertEqual(code,0);self.assertEqual(r['technical_pairs'],5)
        self.assertEqual(r['validated_pairs'],4);self.assertEqual(r['held_pairs'],1)
        self.assertFalse(r['scientific_pass']);self.assertFalse(r['adaptive_benefit_tested'])

    def test_collector_detects_native_hash_drift(self):
        code,r=self.collector_fixture(tamper=True)
        self.assertEqual(code,1);self.assertEqual(r['technical_pairs'],4)
        self.assertTrue((self.root/'downloads/FAILED.json').exists())


if __name__=='__main__':unittest.main(verbosity=2)
