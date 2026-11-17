"""Synthetic workflow/guard regressions. All subprocess/MD actions mocked."""
from contextlib import ExitStack
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_benchmark_runtime as runtime
from titan_path_benchmark_target import persistence_step


class WorkerTests(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(prefix='B1-WORKER-SYNTHETIC-');self.addCleanup(t.cleanup);self.root=Path(t.name)
        self.cell={'cell':'synthetic','trajectory_id':'synthetic','method':'adaptive','route_label':'blue','block':0,'source_branch':'branch_a'}
        self.initial=self.root/'generated/initial/synthetic';(self.initial/'prepared').mkdir(parents=True)
        self.context={'segments':[],'state':{'record_sha256':'a'*64},'delta':[.25,0,0]}
        self.receipt={'tpr_sha256':'a'*64,'rng_lock':{'ld_seed':123}}

    def simulate(self,*,stop_at=None,hold_at=None,target_at=None,md_error_at=None,ready=True):
        calls={'md':0,'selection':0}
        def read(path):
            name=Path(path).name
            values={'SUBMISSION.json':{'array_job':'123'},'PREPARE_READY.json':{'status':'PASS_PREPARE_NO_MD' if ready else 'FAILED','build_identity_sha256':'a'*64,'initial_artifacts':{}},
                    'BUILD_IDENTITY.json':{'binary_sha256':{}},'TARGET.json':{},'FRAME_DEFINITION.json':{'selected_source_ids':[],'reference_xyz_A':[]},
                    'SELECTION_RESULT.json':{'accepted':stop_at!=0},'PREPARED.json':self.receipt,'BOUNDARY.json':{'reasons':[]}}
            return copy.deepcopy(values[name])
        def proposal(*args,**kwargs):
            calls['selection']+=1
            return {'accepted':stop_at!=calls['selection']}
        def prepare(root,cell,context,segment):segment.mkdir();return self.receipt
        def md(root,cell,context,segment,receipt,**kwargs):
            self.assertEqual(kwargs['observation_policy'].max_segments,100)
            calls['md']+=1
            if md_error_at==calls['md']:raise ValueError('SYNTHETIC_NATIVE_FAILURE')
            observed=segment/'observed';observed.mkdir()
            (observed/'whole.gro').write_text('SYNTHETIC\n2615\n'+'\n'*2616)
            context['segments'].append({'output_map':{'xtc':'synthetic.xtc'},'end_time_ps':1000+5*calls['md']})
            return 'HOLD' if calls['md']==hold_at else 'READY'
        def persistence(*args,**kwargs):return {'last_time_ps':1000+5*calls['md'],'target_reached':calls['md']==target_at}
        with ExitStack() as stack:
            stack.enter_context(patch.dict(os.environ,{'SLURM_ARRAY_JOB_ID':'123','SLURM_ARRAY_TASK_ID':'0','SLURM_JOB_PARTITION':'batch'}))
            for owner,name,value in [(runtime,'verify',lambda _: {'cells':[self.cell]}),(runtime,'read',read),
                    (runtime,'sha',lambda _: 'a'*64),(runtime,'load_context',lambda _: copy.deepcopy(self.context)),
                    (runtime,'continuation_gate',lambda *a:None),(runtime,'read_gro',lambda _: [{'time_ps':1000+5*calls['md'],'xyz':[],'box':[]}]),
                    (runtime,'measure_target',lambda *a:{'condition_met':False}),(runtime,'persistence_step',persistence),
                    (runtime,'logged',lambda *a,**k:''),(runtime.paired,'native_proposal',proposal),
                    (runtime.paired,'prepare_segment',prepare),(runtime.paired,'run_one',md)]:
                stack.enter_context(patch.object(owner,name,side_effect=value))
            code=runtime.run_trajectory(self.root,0)
        report=json.loads((self.root/'outputs/trajectories/synthetic/TRAJECTORY_RESULT.json').read_text())
        return code,report,calls

    def test_exact_100_segment_budget_not_six(self):
        code,report,calls=self.simulate();self.assertEqual(code,0);self.assertEqual(calls['md'],100);self.assertEqual(report['status'],'BUDGET_EXHAUSTED')

    def test_initial_method_stop_runs_no_md(self):
        code,report,calls=self.simulate(stop_at=0);self.assertEqual(calls['md'],0);self.assertEqual(report['status'],'METHOD_STOP')

    def test_later_method_stop_no_next_md(self):
        code,report,calls=self.simulate(stop_at=2);self.assertEqual(calls['md'],2);self.assertEqual(report['status'],'METHOD_STOP')

    def test_observer_hold_stops_only_owned_trajectory(self):
        code,report,calls=self.simulate(hold_at=3);self.assertEqual(calls['md'],3);self.assertEqual(report['status'],'OBSERVER_HOLD')
        self.assertEqual([p.name for p in (self.root/'outputs/trajectories').iterdir()],['synthetic'])

    def test_target_stops_after_committed_segment(self):
        code,report,calls=self.simulate(target_at=4);self.assertEqual(calls['md'],4);self.assertEqual(report['status'],'TARGET_REACHED')

    def test_native_failure_retains_completed_exposure_no_retry(self):
        code,report,calls=self.simulate(md_error_at=2);self.assertEqual(code,1);self.assertEqual(calls['md'],2)
        self.assertEqual(report['completed_segments'],1);self.assertEqual(report['status'],'TECHNICAL_FAILURE');self.assertFalse(report['auto_retry'])

    def test_failed_prepare_prevents_md(self):
        code,report,calls=self.simulate(ready=False);self.assertEqual(calls['md'],0);self.assertEqual(report['status'],'TECHNICAL_FAILURE')


class SubmissionCollectionTests(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(prefix='B1-CHAIN-SYNTHETIC-');self.addCleanup(t.cleanup);self.root=Path(t.name)

    def test_submit_chain_afterany_and_30_indices(self):
        results=[SimpleNamespace(returncode=0,stdout=str(n),stderr='') for n in (101,102,103)]
        with patch.object(runtime,'verify'),patch.object(runtime.subprocess,'run',side_effect=results) as call:
            runtime.submit(self.root)
        commands=[x.args[0] for x in call.call_args_list]
        self.assertIn('--array=0-29%4',commands[1]);self.assertIn('--dependency=afterany:101',commands[1]);self.assertIn('--dependency=afterany:102',commands[2])
        with patch.object(runtime,'verify'),patch.object(runtime.subprocess,'run') as call:
            with self.assertRaises(FileExistsError):runtime.submit(self.root)
            self.assertEqual(call.call_count,0)

    def test_partial_submission_kept_no_retry(self):
        results=[SimpleNamespace(returncode=0,stdout='101',stderr=''),SimpleNamespace(returncode=1,stdout='',stderr='synthetic')]
        with patch.object(runtime,'verify'),patch.object(runtime.subprocess,'run',side_effect=results) as call:
            with self.assertRaises(ValueError):runtime.submit(self.root)
        self.assertEqual(call.call_count,2);self.assertTrue((self.root/'generated/SUBMIT_build.json').is_file());self.assertFalse((self.root/'generated/SUBMISSION.json').exists())

    def fixture(self):
        cells=[{'cell':m,'method':m} for m in ('adaptive','tangent','fixed')]
        for folder in ('generated','outputs','tryton_logs'): (self.root/folder).mkdir()
        for name in ('PLAN.json','PROTOCOL.json','RUNTIME_ROUTES.json','PACKAGE_LOCK.json'):(self.root/name).write_text('{}')
        (self.root/'FILES.sha256').write_text('SYNTHETIC')
        (self.root/'generated/SUBMISSION.json').write_text('{"collect_job":"103"}')
        return cells

    def test_missing_workers_still_make_hash_verified_diagnostic_download(self):
        cells=self.fixture()
        with patch.object(runtime,'verify',return_value={'cells':cells}),patch.dict(os.environ,{'SLURM_JOB_ID':'103'}):runtime.collect(self.root)
        ready=json.loads((self.root/'downloads/READY.json').read_text());self.assertEqual(ready['status'],'FAILED_WITH_DIAGNOSTICS')
        runtime.audit(self.root)
        report=json.loads((self.root/'generated/COLLECTION.json').read_text());self.assertEqual(len(report['trajectories']),3)
        self.assertTrue(all(r['status']=='TECHNICAL_FAILURE' for r in report['trajectories']))
        bundle=self.root/'downloads'/ready['bundle'];bundle.write_bytes(bundle.read_bytes()+b'corruption')
        with self.assertRaises(ValueError):runtime.audit(self.root)

    def test_target_first_hit_is_latched(self):
        state={}
        for i in range(26):state=persistence_step(state,i*.2,True)
        self.assertTrue(state['target_reached'])
        later=persistence_step(state,5.2,False)
        self.assertTrue(later['target_reached']);self.assertEqual(later['first_reached_time_ps'],5.)


class ContinuationTests(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(prefix='B1-CONTINUATION-SYNTHETIC-');self.addCleanup(t.cleanup);self.root=Path(t.name)
        sys.path.insert(0,str(Path(__file__).resolve().parent))
        from titan_path_tpr_pair_contract_smoke import PairContract
        fixture=PairContract();fixture.lock={'ld_seed':123}
        self.text=fixture.text().replace('header:','   init-step = 500000\n   simulation-part = 2\nheader:')
        self.a=self.root/'initial';self.b=self.root/'next';self.a.mkdir();self.b.mkdir()
        (self.a/'tpr.dump').write_text(self.text)

    def test_only_clock_part_steering_and_state_may_change(self):
        changed=self.text.replace('500000','502500').replace('simulation-part = 2','simulation-part = 3').replace('vec[0] = 1','vec[0] = 0.9').replace('init = 0.1','init = -5.2').replace('x[1]={1,0,0}','x[1]={1.1,0,0}')
        (self.b/'tpr.dump').write_text(changed);runtime.continuation_gate(self.a,self.b)

    def test_seed_or_thermostat_drift_rejected(self):
        for old,new in [('ld-seed = 123','ld-seed = 124'),('nsttcouple = 50','nsttcouple = 100')]:
            (self.b/'tpr.dump').write_text(self.text.replace(old,new))
            with self.assertRaises(ValueError):runtime.continuation_gate(self.a,self.b)

    def test_missing_state_or_topology_drift_rejected(self):
        for changed in (self.text.replace('   v[1]={0,0,1}\n',''),self.text.replace('parameter = 1','parameter = 2')):
            (self.b/'tpr.dump').write_text(changed)
            with self.assertRaises(ValueError):runtime.continuation_gate(self.a,self.b)


if __name__=='__main__':unittest.main(verbosity=2)
