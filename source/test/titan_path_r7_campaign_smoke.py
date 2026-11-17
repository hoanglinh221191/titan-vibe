import copy,json,os,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_r7_campaign as c
import titan_path_receding_dispatch as d
from titan_path_multiframe_feedback import new_state
from titan_path_feedback import unseal,seal


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        d.write(self.root/'PROTOCOL.json',{'schema':c.SCHEMA,'route_tracking':'branch-aware-rerouting-v1',
            'rerouting_enabled':True,'max_intervals':400,'interval_ps':5.,'total_new_md_cap_ps':30000.})
        sources=[];self.cells=[]
        for b in (2,3,4):
            name=f'blue__block_{b:02d}__adaptive'
            source={'cell':name,'block':b,'rng_lock':{'source_checkpoint_sha256':str(b)*64},'rng_pair_id':name,
                'source_branch':'b','route':'route','route_points':'points','route_id':'r','reference_sha256':'a'*64}
            sources.append(source)
        for seed in c.SEEDS:
            for source in sources:
                row={**source,'source_cell':source['cell'],'cell':source['cell']+f'__seed_{seed:03d}','planner_seed':seed}
                self.cells.append(row)
                for name in (row['cell'],row['source_cell']):
                    p=self.root/'inputs/cells'/name/'baseline_parameters.tsv';p.parent.mkdir(parents=True,exist_ok=True);p.write_text('same')
        d.write(self.root/'R7_SOURCE_PLAN.json',{'cells':sources})
        d.write(self.root/'PLAN.json',{'cells':self.cells})
    def test_full_design(self):c.validate_plan(self.root,self.cells)
    def test_duplicate_seed_rejected(self):
        rows=copy.deepcopy(self.cells);rows[3]['planner_seed']=11
        with self.assertRaises(ValueError):c.validate_plan(self.root,rows)
    def test_md_rng_change_rejected(self):
        rows=copy.deepcopy(self.cells);rows[0]['rng_lock']={'source_checkpoint_sha256':'f'*64}
        with self.assertRaises(ValueError):c.validate_plan(self.root,rows)
    def test_parameter_drift_rejected(self):
        (self.root/'inputs/cells'/self.cells[0]['cell']/'baseline_parameters.tsv').write_text('drift')
        with self.assertRaises(ValueError):c.validate_plan(self.root,self.cells)
    def test_fifteenth_case_is_dispatched(self):
        d.write(self.root/'generated/PREPARE_READY.json',{})
        for row in self.cells[:-1]:d.write(self.root/'outputs'/row['cell']/'TRAJECTORY_RESULT.json',{'status':'MD_BUDGET_EXHAUSTED'})
        with patch.object(d,'submit_job',return_value='99') as submit:
            d.dispatch_after_prepare(self.root)
        self.assertEqual(submit.call_args.args[1:],('md_14',['--array=14']))
    def test_initial_stop_does_not_block_other_seeds(self):
        d.write(self.root/'generated/PREPARE_READY.json',{})
        d.write(self.root/'outputs'/self.cells[0]['cell']/'TRAJECTORY_RESULT.json',{'status':'INITIAL_PLANNER_STOP'})
        with patch.object(d,'submit_job',return_value='99') as submit:d.dispatch_after_prepare(self.root)
        self.assertEqual(submit.call_args.args[1:],('md_01',['--array=1']))
    def test_all_fifteen_terminal_collects(self):
        d.write(self.root/'generated/PREPARE_READY.json',{})
        for row in self.cells:d.write(self.root/'outputs'/row['cell']/'TRAJECTORY_RESULT.json',{'status':'PLANNER_STOP'})
        with patch.object(d,'submit_job',return_value='99') as submit:d.dispatch_after_prepare(self.root)
        self.assertEqual(submit.call_args.args[1],'collect')
    def test_missing_build_never_dispatches_md(self):
        with patch.object(d,'submit_job',return_value='99') as submit:d.dispatch_after_prepare(self.root)
        self.assertEqual(submit.call_args.args[1],'collect')

    def tpr_fixture(self,change):
        from titan_path_tpr_pair_contract_smoke import PairContract
        f=PairContract();f.setUp();self.addCleanup(f.doCleanups)
        a=self.root/'previous';b=self.root/'prepared';a.mkdir();b.mkdir()
        (a/'segment.tpr').write_bytes(b'SYNTHETIC_A');(b/'segment.tpr').write_bytes(b'SYNTHETIC_B')
        (a/'tpr.dump').write_text(f.text());(b/'tpr.dump').write_text(change(f.text()))
        cell={'source_cell':'source','rng_lock':f.lock}
        return cell,a,b
    def test_same_physical_state_with_different_steering_passes(self):
        cell,a,b=self.tpr_fixture(lambda t:t.replace('vec[0] = 1','vec[0] = 0.9').replace('init = 0.1','init = 0.2'))
        with patch('titan_path_benchmark_runtime.continuation_gate',return_value={}),\
             patch('titan_path_restart_policy.for_cell',return_value=None):
            c.initial_tpr_checks(self.root,cell,{'source_prepared':a},b,{'source':a})
        receipt=d.read(b/'SHARED_CHECKPOINT_TPR_PAIR.json')
        self.assertTrue(receipt['physical_state_equality_required'])
    def test_changed_initial_coordinates_rejected(self):
        cell,a,b=self.tpr_fixture(lambda t:t.replace('x[1]={1,0,0}','x[1]={1.1,0,0}'))
        with patch('titan_path_benchmark_runtime.continuation_gate',return_value={}),\
             patch('titan_path_restart_policy.for_cell',return_value=None):
            with self.assertRaisesRegex(ValueError,'PHYSICAL_STATE_MISMATCH'):
                c.initial_tpr_checks(self.root,cell,{'source_prepared':a},b,{'source':a})
    def test_changed_initial_velocity_rejected(self):
        cell,a,b=self.tpr_fixture(lambda t:t.replace('v[1]={0,0,1}','v[1]={0,0,2}'))
        with patch('titan_path_benchmark_runtime.continuation_gate',return_value={}),\
             patch('titan_path_restart_policy.for_cell',return_value=None):
            with self.assertRaisesRegex(ValueError,'PHYSICAL_STATE_MISMATCH'):
                c.initial_tpr_checks(self.root,cell,{'source_prepared':a},b,{'source':a})

    def resume_fixture(self):
        cell=self.cells[0];source=self.root/'inputs/r7_sources'/cell['source_cell'];source.mkdir(parents=True)
        identity={k:'a'*64 for k in ('route_sha256','reference_sha256','frame_definition_sha256','protocol_sha256')}
        identity['route_id']='r'
        before=new_state(identity,'b'*64,10.,[1.]*3,[2.]*3,'Advance','adaptive')
        paths={}
        for key in ('checkpoint','gro','tpr','none','whole','pullx','xtc','tpr_dump','grompp_log'):
            p=self.root/'upstream'/key;p.parent.mkdir(exist_ok=True);p.write_text('SYNTHETIC_'+key)
            paths[key]={'cluster_path':str(p),'sha256':c.sha(p)}
        (source/'snapshot.pdb').write_text('SYNTHETIC_SNAPSHOT')
        d.write(source/'BINDING.json',{'files':paths,'old_direction':[.25,0,0], 'snapshot_sha256':c.sha(source/'snapshot.pdb')})
        d.write(source/'PREVIOUS_STATE.json',before);d.write(source/'HISTORICAL_STATE.json',before)
        d.write(source/'SEGMENT_RESULT.json',{'output_hashes':{'checkpoint':paths['checkpoint']['sha256']},'end_time_ps':15.,'boundary':[15,0,0]})
        (self.root/'route').write_text('SYNTHETIC_ROUTE');d.write(self.root/'points',[[0,0,0],[2,0,0]])
        d.write(self.root/'inputs/FRAME_DEFINITION.json',{'selected_source_ids':[],'reference_xyz_A':[]})
        (self.root/'inputs/cells'/cell['cell']/'baseline_parameters.tsv').write_text('mass_u\n'+'1\n'*2615)
        final=unseal(new_state(identity,paths['checkpoint']['sha256'],15.,[1.5]*3,[1.8]*3,'Advance','adaptive'))
        final.update(frame_segments=[0]*3,frame_progress_intervals_A=[[1.5,1.5]]*3,suggested_mode='Advance',route_tracking='branch-aware-rerouting-v1')
        return cell,source,paths,seal(final)

    def test_resume_reaudits_ready_and_preserves_branch_seed(self):
        cell,source,paths,final=self.resume_fixture()
        with patch.object(c,'read_gro',return_value=[]),patch.object(c,'read_xvg',return_value=[]),\
             patch.object(c,'observe_segment',return_value=(final,{'reasons':[]},[])) as replay:
            context=c.bind_resume(self.root,cell)
        self.assertEqual(replay.call_count,1)
        self.assertEqual(context['state']['segments'],0)
        self.assertEqual(context['state']['frame_segments'],[0]*3)
        self.assertEqual(context['planner_seed'],11)
        self.assertEqual(context['state']['checkpoint_sha256'],paths['checkpoint']['sha256'])
    def test_changed_upstream_checkpoint_blocks_replay(self):
        cell,source,paths,final=self.resume_fixture()
        Path(paths['checkpoint']['cluster_path']).write_text('TAMPERED')
        with patch.object(c,'observe_segment',side_effect=AssertionError('must not run')):
            with self.assertRaisesRegex(ValueError,'R7_UPSTREAM_HASH'):c.bind_resume(self.root,cell)
    def test_reaudit_hold_cannot_be_promoted(self):
        cell,source,paths,final=self.resume_fixture();body=unseal(final);body['status']='HOLD';final=seal(body)
        with patch.object(c,'read_gro',return_value=[]),patch.object(c,'read_xvg',return_value=[]),\
             patch.object(c,'observe_segment',return_value=(final,{'reasons':['COM_MISMATCH']},[])):
            with self.assertRaisesRegex(ValueError,'REAUDIT_HOLD'):c.bind_resume(self.root,cell)

    def test_submit_registers_afterany_prepare_guard(self):
        with patch.object(d,'verify',return_value=self.cells),patch.object(d,'submit_job',side_effect=['100','101']) as submit,\
             patch.object(sys,'argv',['dispatch','submit','--root',str(self.root)]):d.main()
        self.assertEqual(submit.call_args_list[0].args[1:],('build',[]))
        self.assertEqual(submit.call_args_list[1].args[1:],('prepare_advance',['--dependency=afterany:100']))
    def test_killed_prepare_collects_failure(self):
        d.write(self.root/'generated/jobs/prepare_advance.json',{'job_id':'101'})
        with patch.object(d,'verify',return_value=self.cells),patch.object(d,'submit_job',return_value='102') as submit,\
             patch.dict(os.environ,{'SLURM_JOB_ID':'101'}),\
             patch.object(sys,'argv',['dispatch','prepare-advance-stage','--root',str(self.root)]):d.main()
        self.assertEqual(submit.call_args.args[1],'collect')
        self.assertEqual(d.read(self.root/'generated/PREPARE_FAILURE.json')['error'],'BUILD_MISSING_OR_INTERRUPTED')
    def test_md14_registers_guard_before_worker(self):
        d.write(self.root/'generated/jobs/md_14.json',{'job_id':'114'})
        events=[]
        with patch.object(d,'verify',return_value=self.cells),\
             patch.object(d,'submit_job',side_effect=lambda *a:events.append(('guard',a[1],a[2])) or '115'),\
             patch.object(d,'run',side_effect=lambda *a:events.append(('worker',a[1])) or 0),\
             patch.dict(os.environ,{'SLURM_ARRAY_TASK_ID':'14','SLURM_ARRAY_JOB_ID':'114'}),\
             patch.object(sys,'argv',['dispatch','md-stage','--root',str(self.root)]):self.assertEqual(d.main(),0)
        self.assertEqual(events[0],('guard','advance_14',['--dependency=afterany:114','--export=ALL,B1_R5_ADVANCE_INDEX=14']))
        self.assertEqual(events[1],('worker',14))


if __name__=='__main__':unittest.main()
