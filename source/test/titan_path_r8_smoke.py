import json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_exploratory import navigation_warnings
import titan_path_receding_dispatch as dispatch
import titan_path_receding_runtime as runtime
import titan_path_exploratory as explore
from titan_path_feedback import seal,unseal
from titan_path_multiframe_feedback import new_state


class R8Tests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup);self.root=Path(tmp.name)
    def test_large_retreat_is_warning(self):
        r=navigation_warnings({'progress_A':[10]*3,'distance_A':[2]*3},{'progress_A':[4]*3,'distance_A':[9]*3})
        self.assertFalse(r['navigation_stop']);self.assertTrue(r['retreat']);self.assertTrue(r['outside_4A'])
    def test_independent_failure_dispatches_next(self):
        dispatch.write(self.root/'PROTOCOL.json',{'independent_cell_failures':True})
        dispatch.write(self.root/'PLAN.json',{'cells':[{'cell':'a'},{'cell':'b'}]})
        dispatch.write(self.root/'generated/PREPARE_READY.json',{})
        dispatch.write(self.root/'outputs/a/TRAJECTORY_RESULT.json',{'status':'TECHNICAL_FAILURE'})
        with patch.object(dispatch,'submit_job',return_value='1') as s:dispatch.dispatch_after_prepare(self.root)
        self.assertEqual(s.call_args.args[1:],('md_01',['--array=1']))
    def test_legacy_failure_still_collects(self):
        dispatch.write(self.root/'PROTOCOL.json',{})
        dispatch.write(self.root/'PLAN.json',{'cells':[{'cell':'a'},{'cell':'b'}]})
        dispatch.write(self.root/'generated/PREPARE_READY.json',{})
        dispatch.write(self.root/'outputs/a/TRAJECTORY_RESULT.json',{'status':'TECHNICAL_FAILURE'})
        with patch.object(dispatch,'submit_job',return_value='1') as s:dispatch.dispatch_after_prepare(self.root)
        self.assertEqual(s.call_args.args[1],'collect')
    def test_initial_failure_does_not_abort_other_cell(self):
        import titan_path_r8_campaign as campaign
        cells=[{'cell':'a','prior_intervals':2},{'cell':'b','prior_intervals':0}]
        with patch.object(campaign,'bind_resume',side_effect=[ValueError('BAD_SOURCE'),{}]),\
             patch.object(runtime,'propose',return_value=False):
            rows=campaign.build_initial_cells(self.root,cells)
        self.assertEqual(len(rows),2)
        self.assertEqual(dispatch.read(self.root/'outputs/a/TRAJECTORY_RESULT.json')['status'],'PREPARE_CELL_FAILURE')
        self.assertEqual(dispatch.read(self.root/'outputs/b/TRAJECTORY_RESULT.json')['status'],'NO_SAFE_INITIAL_DIRECTION')
    def test_terminal_old_episode_does_not_stop_new_policy(self):
        dispatch.write(self.root/'PROTOCOL.json',{'rerouting_enabled':True,'route_tracking':'branch-aware-rerouting-v1',
                       'navigation_guard_policy':'warning-only-v1'})
        identity={k:'a'*64 for k in ('route_sha256','reference_sha256','frame_definition_sha256','protocol_sha256')};identity['route_id']='x'
        context={'state':new_state(identity,'b'*64,0,[1]*3,[1]*3,'Advance','adaptive'),
                 'reroute_episode':{'status':'STOP_EXCURSION_DISTANCE'}}
        def stop(root,cell,context,work,**kw):
            work.mkdir(parents=True);r={'status':'NO_FEASIBLE_INITIAL_STRING'};dispatch.write(work/'RESULT.json',r);return r
        def safe(root,cell,context,work,**kw):
            work.mkdir(parents=True);r={'status':'FIRST_INTERVAL_PROPOSED','kind':'EXPLORATORY_REROUTE','delta_A':[.25,0,0],'immediate_clearance_A':2}
            dispatch.write(work/'RESULT.json',r);return r
        with patch.object(runtime,'plan_string',side_effect=stop),patch.object(explore,'plan',side_effect=safe),\
             patch.object(runtime,'native_query',return_value=([],[{}])),\
             patch('titan_path_reroute_planner.native_sweep_check',return_value={'passes':True}):
            self.assertTrue(runtime.propose(self.root,{},context,self.root/'selection'))


if __name__=='__main__':unittest.main()
