import json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_receding_runtime as runtime
import titan_path_reroute_planner as fallback
from titan_path_multiframe_feedback import new_state
from titan_path_recovery import start_excursion


class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.sweep=patch.object(fallback,'native_sweep_check',return_value={'passes':True,'complete':True,'clearance_A':2.})
        self.sweep.start();self.addCleanup(self.sweep.stop)
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        (self.root/'PROTOCOL.json').write_text(json.dumps({'rerouting_enabled':True,'direct_recovery_enabled':True,'route_tracking':'branch-aware-rerouting-v1'}))
        identity={'route_id':'SYNTHETIC',**{k:'a'*64 for k in ('route_sha256','reference_sha256','frame_definition_sha256','protocol_sha256')}}
        state=new_state(identity,'b'*64,0.,[1]*3,[2]*3,'Advance','adaptive')
        self.context={'state':state,'planner_seed':11}
    def tearDown(self):self.tmp.cleanup()
    def primary(self,root,cell,context,work,**kwargs):
        work.mkdir(parents=True);r={'status':'NO_FEASIBLE_INITIAL_STRING'}
        (work/'RESULT.json').write_text(json.dumps(r));return r
    def fallback(self,root,cell,context,work,**kwargs):
        work.mkdir(parents=True);e=start_excursion([2]*3,[1]*3)
        e.update(entry_state_sha256=context['state']['record_sha256'],entry_checkpoint_sha256='b'*64,goal_progress_A=[1.5]*3)
        r={'status':'FIRST_INTERVAL_PROPOSED','kind':'BOUNDED_REROUTE','episode':e,'delta_A':[.25,0,0],'immediate_clearance_A':2.}
        (work/'RESULT.json').write_text(json.dumps(r));return r
    def test_fallback_is_explicit_and_episode_persisted(self):
        with patch.object(runtime,'plan_string',side_effect=self.primary),patch.object(fallback,'plan',side_effect=self.fallback),\
             patch.object(runtime,'native_query',return_value=([],[{'total':-1.,'energy_key':1.}])):
            accepted=runtime.propose(self.root,{'cell':'synthetic'},self.context,self.root/'selection')
        self.assertTrue(accepted);self.assertIn('reroute_episode',self.context)
        r=json.loads((self.root/'selection/RESULT.json').read_text());self.assertEqual(r['kind'],'BOUNDED_REROUTE')
        admission=json.loads((self.root/'selection/R5_ADMISSION.json').read_text());self.assertEqual(admission['planner_kind'],'BOUNDED_REROUTE')
    def test_legacy_protocol_does_not_invoke_fallback(self):
        (self.root/'PROTOCOL.json').write_text('{}')
        with patch.object(runtime,'plan_string',side_effect=self.primary),patch.object(fallback,'plan',side_effect=AssertionError('must not run')):
            self.assertFalse(runtime.propose(self.root,{},self.context,self.root/'selection'))
    def test_fallback_requires_tracking_policy(self):
        (self.root/'PROTOCOL.json').write_text('{"rerouting_enabled":true}')
        with self.assertRaises(ValueError):runtime.propose(self.root,{},self.context,self.root/'selection')

    def accepted_primary(self,root,cell,context,work,**kwargs):
        work.mkdir(parents=True)
        r={'status':'FIRST_INTERVAL_PROPOSED','delta_A':[.25,0,0],'immediate_clearance_A':2.}
        (work/'RESULT.json').write_text(json.dumps(r));return r

    def test_primary_success_does_not_reset_active_excursion(self):
        episode=start_excursion([2]*3,[1]*3);episode['intervals']=7
        self.context['reroute_episode']=episode
        with patch.object(runtime,'plan_string',side_effect=self.accepted_primary),\
             patch.object(runtime,'native_query',return_value=([],[{'total':-1.,'energy_key':1.}])):
            self.assertTrue(runtime.propose(self.root,{},self.context,self.root/'selection'))
        self.assertEqual(self.context['reroute_episode']['intervals'],7)
        self.assertFalse((self.root/'selection/REROUTE_EXIT.json').exists())

    def test_observed_return_closes_episode(self):
        episode=start_excursion([2]*3,[1]*3);episode['status']='RETURN_OBSERVED'
        self.context['reroute_episode']=episode
        with patch.object(runtime,'plan_string',side_effect=self.accepted_primary),\
             patch.object(runtime,'native_query',return_value=([],[{'total':-1.,'energy_key':1.}])):
            self.assertTrue(runtime.propose(self.root,{},self.context,self.root/'selection'))
        self.assertNotIn('reroute_episode',self.context)
        self.assertEqual(json.loads((self.root/'selection/REROUTE_EXIT.json').read_text())['reason'],'ACTUAL_RETURN_OBSERVED')

    def test_exhausted_episode_cannot_be_bypassed_by_primary(self):
        episode=start_excursion([2]*3,[1]*3);episode['status']='STOP_RETURN_DEADLINE'
        self.context['reroute_episode']=episode
        with patch.object(runtime,'plan_string',side_effect=AssertionError('must not run')):
            self.assertFalse(runtime.propose(self.root,{},self.context,self.root/'selection'))


if __name__=='__main__':unittest.main()
