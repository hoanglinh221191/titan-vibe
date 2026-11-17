"""R4 synthetic gate, exact index set, fail-closed dispatch and real collector."""
import copy
import json
import os
from pathlib import Path
import unittest
from unittest.mock import patch
import titan_path_benchmark_campaign as campaign
import titan_path_benchmark_continuation as c
import titan_path_benchmark_quota as q
import titan_path_benchmark_runtime as engine
import titan_path_benchmark_campaign_smoke as helpers


class ContinuationTests(unittest.TestCase):
    fake=helpers.CampaignTests.fake
    submissions=helpers.CampaignTests.submissions
    env=helpers.CampaignTests.env

    def setUp(self):
        helpers.CampaignTests.setUp(self)
        q.publish(self.root/'CONTINUATION_POLICY.json',c.POLICY)
        q.publish(self.root/'CAMPAIGN_POLICY.json',c.CAMPAIGN,replace=True)
        protocol=q.read(self.root/'PROTOCOL.json');protocol['campaign_policy_sha256']=q.digest(self.root/'CAMPAIGN_POLICY.json')
        q.publish(self.root/'PROTOCOL.json',protocol,replace=True)
        self.prior=[]
        for i,cell in enumerate(self.cells):
            status='METHOD_STOP' if i<17 else 'TECHNICAL_FAILURE' if i==17 else campaign.NOT_RUN
            self.prior.append({'trajectory_id':cell['cell'],'status':status,'md_executed':False,
                               'completed_segments':0,'scientific_pass':False})
        q.publish(self.root/'inputs/upstream/R3_COLLECTION.json',{'trajectories':self.prior})
        q.publish(self.root/'inputs/upstream/R3_SUBMISSION.json',{'stages':{}})
        q.publish(self.root/'CONTINUATION_BINDING.json',{'r3_bundle_sha256':c.R3_BUNDLE_SHA,
            'authorized_cells':[self.cells[i]['cell'] for i in c.INDICES],
            'files':{p.relative_to(self.root).as_posix():q.digest(p) for p in (self.root/'inputs/upstream').iterdir()}})
        self.assertEqual(campaign.verify(self.root),self.cells)

    def complete(self,key,status='METHOD_STOP',rc=0):
        if key!='build':return helpers.CampaignTests.complete(self,key,status,rc)
        place=q.directory(self.root,key);job=q.read(place/'ACCEPTED.json')['job_id']
        q.publish(self.root/'generated/BUILD_IDENTITY.json',{'synthetic':True})
        q.publish(self.root/'generated/replay_gate/proof.json',{'synthetic':True})
        gate={'status':'PASS_CLUSTER_REPLAY_NO_MD','old_rc':1,'fixed_rc':0,'md_executed':False,
              'artifacts':{'generated/replay_gate/proof.json':q.digest(self.root/'generated/replay_gate/proof.json')},
              'build_identity_sha256':q.digest(self.root/'generated/BUILD_IDENTITY.json'),
              'supervisor_verdict':{'accepted':False,'md_authorized':False,
                  'frames':[{'reason':'STOP_FRAME_RECOVERY_DESCENT'} for _ in range(3)]}}
        q.publish(self.root/'generated/REPLAY_GATE.json',gate)
        q.publish(self.root/'generated/PREPARE_READY.json',{'status':'PASS_PREPARE_NO_MD','md_executed':False,
            'replay_gate_sha256':q.digest(self.root/'generated/REPLAY_GATE.json'),
            'records':[{'cell':self.cells[i]['cell']} for i in c.INDICES]})
        q.publish(place/'DONE.json',{'stage':key,'job_id':job,'returncode':rc,'error':None,
            'scientific_pass':False,'completion_evidence':campaign.completion_evidence(self.root,key)})

    def prefix(self,count):
        q.advance(self.root);self.complete('build')
        for i in c.INDICES[:count]:
            q.advance(self.root);self.complete(f'md_{i:02d}')

    def collect_fixture(self):return helpers.CampaignTests.collect_fixture(self)

    def test_only_fourteen_stage_keys(self):
        self.assertEqual(q.stage_keys(self.root),['build']+[f'md_{i:02d}' for i in range(18,30)]+['collect'])

    def test_all_historical_stage_keys_forbidden(self):
        for i in range(18):
            with self.assertRaisesRegex(ValueError,'INVALID_STAGE'):q.enqueue(self.root,f'md_{i:02d}')
        self.assertFalse(self.submissions())

    def test_build_first_then_md18_only(self):
        self.prefix(0);self.assertIn('SUBMITTED md_18',q.advance(self.root))
        self.assertIn('--array=18-18',self.submissions()[-1])

    def test_exact_twelve_once_then_collect(self):
        self.prefix(12);self.assertIn('SUBMITTED collect',q.advance(self.root))
        arrays=[p for cmd in self.submissions() for p in cmd if p.startswith('--array=')]
        self.assertEqual(arrays,[f'--array={i}-{i}' for i in range(18,30)])

    def test_no_gate_no_md(self):
        self.prefix(0);(self.root/'generated/REPLAY_GATE.json').unlink()
        self.assertIn('SUBMITTED collect',q.advance(self.root))
        self.assertFalse(any('--array=' in p for cmd in self.submissions() for p in cmd))

    def test_gate_artifact_mutation_stops_md(self):
        self.prefix(0);q.publish(self.root/'generated/replay_gate/proof.json',{'changed':True},replace=True)
        self.assertIn('SUBMITTED collect',q.advance(self.root))

    def test_gate_hash_changed_stops_md(self):
        self.prefix(0);q.publish(self.root/'generated/REPLAY_GATE.json',{'fake':True},replace=True)
        self.assertIn('SUBMITTED collect',q.advance(self.root))

    def test_replay_stop_not_permission_to_rewrite_history(self):
        self.prefix(12);q.advance(self.root);report=self.collect_fixture()
        self.assertEqual(report['status'],'PASS_CONTINUATION_12_TECHNICAL_TERMINALS')
        self.assertEqual(report['combined_accounting_status'],'FAILED_WITH_DIAGNOSTICS')
        self.assertEqual(report['accounting']['planned_count'],30)
        self.assertEqual(report['accounting']['technical_terminal_count'],29)
        self.assertEqual(report['trajectories'][17]['status'],'TECHNICAL_FAILURE')
        self.assertEqual(report['continuation']['new_technical_terminals'],12)
        self.assertFalse(report['scientific_pass'])
        for i in range(18):
            row=dict(report['trajectories'][i]);row.pop('evidence_origin')
            self.assertEqual(row,self.prior[i])

    def test_failed_replay_collects_history_plus_twelve_not_run(self):
        q.advance(self.root);self.complete('build',rc=1);q.advance(self.root)
        report=self.collect_fixture()
        self.assertEqual(report['continuation']['new_not_run'],12)
        self.assertEqual(report['status'],'FAILED_CONTINUATION_WITH_DIAGNOSTICS')
        self.assertEqual(report['trajectories'][17]['status'],'TECHNICAL_FAILURE')

    def test_new_failure_stops_next_and_collects(self):
        self.prefix(3);q.advance(self.root);self.complete('md_21','TECHNICAL_FAILURE',1)
        self.assertIn('SUBMITTED collect',q.advance(self.root))
        report=self.collect_fixture()
        self.assertEqual(report['continuation']['new_technical_terminals'],3)
        self.assertEqual(report['continuation']['new_technical_failures'],1)
        self.assertEqual(report['continuation']['new_not_run'],8)

    def test_interrupted_md18_no_retry(self):
        self.prefix(0);q.advance(self.root)
        self.assertIn('SUBMITTED collect',q.advance(self.root))
        self.assertFalse(any('--array=19-19' in cmd for cmd in self.submissions()))

    def test_unknown_acceptance_never_blind_retry(self):
        self.prefix(0);self.behavior='timeout';q.advance(self.root);self.behavior='ok'
        self.assertEqual(q.advance(self.root),'HOLD_ACCEPTANCE_UNKNOWN')

    def test_quota_pause_keeps_gate(self):
        self.prefix(0);self.behavior='quota';self.assertIn('WAIT_QUOTA',q.advance(self.root))
        self.behavior='ok';self.assertIn('SUBMITTED md_18',q.advance(self.root))
        self.assertFalse((self.root/'generated/CAMPAIGN_HALT.json').exists())

    def test_worker_cannot_jump_to_md19(self):
        self.prefix(0)
        with self.assertRaisesRegex(ValueError,'OUT_OF_ORDER'):campaign.worker_guard(self.root,'md_19')

    def test_direct_runtime_forbids_old_index(self):
        with patch.object(engine,'verify',return_value={'cells':self.cells}):
            with self.assertRaisesRegex(ValueError,'HISTORICAL_TRAJECTORY'):engine.run_trajectory(self.root,17)
        self.assertFalse((self.root/'outputs').exists())

    def test_history_mutation_blocks_submit(self):
        p=self.root/'inputs/upstream/R3_COLLECTION.json';q.publish(p,{'trajectories':[]},replace=True)
        with self.assertRaisesRegex(ValueError,'UPSTREAM_BINDING'):q.advance(self.root)
        self.assertFalse(self.submissions())

    def test_policy_cannot_expand_md_set(self):
        q.publish(self.root/'CONTINUATION_POLICY.json',{**c.POLICY,'indices':list(range(30))},replace=True)
        with self.assertRaisesRegex(ValueError,'POLICY_CHANGED'):q.advance(self.root)

    def test_worker_wrong_index_rejected(self):
        self.prefix(0);q.advance(self.root)
        env=self.env('md_18','70002');env['SLURM_ARRAY_TASK_ID']='17'
        with patch.dict(os.environ,env,clear=True):
            with self.assertRaisesRegex(ValueError,'ARRAY_INDEX'):q.authorize(self.root)

    def test_actual_prepare_stops_before_grompp_if_replay_fails(self):
        self.prefix(0)
        with patch.object(engine,'verify',return_value={'cells':self.cells}),\
             patch.object(engine,'read',side_effect=lambda p:{'build_job':'12'} if p.name=='SUBMISSION.json' else q.read(p)),\
             patch.dict(os.environ,{'SLURM_JOB_ID':'12'}),patch.object(engine,'logged'),patch.object(engine.replay,'build'),\
             patch.object(c,'replay_gate',side_effect=ValueError('SYNTHETIC_REPLAY_FAIL')),patch.object(engine,'make_snapshot') as snap:
            with self.assertRaisesRegex(ValueError,'REPLAY_FAIL'):engine.prepare(self.root)
        snap.assert_not_called()

    def test_gate_cannot_admit_unexpected_replay_continue(self):
        self.prefix(0)
        gate=q.read(self.root/'generated/REPLAY_GATE.json');gate['supervisor_verdict']['accepted']=True
        q.publish(self.root/'generated/REPLAY_GATE.json',gate,replace=True)
        ready=q.read(self.root/'generated/PREPARE_READY.json');ready['replay_gate_sha256']=q.digest(self.root/'generated/REPLAY_GATE.json')
        with self.assertRaisesRegex(ValueError,'VERDICT_CHANGED'):c.validate_gate(self.root,ready)

    def test_gate_cannot_prepare_old_trajectories(self):
        self.prefix(0);ready=q.read(self.root/'generated/PREPARE_READY.json')
        ready['records']=[{'cell':r['cell']} for r in self.cells]
        with self.assertRaisesRegex(ValueError,'ONLY_UNTOUCHED12'):c.validate_gate(self.root,ready)

    def test_prior_submission_conflict_blocks_new_md(self):
        path=self.root/'inputs/upstream/R3_SUBMISSION.json'
        q.publish(path,{'stages':{'md_18':{'ACCEPTED':{'job_id':'13'}}}},replace=True)
        binding=q.read(self.root/'CONTINUATION_BINDING.json')
        binding['files']['inputs/upstream/R3_SUBMISSION.json']=q.digest(path)
        q.publish(self.root/'CONTINUATION_BINDING.json',binding,replace=True)
        with self.assertRaisesRegex(ValueError,'HISTORICAL_DISPATCH'):q.advance(self.root)

    def test_actual_prepare_visits_only_twelve(self):
        self.prefix(0)
        # Remove only the synthetic completion receipt to exercise the real writer.
        (self.root/'generated/PREPARE_READY.json').unlink()
        q.publish(self.root/'inputs/FRAME_DEFINITION.json',{'selected_source_ids':[],'reference_xyz_A':[]})
        visited=[]
        def snapshot(root,branch,place):place.mkdir(parents=True)
        def context(root,cell):return {'snapshot':self.root/'synthetic_snapshot.pdb','state':{}}
        def proposal(root,cell,*args,**kwargs):visited.append(cell['cell']);return {'accepted':False}
        real_sha=engine.sha
        with patch.object(engine,'verify',return_value={'cells':self.cells}),\
             patch.object(engine,'read',side_effect=lambda p:{'build_job':'12'} if p.name=='SUBMISSION.json' else q.read(p)),\
             patch.dict(os.environ,{'SLURM_JOB_ID':'12','GMX_MPI':'synthetic_gmx'}),\
             patch.object(engine,'sha',side_effect=lambda p:engine.GMX_SHA if str(p)=='synthetic_gmx' else real_sha(p)),\
             patch.object(engine,'logged',return_value='GROMACS version: 2026.2\n'),patch.object(engine.replay,'build'),\
             patch.object(c,'replay_gate'),patch.object(engine,'make_snapshot',side_effect=snapshot),\
             patch.object(engine,'initial_context',side_effect=context),patch.object(engine.paired,'native_proposal',side_effect=proposal),\
             patch.object(engine,'frame_geometry',return_value=([],[],[],[])),patch.object(engine,'read_gro',return_value=[{'box':[]}]),\
             patch.object(engine,'define_target',return_value={}),patch.object(engine,'measure_target',return_value={'condition_met':False}):
            engine.prepare(self.root)
        self.assertEqual(visited,[self.cells[i]['cell'] for i in c.INDICES])
        self.assertEqual(len(q.read(self.root/'generated/PREPARE_READY.json')['records']),12)


if __name__=='__main__':unittest.main()
