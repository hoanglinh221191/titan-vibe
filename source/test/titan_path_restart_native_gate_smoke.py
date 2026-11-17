"""Integration gates using newly produced no-MD native TPRs."""
import copy
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import titan_path_restart_policy as policy
import titan_path_benchmark_runtime as runtime
import titan_path_paired_loop as paired

REPO = Path(__file__).resolve().parents[1]
NATIVE = REPO/'.codex_tmp/b1_restart_fix_native_20260914_r1'
RAW = REPO/'.codex_tmp/b1_r2_native_restart_20260914'


class NativeIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); (self.root/'inputs').mkdir()
        for relative in ('PROTOCOL.json', policy.POLICY_PATH):
            shutil.copyfile(NATIVE/relative, self.root/relative)
        plan = json.loads((NATIVE/'PLAN.json').read_text())
        self.cell = next(c for c in plan['cells'] if c['cell']=='blue__block_03__tangent')
        self.initial = self.root/'initial'; self.current = self.root/'current'
        shutil.copytree(NATIVE/'generated/initial'/self.cell['cell']/'prepared', self.initial)
        shutil.copytree(NATIVE/'generated/retained_snapshot_producers'/self.cell['cell']/'prepared', self.current)
        self.contract = policy.for_cell(self.root, self.cell)

    def test_real_derived_radius_change_passes_only_opt_in(self):
        with self.assertRaisesRegex(ValueError, 'NONSTEERING_DRIFT'):
            runtime.continuation_gate(self.initial, self.current)
        result = runtime.continuation_gate(self.initial, self.current, policy=self.contract)
        self.assertEqual(result['initial']['rlist_nm'], 1.223)
        self.assertEqual(result['current']['rlist_nm'], 1.222)

    def test_seed_change_still_rejected_with_auto_radius_policy(self):
        path = self.current/'tpr.dump'; text = path.read_text(); value = self.cell['rng_lock']['ld_seed']
        old = 'ld-seed                        = '+str(value)
        self.assertIn(old, text); path.write_text(text.replace(old, old+'0', 1))
        with self.assertRaisesRegex(ValueError, 'NONSTEERING_DRIFT'):
            runtime.continuation_gate(self.initial, self.current, policy=self.contract)

    def test_pre_md_receipts_survive_platform_change(self):
        receipt = json.loads((self.initial/'PREPARED.json').read_text())
        policy.validate_pre_md(self.root, self.cell, self.initial, receipt)

    def test_removed_neighbor_proof_blocks_real_run_one_before_mdrun(self):
        receipt = json.loads((self.initial/'PREPARED.json').read_text()); receipt.pop('auto_verlet_proof')
        (self.root/'generated').mkdir(); (self.root/'generated/SUBMISSION.json').write_text('{"array_job":"12"}')
        with patch.dict(os.environ, {'SLURM_JOB_PARTITION':'batch','SLURM_ARRAY_JOB_ID':'12'}), patch.object(paired,'logged') as run:
            with self.assertRaisesRegex(ValueError, 'RESTART_RECEIPT_CHANGED'):
                paired.run_one(self.root, self.cell, {}, self.initial, receipt)
        run.assert_not_called()

    def test_changed_producer_log_invalidates_pre_md_receipt(self):
        receipt = json.loads((self.initial/'PREPARED.json').read_text())
        path = self.initial/'tpr.grompp.log'; path.write_text(path.read_text()+'\nchanged\n')
        with self.assertRaisesRegex(ValueError, 'RESTART_RECEIPT_CHANGED'):
            policy.validate_pre_md(self.root, self.cell, self.initial, receipt)

    def test_forged_margin_receipt_rejected(self):
        receipt = json.loads((self.initial/'PREPARED.json').read_text())
        receipt['pbc_reference_guard']['groups'][0]['remaining_margin_nm'] += .01
        with self.assertRaisesRegex(ValueError, 'GEOMETRY_RECEIPT_VALUE'):
            policy.validate_pre_md(self.root, self.cell, self.initial, receipt)

    def test_actual_pbcatom_change_rejected(self):
        receipt = json.loads((self.initial/'PREPARED.json').read_text())
        path = self.initial/'tpr.dump'; text = path.read_text()
        import re
        changed, count = re.subn(r'(pbcatom\s*=\s*)1646', r'\g<1>1301', text, count=1)
        self.assertEqual(count, 1); path.write_text(changed)
        with self.assertRaisesRegex(ValueError, 'ANCHOR_CHANGED'):
            policy.geometry_guard(changed, self.contract)

    def test_legacy_source_cannot_be_silently_continued_with_new_anchors(self):
        cell = copy.deepcopy(self.cell); cell['source_branch']='branch_a'
        cell['rng_lock']['source_checkpoint_sha256']=policy.b1_policy()['branches']['branch_a']['equilibrium_checkpoint_sha256']
        current = RAW/'outputs/trajectories/orange__block_00__adaptive/segment_002/prepared'
        prior = RAW/'outputs/trajectories/orange__block_00__adaptive/segment_001/prepared'
        checkpoint = prior/'md/invocations/invocation_0000/state.cpt'
        context={'state':{},'checkpoint':checkpoint,'last_tpr':prior/'segment.tpr','snapshot':prior/'observed/endpoint.pdb'}
        boundary={'time_ps':1010.,'checkpoint_sha256':policy.sha(checkpoint)}
        with patch.object(paired,'unseal',return_value=boundary), patch.object(paired,'expected_pull_membership',return_value={}),\
             patch.object(paired,'logged',side_effect=[(current/'input_checkpoint.dump').read_text(),(prior/'tpr.dump').read_text()]) as calls:
            with self.assertRaisesRegex(ValueError, 'ANCHOR_CHANGED'):
                paired.prepare_segment(self.root,cell,context,self.root/'attempt')
        self.assertEqual(calls.call_count,2)
        self.assertTrue(all('dump' in c.args[0] for c in calls.call_args_list))
        self.assertFalse((self.root/'attempt/PREPARED.json').exists())


if __name__ == '__main__':
    unittest.main()
