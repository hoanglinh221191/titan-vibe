"""Positive and adversarial RNG/full-TPR contract tests, no native processes."""
import copy
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_tpr_pair_contract as c
import titan_path_paired_loop as p


class PairContract(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(prefix='TPR-PAIR-SYNTHETIC-');self.addCleanup(t.cleanup);self.root=Path(t.name)
        self.a=self.root/'a.tpr';self.b=self.root/'b.tpr';self.a.write_bytes(b'SYNTHETIC');self.b.write_bytes(b'SYNTHETIC')
        self.lock=c.seed_lock('pair0','a'*64)

    def text(self):
        return f'''file.tpr:
inputrec:
   integrator = md
   continuation = true
   ld-seed = {self.lock['ld_seed']}
   tcoupl = V-rescale
   pcoupl = C-rescale
   pcoupltype = Isotropic
   nsttcouple = 50
   nstpcouple = 25
   tau-p = 2
   nstenergy = 100
   ref-t: 300 300
   tau-t: 0.5 0.5
   pull-coord 0:
      vec[0] = 1
      vec[1] = 0
      vec[2] = 0
      init = 0.1
      rate = 0.005
      k = 600
header:
   natoms = 2
topology:
   parameter = 1.00000
box (3x3):
   box[0]={{1,0,0}}
   box[1]={{0,1,0}}
   box[2]={{0,0,1}}
x (2x3):
   x[0]={{0,0,0}}
   x[1]={{1,0,0}}
v (2x3):
   v[0]={{0,1,0}}
   v[1]={{0,0,1}}
Group statistics
T-Coupling : 2 (total 2 atoms)
'''

    def test_seed_is_deterministic_positive_and_changes_with_pair_identity(self):
        self.assertEqual(self.lock,c.seed_lock('pair0','a'*64))
        self.assertTrue(1<=self.lock['ld_seed']<=2147483647)
        self.assertNotEqual(self.lock['ld_seed'],c.seed_lock('pair1','a'*64)['ld_seed'])
        self.assertNotEqual(self.lock['ld_seed'],c.seed_lock('pair0','b'*64)['ld_seed'])

    def test_invalid_or_tampered_lock(self):
        for value in (-1,0,True,2**40,None):
            lock=copy.deepcopy(self.lock);lock['ld_seed']=value
            with self.subTest(value=value),self.assertRaises(c.NativeContractError):c.validate_seed_lock(lock)

    def test_complete_identical_first_pair(self):
        r=c.validate_pair(self.text(),self.text(),self.a,self.b,self.lock,initial_same_state=True)
        self.assertEqual(r['scope'],'initial_binary_identical');self.assertFalse(r['md_executed'])

    def test_same_checkpoint_but_different_rng_rejected(self):
        bad=self.text().replace(str(self.lock['ld_seed']),str(self.lock['ld_seed']+1))
        with self.assertRaisesRegex(c.NativeContractError,'STOCHASTIC'):c.validate_pair(self.text(),bad,self.a,self.b,self.lock,initial_same_state=True)

    def test_duplicate_or_missing_seed_rejected(self):
        for bad in (self.text().replace('   ld-seed','   missing-seed'),self.text().replace('   tcoupl','   ld-seed = 123\n   tcoupl')):
            with self.assertRaises(c.NativeContractError):c.validate_stochastic(bad,self.lock)

    def test_thermostat_pressure_and_coupling_parameters(self):
        for old,new in [('V-rescale','Nose-Hoover'),('C-rescale','Berendsen'),('nsttcouple = 50','nsttcouple = 100'),
                        ('nstpcouple = 25','nstpcouple = 50'),('tau-p = 2','tau-p = 3'),('300 300','300 310'),('0.5 0.5','0.4 0.5')]:
            with self.subTest(field=old),self.assertRaises(c.NativeContractError):c.validate_stochastic(self.text().replace(old,new),self.lock)

    def test_dump_equal_but_binary_different_rejected(self):
        self.b.write_bytes(b'SAME_PRINTED_FLOAT_DIFFERENT_BITS')
        with self.assertRaisesRegex(c.NativeContractError,'BINARY'):c.validate_pair(self.text(),self.text(),self.a,self.b,self.lock,initial_same_state=True)

    def test_topology_difference_rejected_in_later_pair(self):
        with self.assertRaisesRegex(c.NativeContractError,'TOPOLOGY'):c.validate_pair(self.text(),self.text().replace('parameter = 1','parameter = 2'),self.a,self.b,self.lock,initial_same_state=False)

    def test_nonsteering_field_and_unknown_field_not_ignored(self):
        for bad in (self.text().replace('nstenergy = 100','nstenergy = 50'),self.text().replace('header:','new-option = 1\nheader:')):
            with self.assertRaisesRegex(c.NativeContractError,'INPUTREC'):c.validate_pair(self.text(),bad,self.a,self.b,self.lock,initial_same_state=False)

    def test_authorized_steering_and_evolving_state_only(self):
        changed=self.text().replace('vec[0] = 1','vec[0] = 0.9').replace('init = 0.1','init = 0.2').replace('x[1]={1,0,0}','x[1]={1.1,0,0}')
        r=c.validate_pair(self.text(),changed,self.a,self.b,self.lock,initial_same_state=False)
        self.assertFalse(r['state_equality_required']);self.assertTrue(r['later_text_precision_not_binary_equivalence'])
        with self.assertRaises(c.NativeContractError):c.validate_pair(self.text(),changed,self.a,self.b,self.lock,initial_same_state=True)

    def test_rate_or_spring_change_not_authorized(self):
        for old,new in [('rate = 0.005','rate = 0.006'),('k = 600','k = 601')]:
            with self.assertRaises(c.NativeContractError):c.validate_pair(self.text(),self.text().replace(old,new),self.a,self.b,self.lock,initial_same_state=False)

    def test_incomplete_duplicate_nonfinite_state_rejected(self):
        for bad in (self.text().replace('   v[1]={0,0,1}\n',''),self.text().replace('x[1]','x[0]'),
                    self.text().replace('v[1]={0,0,1}','v[1]={0,nan,1}'),self.text().replace('Group statistics','No summary'),
                    self.text().replace('x (2x3)','x (3x3)')):
            with self.assertRaises(c.NativeContractError):c.split_full_dump(bad)

    def test_old_receipt_blocks_actual_dispatch(self):
        (self.root/'generated').mkdir();(self.root/'generated/SUBMISSION.json').write_text('{"array_job":"123"}')
        # Support both the historical single receipt and the R5/R6 per-task
        # authorization adapter without bypassing either production check.
        (self.root/'generated/jobs').mkdir()
        (self.root/'generated/jobs/md_00.json').write_text('{"job_id":"123"}')
        (self.root/'PLAN.json').write_text('{"cells":[{"cell":"synthetic"}]}')
        with patch.dict(os.environ,{'SLURM_JOB_PARTITION':'batch','SLURM_ARRAY_JOB_ID':'123','SLURM_ARRAY_TASK_ID':'0'},clear=True),\
             patch.object(p,'logged',side_effect=AssertionError('must not execute')) as command:
            with self.assertRaisesRegex(c.NativeContractError,'RNG_FULL_TPR_RECEIPT_REQUIRED'):p.run_one(self.root,{}, {},self.root,{})
        self.assertEqual(command.call_count,0)

    def test_box_invalid_missing_or_nonfinite_rejected(self):
        for bad in (self.text().replace('box[0]={1,0,0}','box[0]={0,0,0}'),self.text().replace('box[0]={1,0,0}','box[0]={nan,0,0}'),self.text().replace('   box[2]={0,0,1}\n','')):
            with self.assertRaises(c.NativeContractError):c.split_full_dump(bad)

    def test_group_statistics_not_masked_as_evolving_state(self):
        changed=self.text().replace('T-Coupling : 2','T-Coupling : 1')
        with self.assertRaisesRegex(c.NativeContractError,'GROUP_STATISTICS'):c.validate_pair(self.text(),changed,self.a,self.b,self.lock,initial_same_state=False)

    def test_receipt_hash_drift(self):
        segment=self.root/'segment';segment.mkdir();(segment/'segment.tpr').write_bytes(b'A');(segment/'tpr.dump').write_text(self.text())
        receipt={'rng_lock':self.lock,'tpr_sha256':c.sha(segment/'segment.tpr'),'native_dump_sha256':c.sha(segment/'tpr.dump')}
        c.validate_receipt(segment,receipt)
        (segment/'segment.tpr').write_bytes(b'B')
        with self.assertRaisesRegex(c.NativeContractError,'CHANGED'):c.validate_receipt(segment,receipt)


if __name__=='__main__':unittest.main(verbosity=2)
