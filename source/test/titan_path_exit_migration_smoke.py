"""Synthetic unit tests; not a replacement for native checkpoint/TPR gates."""
import copy,json,os,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from titan_path_exit_migration import nearest_vector,migration_header,handoff,policy
from titan_path_multiframe_replay import sha

HEADER='''inputrec:
 init-step = 100
 simulation-part = 10
 dt = 0.002
 ld-seed = 17
 pull-group 0:
  nat = 0
 pull-group 1:
  nat = 2603
  pbcatom = 1604
  ind = 0 1 2
 pull-group 2:
  nat = 12
  pbcatom = 2608
  ind = 2603 2604
 pull-coord 0:
  vec[0] = 1
  vec[1] = 0
  vec[2] = 0
  init = -1
  rate = 0.005
  k = 600
'''

class Migration(unittest.TestCase):
    def test_nearest_periodic_image(self):
        np.testing.assert_allclose(nearest_vector([9.9,.2,0],np.eye(3)*10),[-.1,.2,0],atol=1e-12)
    def test_ambiguous_image_rejected(self):
        with self.assertRaisesRegex(ValueError,'AMBIGUOUS'):nearest_vector([5,0,0],np.eye(3)*10)
    def test_nonfinite_rejected(self):
        with self.assertRaisesRegex(ValueError,'GEOMETRY'):nearest_vector([float('nan'),0,0],np.eye(3))
    def test_singular_rejected(self):
        with self.assertRaisesRegex(ValueError,'GEOMETRY'):nearest_vector([0,0,0],np.zeros((3,3)))
    def test_triclinic_lattice_invariance(self):
        box=np.array([[6,0,0],[1,5,0],[.3,.4,5.]])
        v=np.array([.2,.3,.4]);np.testing.assert_allclose(nearest_vector(v+np.array([3,-4,2])@box,box),v,atol=1e-12)
    def test_only_reference_and_steering_may_change(self):
        changed=HEADER.replace('2603\n','16\n').replace('1604','1235').replace('0 1 2','1183 1185').replace('100','200').replace('= 10','= 11').replace('init = -1','init = -2')
        self.assertEqual(migration_header(HEADER),migration_header(changed))
    def test_ligand_coupling_seed_and_rate_cannot_change(self):
        for old,new in [('nat = 12','nat = 11'),('pbcatom = 2608','pbcatom = 2609'),('2603 2604','2603 2605'),
                        ('dt = 0.002','dt = 0.003'),('ld-seed = 17','ld-seed = 19'),('rate = 0.005','rate = 0.006'),('k = 600','k = 601')]:
            with self.subTest(old=old):self.assertNotEqual(migration_header(HEADER),migration_header(HEADER.replace(old,new)))
    def test_missing_group_rejected(self):
        with self.assertRaisesRegex(ValueError,'GROUP_BLOCK'):migration_header(HEADER.replace('pull-group 1:','pull-group 4:'))
    def test_force_extension_preserved_across_reference_change(self):
        # Old actual-reference = -0.2 nm; new group/direction have different projection.
        ctx={'boundary':[10,1,1.2,1,0,0,0,0,0,1,0,0],'old_direction':[1,0,0],'delta':[0,.25,0]}
        with patch('titan_path_exit_migration.checkpoint_geometry',return_value=(np.array([.2,.7,0]),{})):
            h,proof=handoff(None,None,ctx,'SYNTHETIC',{})
        self.assertAlmostEqual(h.scalar_force_before_kj_mol_nm,h.scalar_force_after_kj_mol_nm)
        self.assertAlmostEqual(h.signed_extension_nm,-.2)
        self.assertTrue(proof['reaction_force_distribution_changed'])
        self.assertFalse(proof['vector_force_continuity_claimed'])
    def test_invalid_old_boundary_rejected(self):
        ctx={'boundary':[10,1,1.2,2,0,0,0,0,0,2,0,0],'old_direction':[1,0,0],'delta':[0,.25,0]}
        with self.assertRaises(ValueError):handoff(None,None,ctx,'SYNTHETIC',{})
    def test_policy_requires_exact_bound_selection(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'inputs').mkdir()
            selection={'schema':'titan-path-exit-reference-v1','group_name':'TITAN_PATH_EXIT_REF',
                'atom_ids_0based':list(range(16)),'pbcatom_0based':1,'residues':[{'residue':['A',n,'X']} for n in (76,77,79,80)]}
            p=root/'inputs/EXIT_REFERENCE.json';p.write_text(json.dumps(selection))
            (root/'PROTOCOL.json').write_text(json.dumps({'exit_reference_sha256':sha(p)}))
            legacy={'policy_sha256':'a'*64};self.assertEqual(policy(root,legacy)['reference_atom_ids'],list(range(16)))
            p.write_text(p.read_text()+' ')
            with self.assertRaisesRegex(ValueError,'BINDING'):policy(root,legacy)
    def test_legacy_policy_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'PROTOCOL.json').write_text('{}');legacy={'policy_sha256':'legacy'}
            self.assertIs(policy(root,legacy),legacy)
    def test_recorded_r9_cell_job_required_before_md(self):
        import titan_path_paired_loop as paired
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'generated/jobs').mkdir(parents=True)
            (root/'PROTOCOL.json').write_text(json.dumps({'schema':'r9-exit-reference-continuation-v1'}))
            (root/'PLAN.json').write_text(json.dumps({'cells':[{'cell':'one'}]}))
            (root/'generated/jobs/md_00.json').write_text(json.dumps({'job_id':'123'}))
            with patch.dict(os.environ,{'SLURM_JOB_PARTITION':'batch','SLURM_ARRAY_TASK_ID':'0','SLURM_ARRAY_JOB_ID':'124'}),\
                 patch.object(paired,'logged',side_effect=AssertionError('must not execute MD')):
                with self.assertRaisesRegex(ValueError,'UNRECORDED_MD_JOB'):paired.run_one(root,{'cell':'one'},{},root,{})
            with patch.dict(os.environ,{'SLURM_JOB_PARTITION':'batch','SLURM_ARRAY_TASK_ID':'0','SLURM_ARRAY_JOB_ID':'123'}),\
                 patch.object(paired,'validate_receipt',side_effect=ValueError('REACHED_NATIVE_GATE')),\
                 patch.object(paired,'logged',side_effect=AssertionError('must not execute MD')):
                with self.assertRaisesRegex(ValueError,'MD_CELL_IDENTITY'):paired.run_one(root,{'cell':'wrong'},{},root,{})
                with self.assertRaisesRegex(ValueError,'REACHED_NATIVE_GATE'):paired.run_one(root,{'cell':'one'},{},root,{})

if __name__=='__main__':unittest.main()
