import hashlib,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from titan_path_package_files import static_files
from titan_path_source_reference import validate_source_reference
from titan_path_rerouting import TrackedRoute,TransitionTrackedRoute
from titan_path_exit_reference import group_com,append_index_group


class Repairs(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(prefix='path-repairs-',dir=ROOT/'.codex_tmp');self.addCleanup(t.cleanup);self.root=Path(t.name)
    def test_inventory_never_descends_runtime(self):
        (self.root/'src').mkdir();(self.root/'src/code.py').write_text('source')
        (self.root/'outputs').mkdir();(self.root/'outputs/data.bin').write_bytes(b'ignore')
        (self.root/'FILES.sha256').write_text('lock')
        original=Path.iterdir
        def guarded(p):
            if p.name in ('outputs','generated','downloads','tryton_logs'):raise AssertionError('runtime traversal')
            return original(p)
        with patch.object(Path,'iterdir',guarded),patch.object(Path,'rglob',side_effect=AssertionError('recursive glob')):
            self.assertEqual(static_files(self.root),{'src/code.py'})
    def reference(self,init=-10.727245447935065):
        path=self.root/'segment.tpr';path.write_bytes(b'SYNTHETIC_TPR')
        (self.root/'PREPARED.json').write_text(json.dumps({'tpr_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'expected_tpr':{'reference_init_nm':init,'rate_nm_per_ps':.005}}))
        return path,{'reference_init_nm':-10.7272,'rate_nm_per_ps':.005}
    def test_real_rounding_boundary(self):
        p,parsed=self.reference()
        r=validate_source_reference(p,parsed,1925.,-1.10225)
        self.assertLess(abs(r['residual_nm']),5e-6)
    def test_wrong_tpr_rejected(self):
        p,parsed=self.reference();p.write_bytes(b'CHANGED')
        with self.assertRaisesRegex(ValueError,'PRODUCER_HASH'):validate_source_reference(p,parsed,1925.,-1.10225)
    def test_real_reference_drift_rejected(self):
        p,parsed=self.reference()
        with self.assertRaisesRegex(ValueError,'PRODUCER_MISMATCH'):validate_source_reference(p,parsed,1925.,-1.1)
    def test_dump_drift_rejected(self):
        p,parsed=self.reference();parsed['reference_init_nm']=-10.8
        with self.assertRaisesRegex(ValueError,'DUMP_PRODUCER_DRIFT'):validate_source_reference(p,parsed,1925.,-1.10225)
    def test_rate_drift_rejected(self):
        p,parsed=self.reference();parsed['rate_nm_per_ps']=.006
        with self.assertRaisesRegex(ValueError,'RATE_DRIFT'):validate_source_reference(p,parsed,1925.,-1.10225)
    def test_nonfinite_rejected(self):
        p,parsed=self.reference()
        with self.assertRaisesRegex(ValueError,'NONFINITE'):validate_source_reference(p,parsed,float('nan'),-1.10225)
    def test_adjacent_corner_is_label_change_not_motion(self):
        points=[[0,0,0],[5,0,0],[5,10,0]]
        with self.assertRaises(ValueError):TrackedRoute(points).project_tracked([5.1,3,0],4.9,0)
        route=TransitionTrackedRoute(points);p=route.project_tracked([5.1,3,0],4.9,0)
        self.assertEqual(p.segment,1);self.assertTrue(route.tracking['arc_relocalized'])
        self.assertIn('ADJACENT_ARC_RELOCALIZATION',route.tracking['warnings'])
        self.assertEqual(route.project_tracked([5.1,3,0],p.s,p.segment),p)
    def test_same_segment_large_jump_still_rejected(self):
        with self.assertRaises(ValueError):TransitionTrackedRoute([[0,0,0],[10,0,0]]).project_tracked([8,0,0],1.,0)
    def test_branch_progress_mismatch_rejected(self):
        with self.assertRaises(ValueError):TransitionTrackedRoute([[0,0,0],[5,0,0],[5,10,0]]).project_tracked([5,3,0],1.,1)
    def test_exit_com_uses_anchor_images(self):
        xyz=np.array([[.1,0,0],[9.9,0,0],[4.,0,0]])
        self.assertTrue(np.allclose(group_com(xyz,[0,1],[1,1],np.eye(3)*10,0),[0.,0,0],atol=1e-12))
    def test_exit_anchor_membership_rejected(self):
        with self.assertRaises(ValueError):group_com(np.zeros((3,3)),[0,1],[1,1],np.eye(3),2)
    def test_index_preserves_old_groups(self):
        old='[ Protein_REF ]\n1 2 3\n'
        new=append_index_group(old,[4,5,6])
        self.assertTrue(new.startswith(old));self.assertIn('5 6 7',new)
        with self.assertRaises(ValueError):append_index_group(new,[7,8,9])
    def test_status_does_not_verify_runtime(self):
        import titan_path_receding_dispatch as d
        (self.root/'PLAN.json').write_text('{"cells":[{"cell":"sample"}]}')
        with patch.object(d,'verify',side_effect=AssertionError('must not run')),\
             patch.object(Path,'rglob',side_effect=AssertionError('must not walk')),\
             patch.object(sys,'argv',['dispatch','status','--root',str(self.root)]):
            self.assertEqual(d.main(),0)
    def test_bounded_collect_never_walks_runtime(self):
        from titan_path_result_export import collect_bounded
        cells=[{'cell':'sample','remaining_intervals':3}]
        folder=self.root/'outputs/sample';folder.mkdir(parents=True)
        (folder/'TRAJECTORY_RESULT.json').write_text('{"cell":"sample","status":"TECHNICAL_FAILURE","completed_intervals":2,"segments":[]}')
        (self.root/'PLAN.json').write_text(json.dumps({'cells':cells}));(self.root/'PROTOCOL.json').write_text('{}')
        with patch.object(Path,'rglob',side_effect=AssertionError('must not walk')):
            records=collect_bounded(self.root,cells,'SYNTHETIC.tar.gz')
        self.assertEqual(records[0]['status'],'TECHNICAL_FAILURE')
        ready=json.loads((self.root/'downloads/READY.json').read_text());self.assertFalse(ready['scientific_pass'])


if __name__=='__main__':unittest.main()
