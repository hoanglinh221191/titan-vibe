"""R6 continuation must not override a remaining HOLD or upstream hash drift."""
import csv,copy,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_multiframe_feedback_smoke import fixture,CORE,IDENTITY
from titan_path_feedback import seal,unseal,FeedbackError
import titan_path_r6_resume as resume


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.cell={'cell':'case','route':'route','route_points':'points',
            'route_id':'SYNTHETIC','reference_sha256':'a'*64}
        self.before,self.raw,self.whole,self.pull=fixture()
        old=unseal(self.before);old.update(status='HOLD',checkpoint_sha256='c'*64,time_ps=5.,progress_A=None,distance_A=None)
        self.historical=seal(old)
        folder=self.root/'inputs/cells/case';folder.mkdir(parents=True)
        with (folder/'baseline_parameters.tsv').open('w',newline='') as f:
            writer=csv.writer(f,delimiter='\t');writer.writerow(['mass_u']);writer.writerows([[12]]*2615)
        self.writes={}
    def tearDown(self):self.tmp.cleanup()
    def read(self,path):
        if path.name=='BINDING.json':return {'files':{k:{'cluster_path':'/r6/'+k,'sha256':'c'*64} for k in ('checkpoint','gro','tpr','none','whole','pullx')},'old_direction':[1,0,0]}
        if path.name=='SEGMENT_RESULT.json':return {'output_hashes':{'xtc':'d'*64},'boundary':[0.]*12}
        if path.name=='STATE.json':return self.historical
        if path.name=='PREVIOUS_STATE.json':return self.before
        if path.name=='FRAME_DEFINITION.json':return {'selected_source_ids':list(range(81)),'reference_xyz_A':CORE.tolist()}
        if path.name=='PROTOCOL.json':return {'route_tracking':'bounded-adjacent-route-tracking-v1'}
        if path.name=='points':return [[0,0,0],[10,0,0]]
        raise AssertionError(path)
    def call(self,drift=False):
        def sha(p):return ('0' if drift else 'c')*64 if Path(p).as_posix().startswith('/r6/') else 'a'*64
        with patch.object(resume,'read',side_effect=self.read),patch.object(resume,'sha',side_effect=sha),\
             patch.object(resume,'write',side_effect=lambda p,v:self.writes.update({str(p):v})),\
             patch.object(resume,'read_gro',side_effect=lambda p:self.raw if p.name=='none' else self.whole),\
             patch.object(resume,'read_xvg',return_value=self.pull):
            return resume.bind_resume(self.root,self.cell)
    def test_reaudit_records_new_state_not_old_rewrite(self):
        old=copy.deepcopy(self.historical);context=self.call()
        self.assertEqual(context['state']['status'],'READY');self.assertEqual(context['state']['segments'],0)
        self.assertEqual(context['state']['time_ps'],5.)
        self.assertEqual(self.historical,old);self.assertTrue(any(p.endswith('ROWS.json') for p in self.writes))
    def test_bad_com_remains_hold(self):
        self.pull[10,9]+=1
        with self.assertRaises(FeedbackError):self.call()
        self.assertFalse(any(p.endswith('BOUND.json') for p in self.writes))
    def test_image_break_remains_hold(self):
        self.raw[10]['xyz'][0,0]+=100
        with self.assertRaises(FeedbackError):self.call()
    def test_source_hash_drift_blocks(self):
        with self.assertRaises(FeedbackError):self.call(drift=True)
        self.assertFalse(self.writes)
    def test_bad_previous_state_blocks(self):
        body=unseal(self.before);body['status']='HOLD';self.before=seal(body)
        with self.assertRaises(FeedbackError):self.call()


if __name__=='__main__':unittest.main()
