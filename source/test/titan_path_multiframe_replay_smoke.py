"""Synthetic package/collector guard tests, not native replay evidence."""
from pathlib import Path
import contextlib
import io
import json
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tools"))
import titan_path_multiframe_replay as r

class ReplaySmoke(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix="NAV6C-SYNTHETIC-");self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
    def collect(self,records,error=None):
        with contextlib.redirect_stdout(io.StringIO()):return r.collect(self.root,records,error)
    def records(self):return [{"cell":"synthetic_"+str(i),"outcome":"STOP_FRAME_DISAGREEMENT","accepted":False} for i in range(5)]
    def test_all_stops_can_be_technical_complete_not_scientific(self):
        self.assertTrue(self.collect(self.records()))
        data=r.read(self.root/"generated/REPORT.json")
        self.assertEqual(data["accepted_count"],0);self.assertFalse(data["scientific_pass"])
        self.assertTrue((self.root/"downloads/READY.json").is_file())
    def test_partial_results_get_failure_bundle_only(self):
        self.assertFalse(self.collect(self.records()[:2]))
        self.assertTrue((self.root/"downloads/FAILED.json").is_file())
        self.assertFalse((self.root/"downloads/READY.json").exists())
    def test_technical_error_never_ready(self):
        rows=self.records();rows[2]["outcome"]="TECHNICAL_ERROR"
        self.assertFalse(self.collect(rows))
    def test_turn_cap_is_a_valid_terminal_not_frame_disagreement(self):
        rows=self.records();rows[0]['outcome']='STOP_TURN_CAP'
        self.assertTrue(self.collect(rows))
    def test_build_error_collected(self):self.assertFalse(self.collect([],"synthetic build error"))
    def test_duplicate_cell_not_complete(self):
        rows=self.records();rows[-1]["cell"]=rows[0]["cell"]
        self.assertFalse(self.collect(rows))
    def test_unknown_outcome_not_complete(self):
        rows=self.records();rows[0]["outcome"]="MADE_UP_PASS"
        self.assertFalse(self.collect(rows))
    def test_fake_accept_count_not_complete(self):
        rows=self.records();rows[0]["accepted"]=True
        self.assertFalse(self.collect(rows))
    def test_archive_contains_every_manifest_member(self):
        self.collect(self.records());ready=r.read(self.root/"downloads/READY.json")
        with tarfile.open(self.root/"downloads"/ready["bundle"]) as a:
            manifest=json.loads(a.extractfile("MANIFEST.json").read())
            self.assertEqual(set(a.getnames()),set(manifest)|{"MANIFEST.json"})
    def test_duplicate_collection_refuses(self):
        self.collect(self.records())
        with self.assertRaises(FileExistsError):self.collect(self.records())
    def test_path_traversal_rejected(self):
        with self.assertRaises(ValueError):r.path_at(self.root,"../outside")
    def test_absolute_path_rejected(self):
        with self.assertRaises(ValueError):r.path_at(self.root,"C:/outside")
    def test_sweep_interior_not_endpoint_only(self):
        self.assertAlmostEqual(r.sweep([[0,0,0]],[[-2,1,0]],[4,0,0]),1.)
    def test_sweep_reverse(self):
        self.assertEqual(r.sweep([[0,0,0]],[[-2,1,0]],[4,0,0]),r.sweep([[0,0,0]],[[2,1,0]],[-4,0,0]))
    def test_submission_guard_blocks_repeat_without_sbatch(self):
        (self.root/"generated/submission.guard").mkdir(parents=True)
        with patch.object(r,"verify",return_value={}),patch.object(r.subprocess,"run") as run:
            with self.assertRaises(FileExistsError):r.submit(self.root)
            run.assert_not_called()
    def test_no_ready_without_archive_audit(self):
        with patch.object(r,"verify",return_value={}):
            with self.assertRaisesRegex(ValueError,"READY absent"):r.audit(self.root)

if __name__=="__main__":unittest.main()
