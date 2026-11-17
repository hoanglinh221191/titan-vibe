import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_time_match import match_frames


class Times(unittest.TestCase):
    def fixture(self,start):
        p=np.arange(round(start/.002),round((start+5)/.002)+1)*.002
        r=np.array([float(f'{float(np.float32(t)):.5f}') for t in p[::100]])
        return r,p
    def test_boundaries(self):
        for start in (0.,195.,1020.,2045.,4095.,8190.):
            with self.subTest(start=start):
                r,p=self.fixture(start);i,e=match_frames(r,p,start,start+5)
                self.assertEqual(len(i),26);self.assertFalse(e['interpolated'])
    def test_missing_pull(self):
        r,p=self.fixture(2045)
        with self.assertRaises(ValueError):match_frames(r,np.delete(p,12),2045,2050)
    def test_duplicate_pull(self):
        r,p=self.fixture(2045);p[12]=p[11]
        with self.assertRaises(ValueError):match_frames(r,p,2045,2050)
    def test_missing_frame(self):
        r,p=self.fixture(2045)
        with self.assertRaises(ValueError):match_frames(np.delete(r,12),p,2045,2050)
    def test_wrong_step_not_nearest_matched(self):
        r,p=self.fixture(2045);r[5]+=.002
        with self.assertRaises(ValueError):match_frames(r,p,2045,2050)
    def test_offgrid_pull(self):
        r,p=self.fixture(2045);p[5]+=.00001
        with self.assertRaises(ValueError):match_frames(r,p,2045,2050)
    def test_nan(self):
        r,p=self.fixture(2045);r[4]=np.nan
        with self.assertRaises(ValueError):match_frames(r,p,2045,2050)
    def test_time_alias_guard(self):
        r,p=self.fixture(65535)
        with self.assertRaisesRegex(ValueError,'PRECISION'):match_frames(r,p,65535,65540)


if __name__=='__main__':unittest.main()
