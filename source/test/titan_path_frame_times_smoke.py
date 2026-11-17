"""Synthetic parser negatives plus retained real gmx-check fixtures; no MD."""
import os
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from titan_path_frame_times import parse_dump_frames
from titan_path_isolated_segment import check_time_coverage
from titan_path_controller import ControllerError


class Frames(unittest.TestCase):
    def assertRaises(self, exception, *args, **kwargs):
        return super().assertRaises((ValueError, ControllerError) if exception is ValueError else exception, *args, **kwargs)
    def frames(self):
        return [{'index': i, 'time_ps': 20+i*.2, 'step': 10000+100*i} for i in range(26)]

    def check(self, text, frames=None, end=25):
        return check_time_coverage(text, 20, end, complete=True, frame_records=frames)

    def test_no_last_not_extrapolated(self):
        with self.assertRaises(ValueError):
            self.check('Reading frame 0 time 20\nTime 26 0.2\n')

    def test_native_frames_supply_missing_last(self):
        result=self.check('Reading frame 0 time 20\nTime 26 0.2\n', self.frames())
        self.assertEqual((result['frames'],result['last_time_ps']), (26,25))

    def test_count_disagreement(self):
        with self.assertRaises(ValueError):self.check('Reading frame 0 time 20\nTime 27 0.2\n', self.frames())

    def test_explicit_last_disagreement(self):
        with self.assertRaises(ValueError):self.check('Reading frame 0 time 20\nLast frame 24 time 25\n', self.frames())

    def test_duplicate_first(self):
        with self.assertRaises(ValueError):self.check('Reading frame 0 time 20\n'*2, self.frames())

    def test_duplicate_last(self):
        with self.assertRaises(ValueError):self.check('Reading frame 0 time 20\n'+'Last frame 25 time 25\n'*2, self.frames())

    def test_bad_inventories(self):
        for key,value in [('time_ps',float('nan')),('time_ps',float('inf')),('time_ps',19),('index',3),('step',10000)]:
            f=self.frames();f[1][key]=value
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):self.check('Reading frame 0 time 20\n', f)

    def test_truncated_endpoint(self):
        with self.assertRaises(ValueError):self.check('Reading frame 0 time 20\n',self.frames()[:-1])

    def test_bad_interval(self):
        with self.assertRaises(ValueError):self.check('Reading frame 0 time 20\n',self.frames(),end=float('nan'))

    def test_xtc_dump(self):
        text='x.xtc frame 0:\n natoms= 39538 step= 10000 time=2.0000000e+01 prec=1000\n x[ 0]={0,0,0}\nx.xtc frame 1:\n natoms=39538 step=10100 time=20.2 prec=1000\n'
        self.assertEqual(len(parse_dump_frames(text.splitlines(True),'trajectory')),2)
        for bad in (text.replace('frame 1','frame 2'),text+'x.xtc frame 2:\n',text.replace('time=20.2','time=nan'),text.replace('step=10100','step=10000')):
            with self.assertRaises(ValueError):parse_dump_frames(bad.splitlines(True),'trajectory')

    def test_edr_dump(self):
        text='time: 20 step: 10000\nComponent Energy\ntime: 20.2 step: 10100\n'
        self.assertEqual(len(parse_dump_frames(text.splitlines(True),'energy')),2)
        for bad in ('',text.replace('20.2','20'),text+'time: malformed\n'):
            with self.assertRaises(ValueError):parse_dump_frames(bad.splitlines(True),'energy')

    def test_real_retained_check_logs(self):
        directory=os.environ.get('TITAN_FRAME_FIXTURES')
        if not directory:self.skipTest('Set TITAN_FRAME_FIXTURES to the five frozen real logs')
        paths=sorted(Path(directory).glob('*.log'));self.assertEqual(len(paths),5)
        import re
        for p in paths:
            text=p.read_text();start=float(re.findall(r'Reading frame\s+0\s+time\s+(\S+)',text)[0])
            with self.subTest(cell=p.stem):
                with self.assertRaises(ValueError):check_time_coverage(text,start,start+5,complete=True)
                # Synthetic inventory here tests compatibility, not native validity.
                f=[{'index':i,'time_ps':start+i*.2,'step':round(start/.002)+100*i} for i in range(26)]
                self.assertEqual(check_time_coverage(text,start,start+5,complete=True,frame_records=f)['frames'],26)


if __name__=='__main__':unittest.main()
