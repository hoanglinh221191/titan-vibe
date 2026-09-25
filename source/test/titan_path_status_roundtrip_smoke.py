"""Synthetic writer -> production Python parser -> native veto regression.

Not a new TITAN candidate search or evidence replay. Includes one analytically
reconstructed R1 defect vector, always labeled a regression fixture.
"""
import argparse
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import json

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_controller import parse_key_value_file
from titan_path_multiframe_replay import status_displacement

WRITER = VETO = None
CORE = [[x,y,z] for x in (-2.,2.) for y in (-2.,2.) for z in (-2.,2.)]


class RoundtripSmoke(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix='NAV6C-R2-SYNTHETIC-')
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def write_parse(self, vector):
        proc = subprocess.run([str(WRITER),*[format(x,'.17g') for x in vector]],capture_output=True,text=True,timeout=10)
        self.assertEqual(proc.returncode,0,proc.stderr)
        self.assertNotIn(',',proc.stdout)
        path = self.root/'status.txt'
        path.write_text(proc.stdout)
        parsed = parse_key_value_file(path)
        delta = status_displacement(parsed)
        self.assertEqual([x.hex() for x in delta],[x.hex() for x in vector])
        return delta,parsed

    def veto(self,delta,old=(1.,0.,0.),recovery=False,extra=''):
        flat = [*old,*delta,*[x for p in CORE for x in p],*[x for p in CORE for x in p],
                0,0,0,5,0,0,10,0,0,100,100,100,3,5 if recovery else 1,0,12]
        packet = self.root/'packet.txt'
        packet.write_text(f'TITAN_MULTIFRAME_INPUT_V1 8 3 1 1 -1 {int(recovery)}\n'+
                          ' '.join(format(x,'.17g') for x in flat)+'\n'+extra)
        return subprocess.run([str(VETO),str(packet)],capture_output=True,text=True,timeout=10)

    def angle_case(self,cap,offset,accepted):
        angle=math.radians(cap+offset)
        delta,_=self.write_parse([.25*math.cos(angle),-.25*math.sin(angle),0.])
        proc=self.veto(delta,recovery=cap==60)
        self.assertEqual(proc.returncode,0,proc.stderr)
        result=json.loads(proc.stdout)
        self.assertEqual(result['accepted'],accepted)
        self.assertFalse(result['md_authorized'])
        self.assertEqual(result['reason'],'ACCEPT_NATIVE_DIRECTION_NO_MD' if accepted else 'STOP_TURN_CAP')
        return result

    def test_advance_exact_cap(self):self.angle_case(20,0,True)
    def test_recovery_exact_cap(self):self.angle_case(60,0,True)
    def test_advance_inside_tolerance(self):self.angle_case(20,.5e-6,True)
    def test_recovery_inside_tolerance(self):self.angle_case(60,.5e-6,True)
    def test_advance_outside_tolerance(self):self.angle_case(20,2e-6,False)
    def test_recovery_outside_tolerance(self):self.angle_case(60,2e-6,False)

    def test_reconstructed_r1_rounding_defect(self):
        old=[-.21753551255192122,-.71343728176110943,.66609724948545801]
        vector=[.021931896723951003,-.24088733401835033,-.063184525127769]
        delta,_=self.write_parse(vector)
        fixed=json.loads(self.veto(delta,old,True).stdout)
        rounded=[float(format(x,'.8f')) for x in vector]
        legacy=json.loads(self.veto(rounded,old,True).stdout)
        self.assertTrue(fixed['accepted'])
        self.assertEqual(legacy['reason'],'STOP_TURN_CAP')
        self.assertGreater(legacy['turn_deg'],60+1e-6)

    def test_missing_precision_contract_rejected(self):
        _,status=self.write_parse([.25,0.,0.]);del status['NUMERIC_SERIALIZATION']
        with self.assertRaises(ValueError):status_displacement(status)

    def test_nonfinite_components_rejected(self):
        _,status=self.write_parse([.25,0.,0.])
        for value in ('nan','inf','-inf','1e999'):
            with self.subTest(value=value):
                status['NEXT_VECTOR_X_A']=value
                with self.assertRaises(ValueError):status_displacement(status)

    def test_zero_and_wrong_norm_rejected(self):
        for x in (0.,1.):
            status={'NUMERIC_SERIALIZATION':'FLOAT64_MAX_DIGITS10',
                    'NEXT_VECTOR_X_A':str(x),'NEXT_VECTOR_Y_A':'0','NEXT_VECTOR_Z_A':'0'}
            with self.assertRaises(ValueError):status_displacement(status)

    def test_missing_component_rejected(self):
        _,status=self.write_parse([.25,0.,0.]);del status['NEXT_VECTOR_Z_A']
        with self.assertRaises(KeyError):status_displacement(status)

    def test_smallest_normal_component_and_negative_zero_preserved(self):
        self.write_parse([.25,float.fromhex('0x1p-1022'),-0.])

    def test_trailing_packet_data_rejected(self):
        delta,_=self.write_parse([.25,0.,0.])
        self.assertEqual(self.veto(delta,extra='1').returncode,2)

    def test_native_nonfinite_packet_rejected(self):
        self.assertEqual(self.veto([math.nan,0.,0.]).returncode,2)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--writer',required=True,type=Path);parser.add_argument('--veto',required=True,type=Path)
    args=parser.parse_args();WRITER=args.writer.resolve();VETO=args.veto.resolve()
    unittest.main(argv=[sys.argv[0]],verbosity=2)
