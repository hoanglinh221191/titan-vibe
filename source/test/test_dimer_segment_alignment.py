"""Independent small fixtures for scientific alignment rejection rules."""
import unittest
from dimer_segment_alignment import segment_identity,closed_no_exchange,determinant

class AlignmentRules(unittest.TestCase):
 def test_known_segment_offset(self):
  offset,matches=segment_identity([(0.,10.,300.,2.),(20.,11.,301.,3.)],{(10.,300.,2.):[(5,160.)],(11.,301.,3.):[(6,180.)]})
  self.assertEqual(offset,160.);self.assertEqual(matches,[(0.,5,160.),(20.,6,180.)])
 def test_missing_row_rejected(self):
  with self.assertRaises(ValueError):segment_identity([(20.,10.,300.,2.)],{})
 def test_duplicate_value_identity_rejected(self):
  with self.assertRaises(ValueError):segment_identity([(20.,10.,300.,2.)],{(10.,300.,2.):[(1,20.),(2,40.)]})
 def test_mixed_offsets_rejected(self):
  with self.assertRaises(ValueError):segment_identity([(0.,10.,300.,2.),(20.,11.,301.,3.)],{(10.,300.,2.):[(1,160.)],(11.,301.,3.):[(2,200.)]})
 def test_temperature_participates_in_identity(self):
  with self.assertRaises(ValueError):segment_identity([(0.,10.,300.,2.)],{(10.,301.,2.):[(0,0.)]})
 def test_closed_exchange_boundaries(self):
  for times in [[50],[100],[75],[49,50,101]]:self.assertFalse(closed_no_exchange(times,50,100))
  self.assertTrue(closed_no_exchange([49,101],50,100));self.assertTrue(closed_no_exchange([],50,100))
 def test_tilted_cell_volume(self):
  self.assertAlmostEqual(determinant([2.,0,0,1,3.,0,.3,.4,4.]),24.)
 def test_exact_cadence_keeps_holes(self):
  energies={float(t):t for t in range(0,201,20)}
  self.assertEqual([energies.get(t) for t in [0.,50.,100.,150.,200.]],[0,None,100,None,200])

if __name__=='__main__':unittest.main(verbosity=2)
