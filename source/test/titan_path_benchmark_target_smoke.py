"""Independent analytic/brute-force synthetic checks; no MD."""
import itertools
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import titan_path_benchmark_target as target
from titan_path_b1_start_coverage import proof


class TargetTests(unittest.TestCase):
    def test_orthorhombic_exact_and_image_invariance(self):
        box = np.diag([10., 12., 14.]); d = np.array([[9, 0, 0], [0, -11, 0], [0, 0, 13.]])
        expected = np.array([[-1,0,0], [0,1,0], [0,0,-1]])
        np.testing.assert_allclose(target.minimum_image_vectors(d, box), expected)
        np.testing.assert_allclose(target.minimum_image_vectors(d + np.array([3, -2, 5]) @ box, box), expected)

    def test_skew_box_against_independent_brute_force(self):
        rng = np.random.default_rng(7543)
        box = np.array([[10.,0,0], [4.9,9,0], [4.8,4.4,8.]])
        d = rng.uniform(-15,15,(60,3))
        got = np.linalg.norm(target.minimum_image_vectors(d,box),axis=1)
        shifts = np.array(list(itertools.product(range(-5,6),repeat=3))) @ box
        oracle = np.min(np.linalg.norm(d[:,None,:]-shifts[None,:,:],axis=2),axis=1)
        np.testing.assert_allclose(got,oracle,rtol=0,atol=1e-12)

    def test_invalid_box_or_numeric_rejected(self):
        for box in (np.zeros((3,3)), np.diag([-1.,1,1]), np.eye(3)*float('nan')):
            with self.assertRaises(ValueError): target.minimum_image_vectors([[1.,2,3]],box)
        with self.assertRaises(ValueError): target.minimum_image_vectors([[float('nan'),0,0]],np.eye(3))

    def test_persistence_requires_five_ps_not_one_frame(self):
        state = {}
        for i in range(26):
            state = target.persistence_step(state,1000+i*.2,True)
            self.assertEqual(state['target_reached'],i==25)

    def test_persistence_reset_and_duplicate_or_gap_rejected(self):
        state = target.persistence_step({},0.,True)
        state = target.persistence_step(state,.2,False)
        state = target.persistence_step(state,.4,True)
        self.assertEqual(state['condition_since_ps'],.4)
        for t in (.4,.3,.8):
            with self.assertRaises(ValueError): target.persistence_step(state,t,True)

    def test_registration_removes_rigid_motion_but_not_ligand_motion(self):
        xyz = np.zeros((2615,3)); core=np.array([[x,y,z] for x in range(3) for y in range(3) for z in range(9)],dtype=float)
        xyz[:81]=core;xyz[2603:]=[3,2,1]
        a=target.registered_coms(xyz,list(range(81)),core,[12.]*12)
        rotation=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
        b=target.registered_coms(xyz@rotation.T+[30,40,50],list(range(81)),core,[12.]*12)
        np.testing.assert_allclose(a,b,atol=1e-12)
        xyz[2603:,0]+=10
        c=target.registered_coms(xyz,list(range(81)),core,[12.]*12)
        np.testing.assert_allclose(np.linalg.norm(c-a,axis=1),10,atol=1e-12)

    def test_target_not_met_initially_and_met_after_contact_loss(self):
        xyz = np.zeros((2615,3)); core=np.array([[x,y,z] for x in range(3) for y in range(3) for z in range(9)],dtype=float)
        xyz[:81]=core;xyz[2603:]=[3,2,1]
        parameters=[{'source_atom_id_0based':i,'molecule':'mol1' if i<2603 else 'mol2',
                     'mass_u':12. if i<2609 else 1.008} for i in range(2615)]
        frozen=target.define_target(xyz,np.eye(3)*100,list(range(81)),core,parameters)
        initial=target.measure_target(frozen,xyz,np.eye(3)*100,list(range(81)),core)
        self.assertFalse(initial['condition_met']);self.assertEqual(initial['retained_fraction'],1.)
        xyz[2603:,0]+=12
        end=target.measure_target(frozen,xyz,np.eye(3)*100,list(range(81)),core)
        self.assertTrue(end['condition_met']);self.assertEqual(end['retained_fraction'],0.)
        self.assertFalse(frozen['policy']['solvent_unbinding_claimed'])
        self.assertFalse(frozen['policy']['route_progress_used'])


class CoverageTests(unittest.TestCase):
    def test_entire_probe_ball_blocked_and_prefix_diagnostic(self):
        result = proof([[0,0,0],[3,0,0],[4,1,0]],np.array([-1.,.2,0]))
        self.assertTrue(result['all_directions_in_probe_ball_have_zero_progress'])
        self.assertEqual(result['unmodified_tangent_delta_progress_A'],0.)
        self.assertAlmostEqual(result['proposed_2A_prefix_tangent_delta_progress_A'],.25)

    def test_near_head_can_make_progress(self):
        result=proof([[0,0,0],[3,0,0],[4,1,0]],np.array([-.1,.2,0]))
        self.assertFalse(result['all_directions_in_probe_ball_have_zero_progress'])
        self.assertAlmostEqual(result['unmodified_tangent_delta_progress_A'],.15)

    def test_nearby_competing_segment_invalidates_global_proof(self):
        result=proof([[0,0,0],[1,0,0],[-1,.1,0]],np.array([-.5,.2,0]))
        self.assertFalse(result['all_directions_in_probe_ball_have_zero_progress'])

    def test_proof_agrees_with_dense_sphere_oracle(self):
        from titan_path_observation import Polyline
        points=[[0,0,0],[3,0,0],[4,1,0]];com=np.array([-1.,.2,0])
        self.assertTrue(proof(points,com)['all_directions_in_probe_ball_have_zero_progress'])
        rng=np.random.default_rng(2234); directions=rng.normal(size=(1000,3));directions/=np.linalg.norm(directions,axis=1)[:,None]
        line=Polyline(points)
        self.assertTrue(all(line.project(com+.25*d).s==0 for d in directions))


if __name__=='__main__': unittest.main(verbosity=2)
