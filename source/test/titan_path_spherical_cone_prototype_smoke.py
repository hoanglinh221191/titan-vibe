"""Constraint-only spherical-domain invariants, no MD."""
import math
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_spherical_cone_prototype import in_cone,enumerate_domain
from titan_path_corner_reward_prototype import CornerPolicy,Polyline,lookahead_goal


class ConeTests(unittest.TestCase):
    def test_axis_inside(self):self.assertTrue(in_cone([0,0,2],[0,0,1],20))
    def test_opposite_outside(self):self.assertFalse(in_cone([0,0,-1],[0,0,1],80))
    def test_boundary_inclusive(self):self.assertTrue(in_cone([math.sin(math.radians(20)),0,math.cos(math.radians(20))],[0,0,1],20))
    def test_invalid_angle(self):
        for a in (0,-1,181,float('nan')):
            with self.assertRaises(ValueError):in_cone([1,0,0],[1,0,0],a)
    def test_zero_vector(self):
        with self.assertRaises(ValueError):in_cone([0,0,0],[1,0,0],20)
    def test_rotation_covariance_of_constraint(self):
        R=np.array([[0,1,0],[-1,0,0],[0,0,1]])
        for x in ([1,.2,0],[0,1,.1],[-1,0,0]):
            self.assertEqual(in_cone(x,[1,0,0],40),in_cone(R@x,R@np.array([1,0,0]),40))
    def test_empty_intersection_no_fallback(self):
        line=Polyline([[-2,0,0],[0,0,0],[0,3,0]]);routes=[line]*3
        goals=[lookahead_goal(line,2,CornerPolicy())]*3
        result=enumerate_domain(np.array([[0,0,10]]),np.array([[.2,-.7,0]]),[1],routes,[2]*3,[1,0,0],goals,20)
        self.assertTrue(result['pairwise_empty_intersection_proved']);self.assertEqual(result['candidates'],[])
        self.assertFalse(result['md_authorized'])
    def test_domain_candidates_all_constraints(self):
        line=Polyline([[-2,0,0],[0,0,0],[0,3,0]]);routes=[line]*3
        goals=[lookahead_goal(line,2,CornerPolicy())]*3
        result=enumerate_domain(np.array([[0,0,10]]),np.array([[.2,-.7,0]]),[1],routes,[2]*3,[1,0,0],goals,80)
        self.assertGreater(len(result['candidates']),0)
        for c in result['candidates']:
            self.assertLessEqual(c['old_turn_deg'],20+1e-8)
            self.assertLessEqual(max(c['path_cone_angles_deg']),80+1e-8)
            self.assertGreaterEqual(c['clearance_A'],1.5)
            self.assertGreater(min(c['goal_approach_A']),1e-6)
        self.assertFalse(result['progress_reward_used'])


if __name__=='__main__':unittest.main()
