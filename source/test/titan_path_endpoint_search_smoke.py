import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from titan_path_endpoint_search import surface_gaps, endpoint_candidates, endpoint_reached, pdb_box


class EndpointTests(unittest.TestCase):
    def test_surface_not_center(self):
        self.assertAlmostEqual(surface_gaps([[5, 0, 0]], [[0, 0, 0]], [2], np.eye(3)*40)[0], 3)
    def test_periodic_neighbor(self):
        self.assertAlmostEqual(surface_gaps([[39, 0, 0]], [[0, 0, 0]], [2], np.eye(3)*40)[0], -1)
    def test_translation_invariance(self):
        a = surface_gaps([[5, 1, 2]], [[0, 0, 0]], [2], np.eye(3)*40)
        b = surface_gaps([[45, 1, 2]], [[0, 0, 0]], [2], np.eye(3)*40)
        np.testing.assert_allclose(a, b)
    def test_radii(self):
        with self.assertRaises(ValueError): surface_gaps([[0, 0, 0]], [[0, 0, 0]], [-1], np.eye(3)*40)
    def test_nan(self):
        with self.assertRaises(ValueError): surface_gaps([[np.nan, 0, 0]], [[0, 0, 0]], [1], np.eye(3)*40)
    def test_near_original_endpoint(self):
        shells = list(endpoint_candidates([1, 2, 3], [1, 0, 0], max_radius_A=1.))
        self.assertEqual([s for s, _ in shells], [0, .5, 1.])
        for radius, points in shells:
            np.testing.assert_allclose(np.linalg.norm(points-[1, 2, 3], axis=1), radius, atol=1e-12)
            self.assertTrue(np.all(points[:, 0] >= 1))
    def test_projection_only_not_success(self):
        self.assertFalse(endpoint_reached([4, 0, 0], [0, 0, 0], [10, 10, 10], 10))
    def test_all_frames(self):
        self.assertFalse(endpoint_reached([0, 0, 0], [0, 0, 0], [10, 10, 8], 10))
    def test_reached(self):
        self.assertTrue(endpoint_reached([.2, 0, 0], [0, 0, 0], [9.5, 9.8, 10], 10))
    def test_nonfinite_reached(self):
        self.assertFalse(endpoint_reached([np.nan, 0, 0], [0, 0, 0], [10, 10, 10], 10))


if __name__ == '__main__': unittest.main()
