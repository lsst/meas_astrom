
from __future__ import division, print_function, absolute_import

from copy import copy
import unittest

import numpy as np

from lsst.meas.astrom.pessimistic_pattern_matcher_b_3D \
    import PessimisticPatternMatcherB
from lsst.log import Log

__deg_to_rad__ = np.pi/180


class TestPythonOptimisticPatternMatcherB(unittest.TestCase):

    """Unittest suite for the python implimentation of
    Optimistic Pattern Matcher B.
    """

    def setUp(self):
        np.random.seed(12345)

        n_points = 1000
        self.reference_catalog = np.empty((n_points, 4))
        cos_theta_array = np.random.uniform(
            np.cos(np.pi/2 + 0.5*__deg_to_rad__),
            np.cos(np.pi/2 - 0.5*__deg_to_rad__), size=n_points)
        sin_theta_array = np.sqrt(1 - cos_theta_array**2)
        phi_array = np.random.uniform(-0.5, 0.5, size=n_points)*__deg_to_rad__
        self.reference_catalog[:, 0] = sin_theta_array*np.cos(phi_array)
        self.reference_catalog[:, 1] = sin_theta_array*np.sin(phi_array)
        self.reference_catalog[:, 2] = cos_theta_array
        self.reference_catalog[:, 3] = (
            np.random.power(1.2, size=n_points)*4 + 20)
        self.source_catalog = copy(self.reference_catalog)
        self.log = Log()

    def testInit(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)

    def testConstructAndMatchPattern(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)

        (pattern_list, src_cand, shift_rot_matrix, cos_theta,
         sin_theta) = self.pyPPMb._construct_pattern_and_shift_rot_matrix(
            self.source_catalog[:6, :3], 6, np.cos(np.radians(60. / 3600.)),
            np.cos(np.radians(1.0)) ** 2, np.radians(5./3600.))
        self.assertGreater(len(pattern_list), 0)
        self.assertEqual(pattern_list[0], 0)
        self.assertEqual(pattern_list[1], 1)
        self.assertEqual(pattern_list[2], 2)
        self.assertEqual(pattern_list[3], 3)
        self.assertEqual(pattern_list[4], 4)
        self.assertEqual(pattern_list[5], 5)

        (pattern_list, src_cand, shift_rot_matrix, cos_theta,
         sin_theta) = self.pyPPMb._construct_pattern_and_shift_rot_matrix(
            self.source_catalog[:9, :3], 6, np.cos(np.radians(60. / 3600.)),
            np.cos(np.radians(1.0)) ** 2, np.radians(5./3600.))
        self.assertGreater(len(pattern_list), 0)
        self.assertEqual(pattern_list[0], 0)
        self.assertEqual(pattern_list[1], 1)
        self.assertEqual(pattern_list[2], 2)
        self.assertEqual(pattern_list[3], 3)
        self.assertEqual(pattern_list[4], 4)
        self.assertEqual(pattern_list[5], 5)

        (pattern_list, src_cand, shift_rot_matrix, cos_theta,
         sin_theta) = self.pyPPMb._construct_pattern_and_shift_rot_matrix(
            self.source_catalog[[2, 4, 8, 16, 32, 64], :3], 6,
            np.cos(np.radians(60. / 3600.)), np.cos(np.radians(1.0)) ** 2,
            np.radians(5./3600.))
        self.assertEqual(pattern_list[0], 2)
        self.assertEqual(pattern_list[1], 4)
        self.assertEqual(pattern_list[2], 8)
        self.assertEqual(pattern_list[3], 16)
        self.assertEqual(pattern_list[4], 32)
        self.assertEqual(pattern_list[5], 64)

    def testMatchPerfect(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)

        matches, distances, pattern_idx, shift = self.pyPPMb.match(
            source_catalog=self.source_catalog, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)
        self.assertEqual(len(matches), len(self.reference_catalog))
        self.assertTrue(np.all(distances < 0.01/3600.0*__deg_to_rad__))

    def testOptimisticMatch(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)

        matches, distances, pattern_idx, shift = self.pyPPMb.match_optimistic(
            source_catalog=self.source_catalog, n_check=9, n_match=6,
            max_n_patterns=100, max_shift=60., max_rotation=6.0, max_dist=5.,
            min_matches=30, pattern_skip_array=None)
        self.assertEqual(len(matches), len(self.reference_catalog))
        self.assertTrue(np.all(distances < 0.01/3600.0*__deg_to_rad__))

    def testMatchSkip(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)

        matches, distances, pattern_idx, shift = self.pyPPMb.match(
            source_catalog=self.source_catalog, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=np.array([0]))
        self.assertEqual(len(matches), len(self.reference_catalog))
        self.assertTrue(np.all(distances < 0.01/3600.0*__deg_to_rad__))

    def testMatchMoreSources(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:500, :3],
            log=self.log)

        matches, distances, pattern_idx, shift = self.pyPPMb.match(
            source_catalog=self.source_catalog, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60.0, max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)
        self.assertEqual(len(matches), len(self.reference_catalog[:500]))
        self.assertTrue(np.all(distances < 0.01/3600.0*__deg_to_rad__))

    def testMatchMoreReferences(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)

        matches, distances, pattern_idx, shift = self.pyPPMb.match(
            source_catalog=self.source_catalog[:500], n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=1.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)
        self.assertEqual(len(matches), len(self.reference_catalog[:500]))
        self.assertTrue(np.all(distances < 0.01/3600.0*__deg_to_rad__))

    def testShift(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)
        theta = np.radians(45.0 / 3600.)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        theta_rotation = np.array(
            [[cos_theta, -sin_theta, 0.],
             [sin_theta,  cos_theta, 0.],
             [       0.,         0., 1.]])

        self.source_catalog[:, :3] = np.dot(
            theta_rotation, self.source_catalog[:, :3].transpose()).transpose()

        matches, distances, pattern_idx, shift = self.pyPPMb.match(
            source_catalog=self.source_catalog, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60, max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)

        self.assertEqual(len(matches), len(self.reference_catalog))
        self.assertTrue(np.all(distances < 0.1/3600.0 * __deg_to_rad__))

    def testRotation(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)
        phi = 2.5*__deg_to_rad__
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        phi_rotation = np.array(
            [[1.,      0.,       0.],
             [0., cos_phi, -sin_phi],
             [0., sin_phi,  cos_phi]])

        self.source_catalog[:, :3] = np.dot(
            phi_rotation, self.source_catalog[:, :3].transpose()).transpose()

        matches, distances, pattern_idx, shift = self.pyPPMb.match(
            source_catalog=self.source_catalog, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60, max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)

        self.assertEqual(len(matches), len(self.reference_catalog))
        self.assertTrue(np.all(distances < 0.1 / 3600.0 * __deg_to_rad__))

    def testShiftRotation(self):
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_catalog=self.reference_catalog[:, :3],
            log=self.log)
        theta = np.radians(45.0 / 3600.)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        theta_rotation = np.array(
            [[cos_theta, -sin_theta, 0.],
             [sin_theta,  cos_theta, 0.],
             [       0.,         0., 1.]])

        phi = 2.5*__deg_to_rad__
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        phi_rotation = np.array(
            [[1.,      0.,       0.],
             [0., cos_phi, -sin_phi],
             [0., sin_phi,  cos_phi]])

        shift_rot_matrix = np.dot(theta_rotation, phi_rotation)

        self.source_catalog[:, :3] = np.dot(
            shift_rot_matrix,
            self.source_catalog[:, :3].transpose()).transpose()

        matches, distances, pattern_idx, shift = self.pyPPMb.match(
            source_catalog=self.source_catalog, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)

        self.assertEqual(len(matches), len(self.reference_catalog))
        self.assertTrue(np.all(distances < 0.1/3600.0*__deg_to_rad__))


if __name__ == '__main__':
    unittest.main()
