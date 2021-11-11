#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

from copy import copy
import unittest
import logging

import numpy as np

from lsst.meas.astrom.pessimistic_pattern_matcher_b_3D \
    import PessimisticPatternMatcherB

__deg_to_rad__ = np.pi/180


class TestPessimisticPatternMatcherB(unittest.TestCase):

    """Unittest suite for the Pessimistic Pattern Matcher B.
    """

    def setUp(self):
        np.random.seed(12345)

        n_points = 1000
        # reference_obj_array is a number array representing
        # 3D points randomly draw on a 1 sq deg patch.
        self.reference_obj_array = np.empty((n_points, 4))
        cos_theta_array = np.random.uniform(
            np.cos(np.pi/2 + 0.5*__deg_to_rad__),
            np.cos(np.pi/2 - 0.5*__deg_to_rad__), size=n_points)
        sin_theta_array = np.sqrt(1 - cos_theta_array**2)
        phi_array = np.random.uniform(-0.5, 0.5, size=n_points)*__deg_to_rad__
        self.reference_obj_array[:, 0] = sin_theta_array*np.cos(phi_array)
        self.reference_obj_array[:, 1] = sin_theta_array*np.sin(phi_array)
        self.reference_obj_array[:, 2] = cos_theta_array
        self.reference_obj_array[:, 3] = (
            np.random.power(1.2, size=n_points)*4 + 20)

        # Our initial source catalog is a straight copy of the reference
        # array at first. In some of the tests we will add rotations and
        # shifts to the data in order to test the input and outputs of our
        # matcher.
        self.source_obj_array = copy(self.reference_obj_array)
        self.log = logging.getLogger(__name__)

    def testConstructPattern(self):
        """ Test that a specified pattern can be found in the reference
        data and that the explicit ids match.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)

        pattern_struct = self.pyPPMb._construct_pattern_and_shift_rot_matrix(
            self.source_obj_array[:6, :3], 6, np.cos(np.radians(60. / 3600.)),
            np.cos(np.radians(1.0)) ** 2, np.radians(5./3600.))
        pattern_list = pattern_struct.ref_candidates
        self.assertGreater(len(pattern_list), 0)
        self.assertEqual(pattern_list[0], 0)
        self.assertEqual(pattern_list[1], 1)
        self.assertEqual(pattern_list[2], 2)
        self.assertEqual(pattern_list[3], 3)
        self.assertEqual(pattern_list[4], 4)
        self.assertEqual(pattern_list[5], 5)

        pattern_struct = self.pyPPMb._construct_pattern_and_shift_rot_matrix(
            self.source_obj_array[:9, :3], 6, np.cos(np.radians(60. / 3600.)),
            np.cos(np.radians(1.0)) ** 2, np.radians(5./3600.))
        pattern_list = pattern_struct.ref_candidates
        self.assertGreater(len(pattern_list), 0)
        self.assertEqual(pattern_list[0], 0)
        self.assertEqual(pattern_list[1], 1)
        self.assertEqual(pattern_list[2], 2)
        self.assertEqual(pattern_list[3], 3)
        self.assertEqual(pattern_list[4], 4)
        self.assertEqual(pattern_list[5], 5)

        pattern_struct = self.pyPPMb._construct_pattern_and_shift_rot_matrix(
            self.source_obj_array[[2, 4, 8, 16, 32, 64], :3], 6,
            np.cos(np.radians(60. / 3600.)), np.cos(np.radians(1.0)) ** 2,
            np.radians(5./3600.))
        pattern_list = pattern_struct.ref_candidates
        self.assertEqual(pattern_list[0], 2)
        self.assertEqual(pattern_list[1], 4)
        self.assertEqual(pattern_list[2], 8)
        self.assertEqual(pattern_list[3], 16)
        self.assertEqual(pattern_list[4], 32)
        self.assertEqual(pattern_list[5], 64)

    def testMatchPerfect(self):
        """ Input objects that have no shift or rotation to the matcher
        and test that we return a match.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)
        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array))
        self.assertTrue(
            np.all(match_struct.distances_rad < 0.01/3600.0 * __deg_to_rad__))

    def testOptimisticMatch(self):
        """ Test the optimistic mode of the pattern matcher. That is
        the algorithm with the early exit strategy as described in
        Tabur 2007.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array, n_check=9, n_match=6,
            n_agree=1, max_n_patterns=100, max_shift=60., max_rotation=6.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)
        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array))
        self.assertTrue(
            np.all(match_struct.distances_rad < 0.01/3600.0 * __deg_to_rad__))

    def testMatchSkip(self):
        """ Test the ability to skip specified patterns in the matching
        process.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=np.array([0]))
        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array))
        self.assertTrue(
            np.all(match_struct.distances_rad < 0.01/3600.0 * __deg_to_rad__))

    def testMatchMoreSources(self):
        """ Test the case where we have more sources than references
        but no rotation or shift.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:500, :3],
            log=self.log)

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60.0, max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)
        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array[:500]))
        self.assertTrue(
            np.all(match_struct.distances_rad < 0.01/3600.0 * __deg_to_rad__))

    def testMatchMoreReferences(self):
        """ Test the case where we have more references than sources
        but no rotation or shift.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array[:500], n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=1.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)
        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array[:500]))
        self.assertTrue(
            np.all(match_struct.distances_rad < 0.01/3600.0 * __deg_to_rad__))

    def testShift(self):
        """ Test the matcher when a shift is applied to the data.

        We say shift here as while we are rotating the unit-sphere in 3D, on
        our 'focal plane' this will appear as a shift.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)
        theta = np.radians(45.0 / 3600.)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        theta_rotation = self.pyPPMb._create_spherical_rotation_matrix(
            np.array([0, 0, 1]), cos_theta, sin_theta)

        self.source_obj_array[:, :3] = np.dot(
            theta_rotation,
            self.source_obj_array[:, :3].transpose()).transpose()

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60, max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)

        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array))
        self.assertTrue(
            np.all(match_struct.distances_rad < 0.01/3600.0 * __deg_to_rad__))

    def testRotation(self):
        """ Test the matcher for when a roation is applied to the data.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)
        phi = 2.5*__deg_to_rad__
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        phi_rotation = self.pyPPMb._create_spherical_rotation_matrix(
            np.array([1, 0, 0]), cos_phi, sin_phi)

        self.source_obj_array[:, :3] = np.dot(
            phi_rotation, self.source_obj_array[:, :3].transpose()).transpose()

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60, max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)

        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array))
        self.assertTrue(
            np.all(match_struct.distances_rad < 0.01/3600.0 * __deg_to_rad__))

    def testShiftRotation(self):
        """ Test both a shift and rotation being applied to the data.
        """
        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)
        theta = np.radians(45.0 / 3600.)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        theta_rotation = self.pyPPMb._create_spherical_rotation_matrix(
            np.array([0, 0, 1]), cos_theta, sin_theta)

        phi = 2.5 * __deg_to_rad__
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        phi_rotation = self.pyPPMb._create_spherical_rotation_matrix(
            np.array([1, 0, 0]), cos_phi, sin_phi)

        shift_rot_matrix = np.dot(theta_rotation, phi_rotation)

        self.source_obj_array[:, :3] = np.dot(
            shift_rot_matrix,
            self.source_obj_array[:, :3].transpose()).transpose()

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)

        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array))
        self.assertTrue(
            np.all(match_struct.distances_rad < 0.01/3600.0 * __deg_to_rad__))

    def testLinearDistortion(self):
        """ Create a simple linear distortion and test that the correct
        references are still matched.
        """

        self.pyPPMb = PessimisticPatternMatcherB(
            reference_array=self.reference_obj_array[:, :3],
            log=self.log)

        max_z = np.cos(np.pi/2 + 0.5 * __deg_to_rad__)
        min_z = np.cos(np.pi/2 - 0.5 * __deg_to_rad__)
        # The max shift in position we add to the position will be 25
        # arcseconds.
        max_distort = 25.0 / 3600. * __deg_to_rad__
        self.source_obj_array[:, 2] = (
            self.source_obj_array[:, 2]
            - max_distort * (self.source_obj_array[:, 2] - min_z)
            / (max_z - min_z))
        # Renomalize the 3 vectors to be unit length.
        distorted_dists = np.sqrt(self.source_obj_array[:, 0] ** 2
                                  + self.source_obj_array[:, 1] ** 2
                                  + self.source_obj_array[:, 2] ** 2)
        self.source_obj_array[:, 0] /= distorted_dists
        self.source_obj_array[:, 1] /= distorted_dists
        self.source_obj_array[:, 2] /= distorted_dists

        match_struct = self.pyPPMb.match(
            source_array=self.source_obj_array, n_check=9, n_match=6,
            n_agree=2, max_n_patterns=100, max_shift=60., max_rotation=5.0,
            max_dist=5., min_matches=30, pattern_skip_array=None)

        self.assertEqual(len(match_struct.match_ids),
                         len(self.reference_obj_array))
        self.assertTrue(
            np.all(match_struct.distances_rad < 10 / 3600.0 * __deg_to_rad__))

    def testNoReferenceSources(self):
        """Check that we get a helpful error when no reference objects are
        supplied.
        """
        with self.assertRaisesRegex(ValueError, "No reference objects supplied"):
            PessimisticPatternMatcherB(np.ndarray((0, 3)), self.log)


if __name__ == '__main__':
    unittest.main()
