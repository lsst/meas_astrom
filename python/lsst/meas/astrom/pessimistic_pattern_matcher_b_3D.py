
__all__ = ["PessimisticPatternMatcherB"]

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.stats import sigmaclip

import lsst.pipe.base as pipeBase

from .pessimisticPatternMatcherUtils import construct_pattern_and_shift_rot_matrix


def _rotation_matrix_chi_sq(flattened_rot_matrix,
                            pattern_a,
                            pattern_b,
                            max_dist_rad):
    """Compute the squared differences for least squares fitting.

    Given a flattened rotation matrix, one N point pattern and another N point
    pattern to transform into to, compute the squared differences between the
    points in the two patterns after the rotation.

    Parameters
    ----------
    flattened_rot_matrix : `numpy.ndarray`, (9, )
        A flattened array representing a 3x3 rotation matrix. The array is
        flattened to comply with the API of scipy.optimize.least_squares.
        Flattened elements are [[0, 0], [0, 1], [0, 2], [1, 0]...]
    pattern_a : `numpy.ndarray`, (N, 3)
        A array containing N, 3 vectors representing the objects we would like
        to transform into the frame of pattern_b.
    pattern_b : `numpy.ndarray`, (N, 3)
        A array containing N, 3 vectors representing the reference frame we
        would like to transform pattern_a into.
    max_dist_rad : `float`
        The maximum distance allowed from the pattern matching. This value is
        used as the standard error for the resultant chi values.

    Returns
    -------
    noralized_diff : `numpy.ndarray`, (9,)
        Array of differences between the vectors representing of the source
        pattern rotated into the reference frame and the converse. This is
        used to minimize in a least squares fitter.
    """
    # Unflatten the rotation matrix
    rot_matrix = flattened_rot_matrix.reshape((3, 3))
    # Compare the rotated source pattern to the references.
    rot_pattern_a = np.dot(rot_matrix, pattern_a.transpose()).transpose()
    diff_pattern_a_to_b = rot_pattern_a - pattern_b
    # Return the flattened differences and length tolerances for use in a least
    # squares fitter.
    return diff_pattern_a_to_b.flatten() / max_dist_rad


class PessimisticPatternMatcherB:
    """Class implementing a pessimistic version of Optimistic Pattern Matcher
    B (OPMb) from `Tabur (2007)`_, as described in `DMTN-031`_

    Parameters
    ----------
    reference_array : `numpy.ndarray`, (N, 3)
        spherical points x, y, z of to use as reference objects for
        pattern matching.
    log : `lsst.log.Log` or `logging.Logger`
        Logger for outputting debug info.

    Notes
    -----
    The class loads and stores the reference object
    in a convenient data structure for matching any set of source objects that
    are assumed to contain each other. The pessimistic nature of the algorithm
    comes from requiring that it discovers at least two patterns that agree on
    the correct shift and rotation for matching before exiting. The original
    behavior of OPMb can be recovered simply. Patterns matched between the
    input datasets are n-spoked pinwheels created from n+1 points. Refer to
    `DMTN-031`_ for more details.

    .. _Tabur (2007): https://doi.org/10.1071/AS07028
    .. _DMTN-031: https://dmtn-031.lsst.io/
    """

    def __init__(self, reference_array, log):
        self._reference_array = reference_array
        self._n_reference = len(self._reference_array)
        self.log = log

        if self._n_reference > 0:
            self._build_distances_and_angles()
        else:
            raise ValueError("No reference objects supplied")

    def _build_distances_and_angles(self):
        """Create the data structures we will use to search for our pattern
        match in.

        Throughout this function and the rest of the class we use id to
        reference the position in the input reference catalog and index to
        'index' into the arrays sorted on distance.
        """
        # Create empty arrays to store our pair information per
        # reference object.
        self._dist_array = np.empty(
            int(self._n_reference * (self._n_reference - 1) / 2),
            dtype="float32")
        self._id_array = np.empty(
            (int(self._n_reference * (self._n_reference - 1) / 2), 2),
            dtype="uint16")

        startIdx = 0
        # Loop over reference objects storing pair distances and ids.
        for ref_id, ref_obj in enumerate(self._reference_array):
            # Set the ending slicing index to the correct length for the
            # pairs we are creating.
            endIdx = startIdx + self._n_reference - 1 - ref_id

            # Reserve and fill the ids of each reference object pair.
            # 16 bit is safe for the id array as the catalog input from
            # MatchPessimisticB is limited to a max length of 2 ** 16.
            self._id_array[startIdx:endIdx, 0] = ref_id
            self._id_array[startIdx:endIdx, 1] = np.arange(ref_id + 1,
                                                           self._n_reference,
                                                           dtype="uint16")

            # Compute the vector deltas for each pair of reference objects.
            # Compute and store the distances.
            self._dist_array[startIdx:endIdx] = np.sqrt(
                ((self._reference_array[ref_id + 1:, :]
                  - ref_obj) ** 2).sum(axis=1))
            # Set startIdx of the slice to the end of the previous slice.
            startIdx = endIdx

        # Sort each array on the pair distances for the initial
        # optimistic pattern matcher lookup.
        sorted_dist_args = self._dist_array.argsort()
        self._dist_array = self._dist_array[sorted_dist_args]
        self._id_array = self._id_array[sorted_dist_args]

    def match(self, source_array, n_check, n_match, n_agree,
              max_n_patterns, max_shift, max_rotation, max_dist,
              min_matches, pattern_skip_array=None):
        """Match a given source catalog into the loaded reference catalog.

        Given array of points on the unit sphere and tolerances, we
        attempt to match a pinwheel like pattern between these input sources
        and the reference objects this class was created with. This pattern
        informs of the shift and rotation needed to align the input source
        objects into the frame of the references.

        Parameters
        ----------
        source_array : `numpy.ndarray`, (N, 3)
            An array of spherical x,y,z coordinates and a magnitude in units
            of objects having a lower value for sorting. The array should be
            of shape (N, 4).
        n_check  : `int`
            Number of sources to create a pattern from. Not all objects may be
            checked if n_match criteria is before looping through all n_check
            objects.
        n_match : `int`
            Number of objects to use in constructing a pattern to match.
        n_agree : `int`
            Number of found patterns that must agree on their shift and
            rotation before exiting. Set this value to 1 to recover the
            expected behavior of Optimistic Pattern Matcher B.
        max_n_patters : `int`
            Number of patterns to create from the input source objects to
            attempt to match into the reference objects.
        max_shift : `float`
            Maximum allowed shift to match patterns in arcseconds.
        max_rotation : `float`
            Maximum allowed rotation between patterns in degrees.
        max_dist : `float`
            Maximum distance in arcseconds allowed between candidate spokes in
            the source and reference objects. Also sets that maximum distance
            in the intermediate verify, pattern shift/rotation agreement, and
            final verify steps.
        pattern_skip_array : `int`
            Patterns we would like to skip. This could be due to the pattern
            being matched on a previous iteration that we now consider invalid.
            This assumes the ordering of the source objects is the same
            between different runs of the matcher which, assuming no object
            has been inserted or the magnitudes have changed, it should be.

        Returns
        -------
        output_struct : `lsst.pipe.base.Struct`
            Result struct with components

            - ``matches`` : (N, 2) array of matched ids for pairs. Empty list if no
              match found (`numpy.ndarray`, (N, 2) or `list`)
            - ``distances_rad`` : Radian distances between the matched objects.
              Empty list if no match found (`numpy.ndarray`, (N,))
            - ``pattern_idx``: Index of matched pattern. None if no match found
              (`int`).
            - ``shift`` : Magnitude for the shift between the source and reference
              objects in arcseconds. None if no match found (`float`).
        """
        # Given our input source_array we sort on magnitude.
        sorted_source_array = source_array[source_array[:, -1].argsort(), :3]
        n_source = len(sorted_source_array)

        # Initialize output struct.
        output_match_struct = pipeBase.Struct(
            match_ids=[],
            distances_rad=[],
            pattern_idx=None,
            shift=None,
            max_dist_rad=None,)

        if n_source <= 0:
            self.log.warning("Source object array is empty. Unable to match. Exiting matcher.")
            return None

        # To test if the shifts and rotations we find agree with each other,
        # we first create two test points situated at the top and bottom of
        # where the z axis on the sphere bisects the source catalog.
        test_vectors = self._compute_test_vectors(source_array[:, :3])

        # We now create an empty list of our resultant rotated vectors to
        # compare the different rotations we find.
        rot_vect_list = []

        # Convert the tolerances to values we will use in the code.
        max_cos_shift = np.cos(np.radians(max_shift / 3600.))
        max_cos_rot_sq = np.cos(np.radians(max_rotation)) ** 2
        max_dist_rad = np.radians(max_dist / 3600.)

        # Loop through the sources from brightest to faintest, grabbing a
        # chunk of n_check each time.
        for pattern_idx in range(np.min((max_n_patterns,
                                         n_source - n_match))):

            # If this pattern is one that we matched on the past but we
            # now want to skip, we do so here.
            if pattern_skip_array is not None and \
               np.any(pattern_skip_array == pattern_idx):
                self.log.debug(
                    "Skipping previously matched bad pattern %i...",
                    pattern_idx)
                continue
            # Grab the sources to attempt to create this pattern.
            pattern = sorted_source_array[
                pattern_idx: np.min((pattern_idx + n_check, n_source)), :3]

            # Construct a pattern given the number of points defining the
            # pattern complexity. This is the start of the primary tests to
            # match our source pattern into the reference objects.
            construct_return_struct = \
                self._construct_pattern_and_shift_rot_matrix(
                    pattern, n_match, max_cos_shift, max_cos_rot_sq,
                    max_dist_rad)

            # Our struct is None if we could not match the pattern.
            if construct_return_struct.ref_candidates is None or \
               construct_return_struct.shift_rot_matrix is None or \
               construct_return_struct.cos_shift is None or \
               construct_return_struct.sin_rot is None:
                continue

            # Grab the output data from the Struct object.
            ref_candidates = construct_return_struct.ref_candidates
            shift_rot_matrix = construct_return_struct.shift_rot_matrix
            cos_shift = construct_return_struct.cos_shift
            sin_rot = construct_return_struct.sin_rot

            # If we didn't match enough candidates we continue to the next
            # pattern.
            if len(ref_candidates) < n_match:
                continue

            # Now that we know our pattern and shift/rotation are valid we
            # store the the rotated versions of our test points for later
            # use.
            tmp_rot_vect_list = []
            for test_vect in test_vectors:
                tmp_rot_vect_list.append(np.dot(shift_rot_matrix, test_vect))
            # Since our test point are in the source frame, we can test if
            # their lengths are mostly preserved under the transform.
            if not self._test_pattern_lengths(np.array(tmp_rot_vect_list),
                                              max_dist_rad):
                continue

            tmp_rot_vect_list.append(pattern_idx)
            rot_vect_list.append(tmp_rot_vect_list)

            # Test if we have enough rotations, which agree, or if we
            # are in optimistic mode.
            if self._test_rotation_agreement(rot_vect_list, max_dist_rad) < \
               n_agree - 1:
                continue

            # Run the final verify step.
            match_struct = self._final_verify(source_array[:, :3],
                                              shift_rot_matrix,
                                              max_dist_rad,
                                              min_matches)
            if match_struct.match_ids is None or \
               match_struct.distances_rad is None or \
               match_struct.max_dist_rad is None:
                continue

            # Convert the observed shift to arcseconds
            shift = np.degrees(np.arccos(cos_shift)) * 3600.
            # Print information to the logger.
            self.log.debug("Succeeded after %i patterns.", pattern_idx)
            self.log.debug("\tShift %.4f arcsec", shift)
            self.log.debug("\tRotation: %.4f deg",
                           np.degrees(np.arcsin(sin_rot)))

            # Fill the struct and return.
            output_match_struct.match_ids = \
                match_struct.match_ids
            output_match_struct.distances_rad = \
                match_struct.distances_rad
            output_match_struct.pattern_idx = pattern_idx
            output_match_struct.shift = shift
            output_match_struct.max_dist_rad = match_struct.max_dist_rad
            return output_match_struct

        self.log.debug("Failed after %i patterns.", pattern_idx + 1)
        return output_match_struct

    def _compute_test_vectors(self, source_array):
        """Compute spherical 3 vectors at the edges of the x, y, z extent
        of the input source catalog.

        Parameters
        ----------
        source_array : `numpy.ndarray`, (N, 3)
            array of 3 vectors representing positions on the unit
            sphere.

        Returns
        -------
        test_vectors : `numpy.ndarray`, (N, 3)
            Array of vectors representing the maximum extents in x, y, z
            of the input source array. These are used with the rotations
            the code finds to test for agreement from different patterns
            when the code is running in pessimistic mode.
        """
        # Get the center of source_array.
        if np.any(np.logical_not(np.isfinite(source_array))):
            self.log.warning("Input source objects contain non-finite values. "
                             "This could end badly.")
        center_vect = np.nanmean(source_array, axis=0)

        # So that our rotation test works over the full sky we compute
        # the max extent in each Cartesian direction x,y,z.
        xbtm_vect = np.array([np.min(source_array[:, 0]), center_vect[1],
                              center_vect[2]], dtype=np.float64)
        xtop_vect = np.array([np.max(source_array[:, 0]), center_vect[1],
                              center_vect[2]], dtype=np.float64)
        xbtm_vect /= np.sqrt(np.dot(xbtm_vect, xbtm_vect))
        xtop_vect /= np.sqrt(np.dot(xtop_vect, xtop_vect))

        ybtm_vect = np.array([center_vect[0], np.min(source_array[:, 1]),
                              center_vect[2]], dtype=np.float64)
        ytop_vect = np.array([center_vect[0], np.max(source_array[:, 1]),
                              center_vect[2]], dtype=np.float64)
        ybtm_vect /= np.sqrt(np.dot(ybtm_vect, ybtm_vect))
        ytop_vect /= np.sqrt(np.dot(ytop_vect, ytop_vect))

        zbtm_vect = np.array([center_vect[0], center_vect[1],
                              np.min(source_array[:, 2])], dtype=np.float64)
        ztop_vect = np.array([center_vect[0], center_vect[1],
                              np.max(source_array[:, 2])], dtype=np.float64)
        zbtm_vect /= np.sqrt(np.dot(zbtm_vect, zbtm_vect))
        ztop_vect /= np.sqrt(np.dot(ztop_vect, ztop_vect))

        # Return our list of vectors for later rotation testing.
        return np.array([xbtm_vect, xtop_vect, ybtm_vect, ytop_vect,
                         zbtm_vect, ztop_vect])

    def _construct_pattern_and_shift_rot_matrix(self, src_pattern_array,
                                                n_match, max_cos_theta_shift,
                                                max_cos_rot_sq, max_dist_rad):
        """Test an input source pattern against the reference catalog.
        Returns the candidate matched patterns and their
        implied rotation matrices or None.

        Parameters
        ----------
        src_pattern_array : `numpy.ndarray`, (N, 3)
            Sub selection of source 3 vectors to create a pattern from
        n_match : `int`
            Number of points to attempt to create a pattern from. Must be
            >= len(src_pattern_array)
        max_cos_theta_shift : `float`
            Maximum shift allowed between two patterns' centers.
        max_cos_rot_sq : `float`
            Maximum rotation between two patterns that have been shifted
            to have their centers on top of each other.
        max_dist_rad : `float`
            Maximum delta distance allowed between the source and reference
            pair distances to consider the reference pair a candidate for
            the source pair. Also sets the tolerance between the opening
            angles of the spokes when compared to the reference.

        Return
        -------
        output_matched_pattern : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``ref_candidates`` : ids of the matched pattern in the internal
              reference_array object (`list` of `int`).
            - ``src_candidates`` : Pattern ids of the sources matched
              (`list` of `int`).
            - ``shift_rot_matrix_src_to_ref`` : 3x3 matrix specifying the full
              shift and rotation between the reference and source objects.
              Rotates source into reference frame. `None` if match is not
              found. (`numpy.ndarray`, (3, 3))
            - ``shift_rot_matrix_ref_to_src`` : 3x3 matrix specifying the full
              shift and rotation of the reference and source objects. Rotates
              reference into source frame. None if match is not found
              (`numpy.ndarray`, (3, 3)).
            - ``cos_shift`` : Magnitude of the shift found between the two
              patten centers. `None` if match is not found (`float`).
            - ``sin_rot`` : float value of the rotation to align the already
              shifted source pattern to the reference pattern. `None` if no match
              found (`float`).
        """
        # Create the delta vectors and distances we will need to assemble the
        # spokes of the pattern.
        src_delta_array = np.empty((len(src_pattern_array) - 1, 3))
        src_delta_array[:, 0] = (src_pattern_array[1:, 0]
                                 - src_pattern_array[0, 0])
        src_delta_array[:, 1] = (src_pattern_array[1:, 1]
                                 - src_pattern_array[0, 1])
        src_delta_array[:, 2] = (src_pattern_array[1:, 2]
                                 - src_pattern_array[0, 2])
        src_dist_array = np.sqrt(src_delta_array[:, 0]**2
                                 + src_delta_array[:, 1]**2
                                 + src_delta_array[:, 2]**2)

        pattern_result = construct_pattern_and_shift_rot_matrix(
            src_pattern_array, src_delta_array, src_dist_array,
            self._dist_array, self._id_array, self._reference_array, n_match,
            max_cos_theta_shift, max_cos_rot_sq, max_dist_rad)

        if pattern_result.success:
            candidate_array = np.array(pattern_result.candidate_pairs)
            fit_shift_rot_matrix = self._intermediate_verify(
                src_pattern_array[candidate_array[:, 1]],
                self._reference_array[candidate_array[:, 0]],
                pattern_result.shift_rot_matrix, max_dist_rad)
            return pipeBase.Struct(ref_candidates=candidate_array[:, 0].tolist(),
                                   src_candidates=candidate_array[:, 1].tolist(),
                                   shift_rot_matrix=fit_shift_rot_matrix,
                                   cos_shift=pattern_result.cos_shift,
                                   sin_rot=pattern_result.sin_rot)
        return pipeBase.Struct(ref_candidates=[],
                               src_candidates=[],
                               shift_rot_matrix=None,
                               cos_shift=None,
                               sin_rot=None)

    def _intermediate_verify(self, src_pattern, ref_pattern, shift_rot_matrix,
                             max_dist_rad):
        """ Perform an intermediate verify step.
        Rotate the matches references into the source frame and test their
        distances against tolerance. Only return true if all points are within
        tolerance.

        Parameters
        ----------
        src_pattern : `numpy.ndarray`, (N,3)
            Array of 3 vectors representing the points that make up our source
            pinwheel pattern.
        ref_pattern : `numpy.ndarray`, (N,3)
            Array of 3 vectors representing our candidate reference pinwheel
            pattern.
        shift_rot_matrix : `numpy.ndarray`, (3,3)
            3x3 rotation matrix that takes the source objects and rotates them
            onto the frame of the reference objects
        max_dist_rad : `float`
            Maximum distance allowed to consider two objects the same.

        Returns
        -------
        fit_shift_rot_matrix : `numpy.ndarray`, (3,3)
           Fitted shift/rotation matrix if all of the points in our
           source pattern are within max_dist_rad of their matched reference
           objects. Returns None if this criteria is not satisfied.
        """
        if len(src_pattern) != len(ref_pattern):
            raise ValueError(
                "Source pattern length does not match ref pattern.\n"
                "\t source pattern len=%i, reference pattern len=%i" %
                (len(src_pattern), len(ref_pattern)))

        if self._intermediate_verify_comparison(
                src_pattern, ref_pattern, shift_rot_matrix, max_dist_rad):
            # Now that we know our initial shift and rot matrix is valid we
            # want to fit the implied transform using all points from
            # our pattern. This is a more robust rotation matrix as our
            # initial matrix only used the first 2 points from the source
            # pattern to estimate the shift and rotation. The matrix we fit
            # are allowed to be non unitary but need to preserve the length of
            # the rotated vectors to within the match tolerance.
            fit_shift_rot_matrix = least_squares(
                _rotation_matrix_chi_sq,
                x0=shift_rot_matrix.flatten(),
                args=(src_pattern, ref_pattern, max_dist_rad)
            ).x.reshape((3, 3))
            # Do another verify in case the fits have wondered off.
            if self._intermediate_verify_comparison(
                    src_pattern, ref_pattern, fit_shift_rot_matrix,
                    max_dist_rad):
                return fit_shift_rot_matrix

        return None

    def _create_spherical_rotation_matrix(self, rot_axis, cos_rotation,
                                          sin_rotion):
        """Construct a generalized 3D rotation matrix about a given
        axis.

        Parameters
        ----------
        rot_axis : `numpy.ndarray`, (3,)
            3 vector defining the axis to rotate about.
        cos_rotation : `float`
            cosine of the rotation angle.
        sin_rotation : `float`
            sine of the rotation angle.

        Return
        ------
        shift_matrix : `numpy.ndarray`, (3, 3)
            3x3 spherical, rotation matrix.
        """
        rot_cross_matrix = np.array(
            [[0., -rot_axis[2], rot_axis[1]],
             [rot_axis[2], 0., -rot_axis[0]],
             [-rot_axis[1], rot_axis[0], 0.]], dtype=np.float64)
        shift_matrix = (cos_rotation*np.identity(3)
                        + sin_rotion*rot_cross_matrix
                        + (1. - cos_rotation)*np.outer(rot_axis, rot_axis))

        return shift_matrix

    def _intermediate_verify_comparison(self, pattern_a, pattern_b,
                                        shift_rot_matrix, max_dist_rad):
        """Test the input rotation matrix against one input pattern and
        a second one.

        If every point in the pattern after rotation is within a distance of
        max_dist_rad to its candidate point in the other pattern, we return
        True.

        Parameters
        ----------
        pattern_a :  `numpy.ndarray`, (N,3)
            Array of 3 vectors representing the points that make up our source
            pinwheel pattern.
        pattern_b :  `numpy.ndarray`, (N,3)
            Array of 3 vectors representing our candidate reference pinwheel
            pattern.
        shift_rot_matrix :  `numpy.ndarray`, (3,3)
            3x3 rotation matrix that takes the source objects and rotates them
            onto the frame of the reference objects
        max_dist_rad : `float`
            Maximum distance allowed to consider two objects the same.


        Returns
        -------
        bool
            True if all rotated source points are within max_dist_rad of
            the candidate references matches.
        """
        shifted_pattern_a = np.dot(shift_rot_matrix,
                                   pattern_a.transpose()).transpose()
        tmp_delta_array = shifted_pattern_a - pattern_b
        tmp_dist_array = (tmp_delta_array[:, 0] ** 2
                          + tmp_delta_array[:, 1] ** 2
                          + tmp_delta_array[:, 2] ** 2)
        return np.all(tmp_dist_array < max_dist_rad ** 2)

    def _test_pattern_lengths(self, test_pattern, max_dist_rad):
        """ Test that the all vectors in a pattern are unit length within
        tolerance.

        This is useful for assuring the non unitary transforms do not contain
        too much distortion.

        Parameters
        ----------
        test_pattern : `numpy.ndarray`, (N, 3)
            Test vectors at the maximum and minimum x, y, z extents.
        max_dist_rad : `float`
            maximum distance in radians to consider two points "agreeing" on
            a rotation

        Returns
        -------
        test : `bool`
            Length tests pass.
        """
        dists = (test_pattern[:, 0] ** 2
                 + test_pattern[:, 1] ** 2
                 + test_pattern[:, 2] ** 2)
        return np.all(
            np.logical_and((1 - max_dist_rad) ** 2 < dists,
                           dists < (1 + max_dist_rad) ** 2))

    def _test_rotation_agreement(self, rot_vects, max_dist_rad):
        """ Test this rotation against the previous N found and return
        the number that a agree within tolerance to where our test
        points are.

        Parameters
        ----------
        rot_vects : `numpy.ndarray`, (N, 3)
            Arrays of rotated 3 vectors representing the maximum x, y,
            z extent on the unit sphere of the input source objects rotated by
            the candidate rotations into the reference frame.
        max_dist_rad : `float`
            maximum distance in radians to consider two points "agreeing" on
            a rotation

        Returns
        -------
        tot_consent : `int`
            Number of candidate rotations that agree for all of the rotated
            test 3 vectors.
        """
        self.log.debug("Comparing pattern %i to previous %i rotations...",
                       rot_vects[-1][-1], len(rot_vects) - 1)

        tot_consent = 0
        for rot_idx in range(max((len(rot_vects) - 1), 0)):
            tmp_dist_list = []
            for vect_idx in range(len(rot_vects[rot_idx]) - 1):
                tmp_delta_vect = (rot_vects[rot_idx][vect_idx]
                                  - rot_vects[-1][vect_idx])
                tmp_dist_list.append(
                    np.dot(tmp_delta_vect, tmp_delta_vect))
            if np.all(np.array(tmp_dist_list) < max_dist_rad ** 2):
                tot_consent += 1
        return tot_consent

    def _final_verify(self,
                      source_array,
                      shift_rot_matrix,
                      max_dist_rad,
                      min_matches):
        """Match the all sources into the reference catalog using the shift/rot
        matrix.

        After the initial shift/rot matrix is found, we refit the shift/rot
        matrix using the matches the initial matrix produces to find a more
        stable solution.

        Parameters
        ----------
        source_array : `numpy.ndarray` (N, 3)
            3-vector positions on the unit-sphere representing the sources to
            match
        shift_rot_matrix : `numpy.ndarray` (3, 3)
            Rotation matrix representing inferred shift/rotation of the
            sources onto the reference catalog. Matrix need not be unitary.
        max_dist_rad : `float`
            Maximum distance allowed for a match.
        min_matches : `int`
            Minimum number of matched objects required to consider the
            match good.

        Returns
        -------
        output_struct : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``match_ids`` : Pairs of indexes into the source and reference
              data respectively defining a match (`numpy.ndarray`, (N, 2)).
            - ``distances_rad`` : distances to between the matched objects in
              the shift/rotated frame. (`numpy.ndarray`, (N,)).
            - ``max_dist_rad`` : Value of the max matched distance. Either
              returning the input value of the 2 sigma clipped value of the
              shift/rotated distances. (`float`).
        """
        output_struct = pipeBase.Struct(
            match_ids=None,
            distances_rad=None,
            max_dist_rad=None,
        )

        # Perform an iterative final verify.
        match_sources_struct = self._match_sources(source_array,
                                                   shift_rot_matrix)
        cut_ids = match_sources_struct.match_ids[
            match_sources_struct.distances_rad < max_dist_rad]

        n_matched = len(cut_ids)
        clipped_struct = self._clip_distances(
            match_sources_struct.distances_rad)
        n_matched_clipped = clipped_struct.n_matched_clipped

        if n_matched < min_matches or n_matched_clipped < min_matches:
            return output_struct

        # The shift/rotation matrix returned by
        # ``_construct_pattern_and_shift_rot_matrix``, above, was
        # based on only six points. Here, we refine that result by
        # using all of the good matches from the “final verification”
        # step, above. This will produce a more consistent result.
        fit_shift_rot_matrix = least_squares(
            _rotation_matrix_chi_sq,
            x0=shift_rot_matrix.flatten(),
            args=(source_array[cut_ids[:, 0], :3],
                  self._reference_array[cut_ids[:, 1], :3],
                  max_dist_rad)
        ).x.reshape((3, 3))

        # Redo the matching using the newly fit shift/rotation matrix.
        match_sources_struct = self._match_sources(
            source_array, fit_shift_rot_matrix)

        # Double check the match distances to make sure enough matches
        # survive still. We'll just overwrite the previously used variables.
        n_matched = np.sum(
            match_sources_struct.distances_rad < max_dist_rad)
        clipped_struct = self._clip_distances(
            match_sources_struct.distances_rad)
        n_matched_clipped = clipped_struct.n_matched_clipped
        clipped_max_dist = clipped_struct.clipped_max_dist

        if n_matched < min_matches or n_matched_clipped < min_matches:
            return output_struct

        # Store our matches in the output struct. Decide between the clipped
        # distance and the input max dist based on which is smaller.
        output_struct.match_ids = match_sources_struct.match_ids
        output_struct.distances_rad = match_sources_struct.distances_rad
        if clipped_max_dist < max_dist_rad:
            output_struct.max_dist_rad = clipped_max_dist
        else:
            output_struct.max_dist_rad = max_dist_rad

        return output_struct

    def _clip_distances(self, distances_rad):
        """Compute a clipped max distance and calculate the number of pairs
        that pass the clipped dist.

        Parameters
        ----------
        distances_rad : `numpy.ndarray`, (N,)
            Distances between pairs.

        Returns
        -------
        output_struct : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``n_matched_clipped`` : Number of pairs that survive the
              clipping on distance. (`float`)
            - ``clipped_max_dist`` : Maximum distance after clipping.
              (`float`).
        """
        clipped_dists, _, clipped_max_dist = sigmaclip(
            distances_rad,
            low=100,
            high=2)
        # Check clipped distances. The minimum value here
        # prevents over convergence on perfect test data.
        if clipped_max_dist < 1e-16:
            clipped_max_dist = 1e-16
            n_matched_clipped = np.sum(distances_rad < clipped_max_dist)
        else:
            n_matched_clipped = len(clipped_dists)

        return pipeBase.Struct(n_matched_clipped=n_matched_clipped,
                               clipped_max_dist=clipped_max_dist)

    def _match_sources(self,
                       source_array,
                       shift_rot_matrix):
        """ Shift both the reference and source catalog to the the respective
        frames and find their nearest neighbor using a kdTree.

        Removes all matches who do not agree when either the reference or
        source catalog is rotated. Cuts on a maximum distance are left to an
        external function.

        Parameters
        ----------
        source_array : `numpy.ndarray`, (N, 3)
            array of 3 vectors representing the source objects we are trying
            to match into the source catalog.
        shift_rot_matrix : `numpy.ndarray`, (3, 3)
            3x3 rotation matrix that performs the spherical rotation from the
            source frame into the reference frame.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``matches`` : array of integer ids into the source and
              reference arrays. Matches are only returned for those that
              satisfy the distance and handshake criteria
              (`numpy.ndarray`, (N, 2)).
            - ``distances`` : Distances between each match in radians after
              the shift and rotation is applied (`numpy.ndarray`, (N)).
        """
        shifted_references = np.dot(
            np.linalg.inv(shift_rot_matrix),
            self._reference_array.transpose()).transpose()
        shifted_sources = np.dot(
            shift_rot_matrix,
            source_array.transpose()).transpose()

        ref_matches = np.empty((len(shifted_references), 2),
                               dtype="uint16")
        src_matches = np.empty((len(shifted_sources), 2),
                               dtype="uint16")

        ref_matches[:, 1] = np.arange(len(shifted_references),
                                      dtype="uint16")
        src_matches[:, 0] = np.arange(len(shifted_sources),
                                      dtype="uint16")

        ref_kdtree = cKDTree(self._reference_array)
        src_kdtree = cKDTree(source_array)

        ref_to_src_dist, tmp_ref_to_src_idx = \
            src_kdtree.query(shifted_references)
        src_to_ref_dist, tmp_src_to_ref_idx = \
            ref_kdtree.query(shifted_sources)

        ref_matches[:, 0] = tmp_ref_to_src_idx
        src_matches[:, 1] = tmp_src_to_ref_idx

        handshake_mask = self._handshake_match(src_matches, ref_matches)
        return pipeBase.Struct(
            match_ids=src_matches[handshake_mask],
            distances_rad=src_to_ref_dist[handshake_mask],)

    def _handshake_match(self, matches_src, matches_ref):
        """Return only those matches where both the source
        and reference objects agree they they are each others'
        nearest neighbor.

        Parameters
        ----------
        matches_src : `numpy.ndarray`, (N, 2)
            int array of nearest neighbor matches between shifted and
            rotated reference objects matched into the sources.
        matches_ref : `numpy.ndarray`, (N, 2)
            int array of nearest neighbor matches between shifted and
            rotated source objects matched into the references.

        Return
        ------
        handshake_mask_array : `numpy.ndarray`, (N,)
            Array positions where the two match catalogs agree.
        """
        handshake_mask_array = np.zeros(len(matches_src), dtype=bool)

        for src_match_idx, match in enumerate(matches_src):
            ref_match_idx = np.searchsorted(matches_ref[:, 1], match[1])
            if match[0] == matches_ref[ref_match_idx, 0]:
                handshake_mask_array[src_match_idx] = True
        return handshake_mask_array
