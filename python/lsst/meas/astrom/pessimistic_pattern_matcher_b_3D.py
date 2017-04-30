
from __future__ import division, print_function, absolute_import

from buildtins import xrange

import numpy as np
from scipy.spatial import cKDTree

import lsst.pipe.base as pipeBase


class PessimisticPatternMatcherB(object):
    """ Class implementing a pessimistic version of Optimsitic Pattern Matcher
    B (OPMb) from Tabur 2007. The class loads and stores the reference object
    in a convienent data structure for matching any set of source obejcts that
    are assumed to contain each other. The pessimistic nature of the algorithm
    comes from requiring that it discovers at least two patterns that agree on
    the correct shift and rotation for matching before exiting. The original
    behavior of OPMb can be recovered simply. Patterns matched between the
    input datasets are n-spoked pinwheels created from n+1 points. Refer to
    DMTN #031 for more details. http://github.com/lsst-dm/dmtn-031
    ----------------------------------------------------------------------------
    Attributes:
        reference_array : float array
            spherical points x, y, z of to use as reference objects for pattern
            matching.
        log : an lsst.log instance
        pair_id_array : int array
           Internal lookup table. Given an id in the reference array, return an
           array of the id pair that contains this object's id sorted on the
           distance to the pairs.
        pair_delta_array : float array
           Internal lookup table. Given an id in the reference array, return an
           array of the 3 vector deltas of all other pairs sorted on their
           distance.
        pair_dist_array : float array
           Internal lookup table. Given an id in the reference return an array
           of pair distances of all other pairs sorted on distance.
        dist_array : float array
            Array of all pairs of objects in the reference array sorted on
            distance.
        id_array : int array
            Array of id pairs that lookup into the reference array sorted
            on pair distance.
        delta_array : float array
           Array of 3 vector deltas for each pair in the reference array sorted
           on pair distance.
    """

    def __init__(self, reference_array, log):
        """
        Arguments
        ---------
        reference_array : float array
            Array of spherical points x, y, z to use as reference objects.
        log : lsst.log
            logger object for reporting warnings and failures.
        """
        self._reference_array = reference_array
        self._n_reference = len(self._reference_array)
        self.log = log

        self._build_distances_and_angles()

    def _build_distances_and_angles(self):
        """ Create the data structures we will use to search for our pattern match in.

        Throughout this function and the rest of the

        class we use id to reference the position in the input reference
        catalog and index to 'index' into the arrays sorted on distance.

        Parameters
        ----------
        self

        Returns
        -------
        None
        """

        # Initialize the arrays we will need for quick look up of pairs once
        # have a candidate spoke center.
        self._pair_id_array = np.empty(
            (self._n_reference, self._n_reference - 1),
            dtype=np.uint32)
        self._pair_delta_array = np.empty(
            (self._n_reference, self._n_reference - 1, 3),
            dtype=np.float64)
        self._pair_dist_array = np.empty(
            (self._n_reference, self._n_reference - 1),
            dtype=np.float64)

        # Create empty lists to temporarially store our pair information per
        # reference object. These will be concatentated into our final arrays.
        sub_id_array_list = []
        sub_delta_array_list = []
        sub_dist_array_list = []

        # Loop over reference objects and store pair distances, ids, and
        # 3 vector deltas.
        for ref_id, ref_obj in enumerate(self._reference_array):

            # Reserve and fill the ids of each reference object pair.
            sub_id_array = np.zeros((self._n_reference - 1 - ref_id, 2),
                                    dtype=np.uint32)
            sub_id_array[:, 0] = ref_id
            sub_id_array[:, 1] = np.arange(ref_id + 1, self._n_reference,
                                           dtype=np.uint32)

            # Compute the vector deltas for each pair of reference objects
            # and compute and store the distances.
            sub_delta_array = (self._reference_array[ref_id + 1:, :] -
                               ref_obj)
            sub_dist_array = np.sqrt(sub_delta_array[:, 0] ** 2 +
                                     sub_delta_array[:, 1] ** 2 +
                                     sub_delta_array[:, 2] ** 2)

            # Append to our arrays to the output lists for later concatenation.
            sub_id_array_list.append(sub_id_array)
            sub_delta_array_list.append(sub_delta_array)
            sub_dist_array_list.append(sub_dist_array)

            # Fill the pair look up arrays row wise and then column wise.
            self._pair_id_array[ref_id, ref_id:] = sub_id_array[:, 1]
            self._pair_delta_array[ref_id, ref_id:, :] = sub_delta_array
            self._pair_dist_array[ref_id, ref_id:] = sub_dist_array

            # Don't fill the array column wise if we are on the last object
            # to avoid array overrun.
            if ref_id < self._n_reference - 1:
                self._pair_id_array[ref_id + 1:, ref_id] = sub_id_array[:, 1]
                self._pair_delta_array[ref_id + 1:, ref_id, :] = \
                    sub_delta_array
                self._pair_dist_array[ref_id + 1:, ref_id] = sub_dist_array

            # Sort each row on distance for fast look up of pairs given
            # the id of one of the objects in the pair.
            sorted_pair_dist_args = self._pair_dist_array[ref_id, :].argsort()
            self._pair_dist_array[ref_id, :] = self._pair_dist_array[
                ref_id, sorted_pair_dist_args]
            self._pair_id_array[ref_id, :] = self._pair_id_array[
                ref_id, sorted_pair_dist_args]
            self._pair_delta_array[ref_id, :, :] = self._pair_delta_array[
                ref_id, sorted_pair_dist_args, :]

        # Concantent our arrays together.
        unsorted_id_array = np.concatenate(sub_id_array_list)
        unsorted_delta_array = np.concatenate(sub_delta_array_list)
        unsorted_dist_array = np.concatenate(sub_dist_array_list)

        # Sort each array on the pair distances for the inital
        # optimistic pattern matcher lookup.
        sorted_dist_args = unsorted_dist_array.argsort()
        self._dist_array = unsorted_dist_array[sorted_dist_args]
        self._id_array = unsorted_id_array[sorted_dist_args]
        self._delta_array = unsorted_delta_array[sorted_dist_args]

        return None

    def match(self, source_array, n_check, n_match, n_agree,
              max_n_patterns, max_shift, max_rotation, max_dist,
              min_matches, pattern_skip_array=None):
        """Match a given source catalog into the loaded reference catalog.

        Given array of points on the unit sphere and tolerances, we
        attempt to match a pinwheel like pattern between these input sources
        and the reference objects this class was created with. This pattern
        informs of the shift and rotation needed to align the input source
        objects into the frame of the refernces.

        Parameters
        ----------
        source_array: float array
            An array of spherical x,y,z coordinates and a magnitude in units of
            objects having a lower value for sorting. The array should be of
            shape (N, 4).
        n_check : int value
            Number of sources to create a pattern from. Not all objects may be
            checked if n_match criteria is before looping through all n_check
            objects.
        n_match : int value
            Number of objects to use in constructing a pattern to match.
        n_agree: int value
            Number of found patterns that must agree on their shift and
            rotation before exiting. Set this value to 1 to recover the
            expected behavior of Optimistic Pattern Matcher B.
        max_n_patters : int value
            Number of patterns to create from the input source objects to
            attempt to match into the reference objects.
        max_shift: float value
            Maximum allowed shift to match patterns in arcseconds.
        max_rotation: float value
            Maximum allowed rotation between patterns in degrees.
        max_dist: float value
            Maximum distance in arcseconds allowed between candidate spokes in
            the source and reference objects. Also sets that maximum distance
            in the intermediate verify, pattern shift/rotation agreement, and
            final verify steps.
        pattern_skip_array: int array
            Patterns we would like to skip. This could be due to the pattern
            being matched on a pervious iteration that we now consider invalid.
            This assumes the ordering of the source objects is the same between
            different runs of the matcher which, assuming no object has been
            inserted or the magnitudes have changed, it should be.

        Returns
        -------
        output_struct : pipe.base.struct
            A pipebase struct containing the following outputs.

                matches : int array
                    (N, 2) array of matched ids for pairs. Empty list if no
                    match found.
                distances : float array
                    Radian distances between the matched objects. Empty list if
                    no match found.
                pattern_idx : int
                    Index of matched pattern. None if no match found.
                shift : float
                    Magnitude for the shift between the source and
                    reference objects in arcseconds. None if no match found.
        """

        # Given our input source_array we sort on magnitude.
        sorted_source_array = source_array[source_array[:, -1].argsort(), :3]
        n_source = len(sorted_source_array)

        # Initialize output struct.
        output_match_struct = pipeBase.Struct(
            matches=[],
            distances=[],
            pattern_idx=None,
            shift=None,)

        if n_source <= 0:
            self.log.warn("Source object array is empty. Unable to match. "
                          "Exiting matcher.")
            return None

        # To test if the shifts and rotations we find agree with each other,
        # we first create two test points situated at the top and bottom of
        # where the z axis on the sphere bisects the source catalog.
        test_vect_list = self._compute_test_vectors(source_array[:, :3])

        # We now create an empty list of our resultant rotated vectors to
        # compare the different rotations we find.
        rot_vect_list = []

        # Convert the tolerances to values we will use in the code.
        max_cos_shift = np.cos(np.radians(max_shift / 3600.))
        max_cos_rot_sq = np.cos(np.radians(max_rotation)) ** 2
        max_dist_rad = np.radians(max_dist / 3600.)

        # Loop through the sources from brightest to faintest, grabbing a
        # chunk of n_check each time.
        for pattern_idx in xrange(np.min((max_n_patterns,
                                          n_source - n_match))):

            # If this pattern is one that we matched on the past but we
            # now want to skip, we do so here.
            if pattern_skip_array is not None and \
               np.any(pattern_skip_array == pattern_idx):
                self.log.debug("Skipping previously matched bad pattern %i..." %
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
               construct_return_struct.shift_rot_matrix is None or\
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
            for test_vect in test_vect_list:
                tmp_rot_vect_list.append(np.dot(shift_rot_matrix, test_vect))
            tmp_rot_vect_list.append(pattern_idx)
            rot_vect_list.append(tmp_rot_vect_list)

            # Test if we have enough rotations, which agree, or if we
            # are in optimistic mode.
            if self._test_rotation_agreement(rot_vect_list, max_dist_rad) < \
               n_agree - 1:
                continue

            # Perform final verify.
            match_sources_struct = self._match_sources(
                    source_array[:, :3], shift_rot_matrix, max_dist_rad)

            # Check that we have enough matches.
            if len(match_sources_struct.matches) >= min_matches:
                # Convert the observed shift to arcseconds
                shift = np.degrees(np.arccos(cos_shift)) * 3600.
                # Print information to the logger.
                self.log.debug("Succeeded after %i patterns." % pattern_idx)
                self.log.debug("\tShift %.4f arcsec" % shift)
                self.log.debug("\tRotation: %.4f deg" %
                               np.degrees(np.arcsin(sin_rot)))
                # Return.
                output_match_struct.matches = match_sources_struct.matches
                output_match_struct.distances = match_sources_struct.distances
                output_match_struct.pattern_idx = pattern_idx
                output_match_struct.shift = shift

                return output_match_struct

        self.log.warn("Failed after %i patterns." % (pattern_idx + 1))
        return output_match_struct

    def _compute_test_vectors(self, source_array):
        """Compute spherical 3 vectors at the edges of the x, y, z extent
        of the input source catalog.

        Parameters
        ----------
        source_array : float array (N, 3)
            array of 3 vectors representing possitions on the unit
            sphere.

        Returns
        -------
        list of 3 vectors
            list of vectors representing the maxinum extents in x, y, z
            of the input source array. These are used with the rotations
            the code finds to test for agreement from different patterns
            when the code is running in pessimistic mode.
        """

        # Get the center of source_array.
        if np.any(np.logical_not(np.isfinite(source_array))):
            self.log.warn("Input source objects contain non-finite values. "
                          "This could end badly.")
        center_vect = np.nanmean(source_array, axis=0)

        # So that our rotation test works over the full sky we compute
        # the max extent in each cartesian direction x,y,z.
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
        return [xbtm_vect, xtop_vect, ybtm_vect, ytop_vect,
                zbtm_vect, ztop_vect]

    def _construct_pattern_and_shift_rot_matrix(self, src_pattern_array,
                                                n_match, max_cos_theta_shift,
                                                max_cos_rot_sq, max_dist_rad):
        """Test an input source pattern against the reference catalog.

        Returns the candidate matched patterns and their
        implied rotation matrices or None.

        Parameters
        ----------
        src_pattern_array : float array
            Sub selection of source 3 vectors to create a pattern from
        n_match : int
            Number of points to attempt to create a pattern from. Must be
            >= len(src_pattern_array)
        max_cos_theta_shift : float
            Maximum shift allowed between two patterns' centers.
        max_cos_rot_sq : float
            Maximum rotation beteween two patterns that have been shifted
            to have their centers on top of each other.
        max_dist_rad : float
            Maximum delta distance allowed between the source and reference
            pair distances to consider the reference pair a candidate for
            the source pair. Also sets the tolerance between the opening
            angles of the spokes when compared to the reference.

        Return
        -------
        lsst.pipe.base.Struct
            Return a Struct containing the following data:

            ref_candidates : list of ints
                ids of the matched pattern in the internal reference_array
                object.
            src_candidates : list of ints
                Pattern ids of the sources matched.
            shift_rot_matrix : float array
                3x3 matrix specifying the full shift
                and rotation between the reference and source objects. Rotates
                source into reference frame. None if match is not found.
            cos_shift : float
                Magnitude of the shift found between the two patten
                centers. None if match is not found.
            sin_rot : float value of the rotation to align the already shifted
               source pattern to the reference pattern. None if no match found.
        """

        # Create our place holder variables for the matched sources and
        # referces. The source list starts with the 0th and first indexed
        # objects as we are garunteed to use those and these define both
        # the shift and rotation of the final pattern.
        output_matched_pattern = pipeBase.Struct(
            ref_candidates=[],
            src_candidates=[0, 1],
            shift_rot_matrix=None,
            cos_shift=None,
            sin_rot=None)

        # Create the delta vectors and distances we will need to assemble the
        # spokes of the pattern.
        src_delta_array = np.empty((len(src_pattern_array) - 1, 3))
        src_delta_array[:, 0] = (src_pattern_array[1:, 0] -
                                 src_pattern_array[0, 0])
        src_delta_array[:, 1] = (src_pattern_array[1:, 1] -
                                 src_pattern_array[0, 1])
        src_delta_array[:, 2] = (src_pattern_array[1:, 2] -
                                 src_pattern_array[0, 2])
        src_dist_array = np.sqrt(src_delta_array[:, 0]**2 +
                                 src_delta_array[:, 1]**2 +
                                 src_delta_array[:, 2]**2)

        # Our first test. We search the reference dataset for pairs
        # that have the same length as our first source pairs to with
        # plus/minus the max_dist tolerance.
        ref_dist_index_array = self._find_candidate_reference_pairs(
            src_dist_array[0], self._dist_array, max_dist_rad)

        # Start our loop over the candidate reference objects.
        for ref_dist_idx in ref_dist_index_array:
            # We have two candidates for which reference object corresponds
            # with the source at the center of our pattern. As such we loop
            # over and test both possiblities.
            tmp_ref_pair_list = self._id_array[ref_dist_idx]
            for pair_idx, ref_id in enumerate(tmp_ref_pair_list):
                output_matched_pattern.src_candidates = [0, 1]
                output_matched_pattern.ref_candidates = []
                output_matched_pattern.shift_rot_matrix = None
                output_matched_pattern.cos_shift = None
                output_matched_pattern.sin_rot = None
                # Test the angle between our candidate ref center and the
                # source center of our pattern. This angular distance also
                # defines the shift we will later use.
                ref_center = self._reference_array[ref_id]
                cos_shift = np.dot(src_pattern_array[0], ref_center)
                if cos_shift < max_cos_theta_shift:
                    continue

                # Store the current cos shift.
                output_matched_pattern.cos_shift = cos_shift

                # We can now append this one as a candaite.
                output_matched_pattern.ref_candidates.append(ref_id)
                ref_delta = self._delta_array[ref_dist_idx]
                # If the candidate reference center we found is second in
                # this pair we need to reverse the direction of the
                # corresponding pair's delta vector.
                if pair_idx == 0:
                    output_matched_pattern.ref_candidates.append(
                        tmp_ref_pair_list[1])
                else:
                    output_matched_pattern.ref_candidates.append(
                        tmp_ref_pair_list[0])
                    ref_delta *= -1

                # For dense fields it will be faster to compute the absolute
                # rotation this pair suggests first rather than saving it
                # after all the spokes are found. We then compute the cos^2
                # of the rotation and first part of the rotation matrix from
                # source to reference frame.
                test_rot_struct = self._test_rotation(
                    src_pattern_array[0], ref_center, src_delta_array[0],
                    ref_delta, cos_shift, max_cos_rot_sq)
                if test_rot_struct.cos_rot_sq is None or \
                   test_rot_struct.shift_matrix is None:
                    continue

                # Get the data from the return struct.
                cos_rot_sq = test_rot_struct.cos_rot_sq
                shift_matrix = test_rot_struct.shift_matrix

                # Now that we have a candidate first spoke and reference
                # pattern center, we mask our future search to only those
                # pairs that contain our candidate reference center.
                ref_dist = self._dist_array[ref_dist_idx]
                tmp_ref_delta_array = self._pair_delta_array[ref_id]
                tmp_ref_dist_arary = self._pair_dist_array[ref_id]
                tmp_ref_id_array = self._pair_id_array[ref_id]

                # Now we feed this sub data to match the spokes of
                # our pattern.
                pattern_spoke_struct = self._create_pattern_spokes(
                    src_pattern_array[0], src_delta_array, src_dist_array,
                    self._reference_array[ref_id], ref_id, ref_delta,
                    ref_dist, tmp_ref_delta_array, tmp_ref_dist_arary,
                    tmp_ref_id_array, max_dist_rad,
                    n_match)

                # If we don't find enough candidates we can continue to the
                # next reference center pair.
                if len(pattern_spoke_struct.ref_spoke_list) < n_match - 2 or \
                   len(pattern_spoke_struct.src_spoke_list) < n_match - 2:
                    continue

                # If we have the right number of matched ids we store these.
                output_matched_pattern.ref_candidates.extend(
                    pattern_spoke_struct.ref_spoke_list)
                output_matched_pattern.src_candidates.extend(
                    pattern_spoke_struct.src_spoke_list)

                # We can now create our full rotation matrix for both the
                # shift and rotation. Reminder shift, aligns the pattern
                # centers, rotation rotates the spokes on top of each other.
                shift_rot_struct = self._create_shift_rot_matrix(
                    cos_rot_sq, shift_matrix, src_delta_array[0],
                    self._reference_array[ref_id], ref_delta)
                # If we fail to create the rotation matrix, continue to the
                # next objects.
                if shift_rot_struct.sin_rot is None or \
                   shift_rot_struct.shift_rot_matrix is None:
                    continue

                # Get the data from the return struct.
                output_matched_pattern.sin_rot = shift_rot_struct.sin_rot
                output_matched_pattern.shift_rot_matrix = \
                    shift_rot_struct.shift_rot_matrix

                # Now that we have enough candidates we test to see if it
                # passes intermediate verify. This shifts and rotates the
                # source pattern into the reference frame and tests that each
                # source/reference object pair is within max_dist.
                if self._intermediate_verify(
                        src_pattern_array[
                            output_matched_pattern.src_candidates],
                        self._reference_array[
                            output_matched_pattern.ref_candidates],
                        output_matched_pattern.shift_rot_matrix, max_dist_rad):
                    return output_matched_pattern

        return output_matched_pattern

    def _find_candidate_reference_pairs(self, src_dist, ref_dist_array,
                                        max_dist_rad):
        """Wrap numpy.searchsorted to find the range of reference spokes within
        a spoke distance tolerance of our source spoke.

        Returns an array sorted from the smallest absolute delta distance
        between source and reference spoke length. This sorting increases the
        speed for the pattern search greatly.

        Parameters
        ----------
        src_dist : float radians
            float value of the distance we would like to search for in
            the reference array in radians.
        ref_dist_array : float array
            sorted array of distances in radians.
        max_dist_rad : float
            maximum plus/minus search to find in the reference array in
            radians.

        Return
        ------
        int array
            indices lookup into the input ref_dist_array sorted by the
            difference in value to the src_dist from absolute value
            smallest to largest.
        """
        # Find the index of the minimum and maximum values that satisfy
        # the tolerance.
        start_idx = np.searchsorted(ref_dist_array, src_dist - max_dist_rad)
        end_idx = np.searchsorted(ref_dist_array, src_dist + max_dist_rad,
                                  side='right')

        # If these are equal there are no candidates and we exit.
        if start_idx == end_idx:
            return []

        # Make sure the endpoints of the input array are respected.
        if start_idx < 0:
            start_idx = 0
        if end_idx > ref_dist_array.shape[0]:
            end_idx = ref_dist_array.shape[0]

        # Now we sort the indices from smallest absolute delta dist difference
        # to the largest and return the vector. This step greatly increases the
        # speed of the algorithm.
        tmp_diff_array = np.fabs(ref_dist_array[start_idx:end_idx] - src_dist)
        return tmp_diff_array.argsort() + start_idx

    def _test_rotation(self, src_center, ref_center, src_delta, ref_delta,
                       cos_shift, max_cos_rot_sq):
        """ Test if the rotation implied between the source
        pattern and reference pattern is within tolerance. To test this
        we need to create the first part of our spherical rotation matrix
        which we also return for use later.

        Parameters
        ----------
        src_center : float array3
            pattern.
        ref_center : float array
            3 vector defining the center of the candidate refence pinwheel
            pattern.
        src_delta : float array
            3 vector delta between the source pattern center and the end of
            the pinwheel spoke.
        ref_delta : float array
            3 vector delta of the candidate matched reference pair
        cos_shift : float
            Cosine of the angle between the source and reference candidate
            centers.
        max_cos_rot_sq : float
            candidate reference pair after shifting the centers on top of each
            other. The function will return None if the rotation implied is
            greater than max_cos_rot_sq.

        Returns
        -------
        lsst.pipe.base.Struct
            Return a pipe.base.Struct containing the following data.

            cos_rot_sq : float
                magnitude of the rotation needed to align the two patterns
                after their centers are shifted on top of each other.
                None if rotation test fails.
            shift_matrix : float array
                3x3 rotation matrix describing the shift needed to align
                the source and candidate reference center.
                None if rotation test fails.
        """

        # Make sure the sine is a real number.
        if cos_shift > 1.0:
            cos_shift = 1.
        elif cos_shift < -1.0:
            cos_shift = -1.
        sin_shift = np.sqrt(1 - cos_shift ** 2)

        # If the sine of our shift is zero we only need to use the idenity
        # matrix for the shift. Else we construct the rotation matrix for
        # shift.
        if sin_shift > 0:
            rot_axis = np.cross(src_center, ref_center)
            rot_axis /= sin_shift
            shift_matrix = self._create_spherical_rotation_matrix(
                rot_axis, cos_shift, sin_shift)
        else:
            shift_matrix = np.identity(3)

        # Now that we have our shift we apply it to the src delta vector
        # and check the rotation.
        rot_src_delta = np.dot(shift_matrix, src_delta)
        cos_rot_sq = (np.dot(rot_src_delta, ref_delta)**2 /
                      (np.dot(rot_src_delta, rot_src_delta) *
                       np.dot(ref_delta, ref_delta)))
        # If the rotation isn't in tolerance return None.
        if cos_rot_sq < max_cos_rot_sq:
            return pipeBase.Struct(
                cos_rot_sq=None,
                shift_matrix=None,)
        return pipeBase.Struct(
            cos_rot_sq=cos_rot_sq,
            shift_matrix=shift_matrix,)

    def _create_spherical_rotation_matrix(self, rot_axis, cos_rotation,
                                          sin_rotion):
        """Construct a generalized 3D rotation matrix about a given
        axis.

        Parameters
        ----------
        rot_axis : float array
            3 vector defining the axis to rotate about.
        cos_rotation : float
            cosine of the rotation angle.
        sin_rotation : float
            sine of the rotation angle.

        Return
        ------
        float array
            3x3 spherical, rotation matrix.
        """

        rot_cross_matrix = np.array(
            [[0., -rot_axis[2], rot_axis[1]],
             [rot_axis[2], 0., -rot_axis[0]],
             [-rot_axis[1], rot_axis[0], 0.]], dtype=np.float64)
        shift_matrix = (cos_rotation*np.identity(3) +
                        sin_rotion*rot_cross_matrix +
                        (1. - cos_rotation)*np.outer(rot_axis, rot_axis))

        return shift_matrix

    def _create_pattern_spokes(self, src_ctr, src_delta_array, src_dist_array,
                               ref_ctr, ref_ctr_id, ref_delta, ref_dist,
                               ref_delta_array, ref_dist_array,
                               ref_id_array, max_dist_rad, n_match):
        """ Create the individual spokes that make up the pattern now that the
        shift and rotation are within tolerance.

        If we can't create a valid pattern we exit early.

        Parameters
        ----------
        src_ctr : float array
            3 vector of the source pinwheel center
        src_delta_array : float array
            Array of 3 vector deltas between the source center and the pairs
            that make up the remaining spokes of the pinwheel
        src_dist_array : float array
            Array of the distances of each src_delta in the pinwheel
        ref_ctr : float array
            3 vector of the candidate refenerce center
        ref_ctr_id : int
            id of the ref_ctr in the master reference array
        ref_delta : float array
            3 vector of the first candidate pair of the pinwheel. This is
            the candidate pair that was matched in the
            main _construct_pattern_and_shift_rot_matrix loop
        ref_dist : float
            Radian distance of the first candidate reference pair
        ref_delta_array : float array
            Array of 3 vector deltas that are have the current candidate
            reference center as part of the pair
        ref_dist_array : float array
            Array of vector distances for each of the reference pairs
        ref_id_array : int array
            Array of id lookups into the master reference array that our
            center id object is paired with.
        max_dist_rad : float
            Maximum search distance
        n_match : int
            Number of source deltas that must be matched into the reference
            deltas in order to consider this a successful pattern match.

        Returns
        -------
        lsst.pipe.base.Struct
            The Struct contains the following data:

            ref_spoke_list : list of ints specifying ids into the master
                reference array
            src_spoke_list : list of ints specifying indices into the current
                source pattern that is being tested.
        """
        # Struct where we will be putting our results.
        output_spokes = pipeBase.Struct(
            ref_spoke_list=[],
            src_spoke_list=[],)

        # Counter for number of spokes we failed to find a reference
        # candidate for. We break the loop if we haven't found enough.
        n_fail = 0
        # Loop over the source pairs.
        for src_idx in xrange(1, len(src_dist_array)):
            if n_fail > len(src_dist_array) - (n_match - 1):
                break

            # Given our length tolerance we can use it to compute a tolernace
            # on the angle between our spoke.
            src_sin_tol = (max_dist_rad /
                           (src_dist_array[src_idx] + max_dist_rad))
            # Test if the small angle approximation will still hold. This is
            # defined as when sin(theta) ~= theta to within 0.1% of each other.
            # This also implicitly sets a minimum spoke length that we can use.
            if src_sin_tol > 0.0447:
                n_fail += 1
                continue

            # Find the reference pairs that include our candidate pattern
            # center and sort them in increasing delta
            ref_dist_idx_array = self._find_candidate_reference_pairs(
                src_dist_array[src_idx], ref_dist_array, max_dist_rad)

            # Test the spokes and return the id of the reference object.
            # Return None if no match is found.
            ref_id = self._test_spoke(
                src_ctr, src_delta_array[src_idx], src_dist_array[src_idx],
                src_delta_array[0], src_dist_array[0], ref_ctr, ref_ctr_id,
                ref_delta, ref_dist, ref_dist_idx_array, ref_delta_array,
                ref_dist_array,
                ref_id_array, src_sin_tol)
            if ref_id is None:
                n_fail += 1
                continue

            # Append the successful indices to our list. The src_idx needs
            # an extra iteration to skip the first and second source objects.
            output_spokes.ref_spoke_list.append(ref_id)
            output_spokes.src_spoke_list.append(src_idx + 1)
            # If we found enough reference objects we can return early. This is
            # n_match - 2 as we already have 2 source objects matched into the
            # reference data.
            if len(output_spokes.ref_spoke_list) >= n_match - 2:
                break
        return output_spokes

    # TODO: Possibly clean up the arguments here by pre computing
    # the source dot and cross products of the opening angle in the
    # previous method.
    def _test_spoke(self, src_ctr, src_delta, src_dist, src_ctr_delta,
                    src_ctr_dist, ref_ctr, ref_ctr_id, ref_delta, ref_dist,
                    ref_dist_idx_array, ref_delta_array, ref_dist_array,
                    ref_id_array, src_sin_tol):
        """Test the opening angle between the first spoke of our pattern
        for the soruce object against the reference object.

        This method makes heavy use of the small angle approximation to perform
        the comparison.

        Parameters
        ----------
        src_ctr : float array
            3 vector of the source pinwheel center
        src_delta : float array
            3 vector delta from the source center and the source object that
            makes up the current spoke of the pinwheel we are testing.
        src_dist : float array
            Distance of the current spoke we are testing
        src_ctr_delta : float array
            3 vector delta between the center of the pattern and the first
            spoke of the pattern. Used to test compute the opening angle
            between the current spoke and the first spoke.
        src_ctr_dist : float
            Distance between the pairs that make up src_ctr_delta
        ref_ctr : float array
            3 vector of the candidate reference center
        ref_ctr_id : int
            id lookup of the ref_ctr into the master reference array
        ref_delta : float array
            3 vector of the first candidate pair of the pinwheel. That is
            the candidate pair that was matched in the
            main _construct_pattern_and_shift_rot_matrix loop
        ref_dist : float
            Radian distance of the first candidate reference pair
        ref_dist_idx_array : int array
            Indices sorted by the delta distance between the source
            spoke we are trying to test and the candidate reference
            spokes.
        ref_delta_array : float array
            Array of 3 vector deltas that are have the current candidate
            reference center as part of the pair
        ref_dist_array : float array
            Array of vector distances for each of the reference pairs
        ref_id_array : int array
            Array of id lookups into the master reference array that our
            center id object is paired with.
        src_sin_tol : float
            Sine of tolerance allowed between source and reference spoke
            opening angles.

        Returns
        -------
        None or int
            If we can not find a candidate spoke we return None else we
            return an int id into the master reference array.
        """

        # Precompute all of the source only cross and dot products so we don't
        # have to do it for each interation in the reference loop.
        cos_theta_src = (np.dot(src_delta, src_ctr_delta) /
                         (src_dist * src_ctr_dist))
        cross_src = (np.cross(src_delta, src_ctr_delta) /
                     (src_dist * src_ctr_dist))
        dot_cross_src = np.dot(cross_src, src_ctr)

        # Loop over our candidate reference objects.
        for ref_dist_idx in ref_dist_idx_array:
            # Check the direction of the delta vector.
            ref_sign = 1
            if ref_id_array[ref_dist_idx] < ref_ctr_id:
                ref_sign = -1

            # Compute the cos between our "center" reference vector and the
            # current reference candidate.
            cos_theta_ref = ref_sign * (
                np.dot(ref_delta_array[ref_dist_idx], ref_delta) /
                (ref_dist_array[ref_dist_idx] * ref_dist))

            # Make sure we can safely make the comparison in case
            # our "center" and candidate vectors are mostly aligned.
            if cos_theta_ref ** 2 < (1 - src_sin_tol ** 2):
                cos_sq_comparison = ((cos_theta_src - cos_theta_ref) ** 2 /
                                     (1 - cos_theta_ref ** 2))
            else:
                cos_sq_comparison = ((cos_theta_src - cos_theta_ref) ** 2 /
                                      src_sin_tol ** 2)
            # Test the difference of the cosine of the reference angle against
            # the source angle. Assumes that the delta between the two is
            # small.
            if cos_sq_comparison > src_sin_tol ** 2:
                continue

            # The cosine tests the manitude of the angle but not
            # its direction. To do that we need to know the sine as well.
            # This cross product calculation does that.
            cross_ref = ref_sign * (
                np.cross(ref_delta_array[ref_dist_idx], ref_delta) /
                (ref_dist_array[ref_dist_idx] * ref_dist))
            dot_cross_ref = np.dot(cross_ref, ref_ctr)

            # Check the value of the cos again to make sure that it is not
            # near zero.
            if abs(cos_theta_src) < src_sin_tol:
                sin_comparison = (dot_cross_src - dot_cross_ref) / src_sin_tol
            else:
                sin_comparison = \
                    (dot_cross_src - dot_cross_ref) / cos_theta_ref

            if abs(sin_comparison) > src_sin_tol:
                continue

            # Return the correct id of the candidate we found.
            return ref_id_array[ref_dist_idx]

        return None

    def _create_shift_rot_matrix(self, cos_rot_sq, shift_matrix, src_delta,
                                 ref_ctr, ref_delta):
        """ Create the final part of our spherical rotation matrix.

        Parameters
        ----------
        cos_rot_sq : float
            cosine of the rotation needed to align our source and reference
            candidate patterns.
        shift_matrix : float array
            3x3 rotation matrix for shifting the source pattern center on top
            of the candidate reference pattern center.
        src_delta : float array
            3 vector delta of representing the first spoke of the source
            pattern
        ref_ctr : float array
            3 vector on the unitsphere representing the center of our
            reference pattern.
        ref_delta : float array
            3 vector delta made by the first pair of the reference pattern.

        Returns
        -------
        lsst.pipe.base.Struct
            Struct object containing the following data:

            sin_rot : float sine of the amount of rotation between the
                source and reference pattern. We use sine here as it is
                signed and tells us the chirality of the rotation.
            shift_rot_matrix : float array representing the 3x3 rotation
                matrix that takes the source patern and shifts and rotates
                it to align with the reference pattern.
        """
        cos_rot = np.sqrt(cos_rot_sq)
        rot_src_delta = np.dot(shift_matrix, src_delta)
        delta_dot_cross = np.dot(np.cross(rot_src_delta, ref_delta), ref_ctr)

        sin_rot = np.sign(delta_dot_cross) * np.sqrt(1 - cos_rot_sq)
        rot_matrix = self._create_spherical_rotation_matrix(
            ref_ctr, cos_rot, sin_rot)

        shift_rot_matrix = np.dot(rot_matrix, shift_matrix)

        return pipeBase.Struct(
            sin_rot=sin_rot,
            shift_rot_matrix=shift_rot_matrix,)

    def _intermediate_verify(self, src_pattern, ref_pattern, shift_rot_matrix,
                             max_dist_rad):
        """ Perform an intermediate verify step. Rotate the matches references
        into the source frame and test their distances against tolerance. Only
        return true if all points are within tolerance.

        Parameters
        ----------
        src_pattern : float array
            Array of 3 vectors representing the points that make up our source
            pinwheel pattern.
        ref_pattern : float array
            Array of 3 vectors representing our candidate reference pinwheel
            pattern.
        shift_rot_matrix : float array
            3x3 rotation matrix that takes the source objects and rotates them
            onto the frame of the reference objects
        max_dist_rad : float
            Maximum distance allowed to consider two objects the same.

        Returns
        -------
        bool
           Return true if all of the points in our source pattern are within
           max_dist_rad of their matched reference objects.
        """
        if len(src_pattern) != len(ref_pattern):
            raise ValueError(
                "Source pattern length does not match ref pattern.\n"
                "\t source pattern len=%i, reference pattern len=%i" %
                (len(src_pattern), len(ref_pattern)))
        shifted_ref_pattern = np.dot(shift_rot_matrix.transpose(),
                                     ref_pattern.transpose()).transpose()
        tmp_delta_array = src_pattern - shifted_ref_pattern
        tmp_dist_array = (tmp_delta_array[:, 0] ** 2 +
                          tmp_delta_array[:, 1] ** 2 +
                          tmp_delta_array[:, 2] ** 2)
        is_good_shift_rot = np.all(tmp_dist_array < max_dist_rad ** 2)
        return is_good_shift_rot

    def _test_rotation_agreement(self, rot_vect_list, max_dist_rad):
        """ Test this rotation against the previous N found and return
        the number that a agree within tolerance to where our test
        points are.

        Parameters
        ----------
        rot_vect_list : list of lists of float arrays
            list of lists of rotated 3 vectors representing the maximum x, y, z
            extent on the unit sphere of the input source objects roated by the
            candidate rotations into the reference frame.
        max_dist_rad : float
            maximum distance in radians to consider two points "agreeing" on
            a rotation

        Returns
        -------
        int
            Number of candidate rotations that agree for all of the rotated
            test 3 vectors.
        """

        self.log.debug("Comparing pattern %i to previous %i rotations..." %
                       (rot_vect_list[-1][-1], len(rot_vect_list) - 1))

        tot_consent = 0
        for rot_idx in xrange(max((len(rot_vect_list) - 1), 0)):
            tmp_dist_list = []
            for vect_idx in xrange(len(rot_vect_list[rot_idx]) - 1):
                tmp_delta_vect = (rot_vect_list[rot_idx][vect_idx] -
                                  rot_vect_list[-1][vect_idx])
                tmp_dist_list.append(
                    np.dot(tmp_delta_vect, tmp_delta_vect))
            if np.all(np.array(tmp_dist_list) < max_dist_rad ** 2):
                tot_consent += 1
        return tot_consent

    def _match_sources(self, source_array, shift_rot_matrix, max_dist_rad):
        """ Shift both the reference and source catalog to the the respective
        frames and find their nearest neighbor using a kdTree. Removes all
        matches who do not agree when either the refernce or source catalog is
        rotated and removes all matches greated than the requested distance.

        Parameters
        ----------
        source_array : float array
            array of 3 vectors representing the source objects we are trying to
            match into the source catalog.
        shift_rot_matrix : float array
            3x3 rotation matrix that performs the spherical rotation from the
            source frame into the reference frame.
        max_dist_rad : float
            maximum distance to consider an object a match after rotating the
            the source and references into the other's frame.

        Returns
        -------
        lsst.pipe.base.Struct
            A Struct object containing the following data
            matches : a (N, 2) array of integer ids into the source and
                reference arrays. Matches are only returned for those that
                satisfy the distance and handshake criteria.
            distances : float array of the distance between each match in
                radians after the shift and rotation is applied.
        """
        shifted_references = np.dot(
            shift_rot_matrix.transpose(),
            self._reference_array.transpose()).transpose()
        shifted_sources = np.dot(
            shift_rot_matrix,
            source_array.transpose()).transpose()

        ref_matches = np.empty((len(shifted_references), 2),
                               dtype=np.uint32)
        src_matches = np.empty((len(shifted_sources), 2),
                               dtype=np.uint32)

        ref_matches[:, 1] = np.arange(len(shifted_references),
                                      dtype=np.uint32)
        src_matches[:, 0] = np.arange(len(shifted_sources),
                                      dtype=np.uint32)

        ref_kdtree = cKDTree(self._reference_array)
        src_kdtree = cKDTree(source_array)

        tmp_src_dist, tmp_src_idx = src_kdtree.query(shifted_references)
        tmp_ref_dist, tmp_ref_idx = ref_kdtree.query(shifted_sources)

        ref_matches[:, 0] = tmp_src_idx
        src_matches[:, 1] = tmp_ref_idx

        handshake_mask = self._handshake_match(ref_matches, src_matches)
        final_mask = np.logical_and(handshake_mask,
                                    tmp_src_dist < max_dist_rad)
        return pipeBase.Struct(
            matches=ref_matches[final_mask],
            distances=tmp_src_dist[final_mask],)

    def _handshake_match(self, matches_ref, matches_src):
        """Return only those matches where both the source
        and reference objects agree they they are each others'
        nearest neighbor.

        Parameters
        ----------
        matches_ref : int array
            (N, 2) int array of nearest neighbor matches between shifted and
            rotated source objects matched into the references.
        matches_src : int array
            (M, 2) int array of nearest neighbor matches between shifted and
            rotated reference objects matched into the sources.

        Return
        ------
        bool array
           Return the array positions where the two match catalogs agree.
        """
        handshake_mask_array = np.zeros(len(matches_ref))

        for ref_match_idx, match in enumerate(matches_ref):
            src_match_idx = np.searchsorted(matches_src[:, 0], match[0])
            if match[1] == matches_src[src_match_idx, 1]:
                handshake_mask_array[ref_match_idx] = True
        return handshake_mask_array
