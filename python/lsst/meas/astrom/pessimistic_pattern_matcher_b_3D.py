from __future__ import division, print_function, absolute_import

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
    input datasets are n-spoked pinwheels. Refer to DMTN #031 for more details.
    http://github.com/lsst-dm/dmtn-031
    ----------------------------------------------------------------------------
    Attributes:
        reference_catalog: Input array of spherical points x, y, z of to use as
        reference objects for pattern matching.
    """

    def __init__(self, reference_catalog, log):

        self._reference_catalog = reference_catalog
        self._n_reference = len(self._reference_catalog)
        self.log = log

        self._build_distances_and_angles()

    def _build_distances_and_angles(self):
        """ Create the data structures we will use to search
        for our match in. Throughout this function and the rest of the class
        we use id to reference the possition in the input reference catalog and
        index to 'index' into the arrays sorted on distance.
        """

        # Initialize the arrays we will need for quick look up of pairs once
        # have a candidate spoke center.
        self._pair_id_array = np.empty(
            (self._n_reference, self._n_reference - 1, 2),
            dtype=np.int)
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
        # vector deltas.
        for ref_id, ref_obj in enumerate(self._reference_catalog):

            # Reserve and fill the ids of each reference object pair.
            sub_id_array = np.zeros((self._n_reference - 1 - ref_id, 2),
                                    dtype=np.int32)
            sub_id_array[:, 0] = ref_id
            sub_id_array[:, 1] = np.arange(ref_id + 1, self._n_reference,
                                           dtype=np.int32)

            # Compute the vector deltas for each pair of reference objects
            # and compute and store the distances.
            sub_delta_array = (self._reference_catalog[ref_id + 1:, :] -
                               ref_obj)
            sub_dist_array = np.sqrt(sub_delta_array[:, 0] ** 2 +
                                     sub_delta_array[:, 1] ** 2 +
                                     sub_delta_array[:, 2] ** 2)

            # Append to our arrays to the output lists for later concatenation.
            sub_id_array_list.append(sub_id_array)
            sub_delta_array_list.append(sub_delta_array)
            sub_dist_array_list.append(sub_dist_array)

            # Fill the pair look up arrays row wise and then column wise.
            self._pair_id_array[ref_id, ref_id:, :] = sub_id_array
            self._pair_delta_array[ref_id, ref_id:, :] = sub_delta_array
            self._pair_dist_array[ref_id, ref_id:] = sub_dist_array

            # Don't fill the array column wise if we are on the last object
            # to avoid overrun.
            if ref_id < self._n_reference - 1:
                self._pair_id_array[ref_id + 1:, ref_id, :] = sub_id_array
                self._pair_delta_array[ref_id + 1:, ref_id, :] = \
                    sub_delta_array
                self._pair_dist_array[ref_id + 1:, ref_id] = sub_dist_array

            # Sort each row on distance for fast look up of pairs given
            # one of the objects present.
            sorted_pair_dist_args = self._pair_dist_array[ref_id, :].argsort()
            self._pair_dist_array[ref_id, :] = self._pair_dist_array[
                ref_id, sorted_pair_dist_args]
            self._pair_id_array[ref_id, :, :] = self._pair_id_array[
                ref_id, sorted_pair_dist_args, :]
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

    def match(self, source_catalog, n_check, n_match, n_agree,
              max_n_patterns, max_shift, max_rotation, max_dist,
              min_matches, pattern_skip_array=None):
        r"""Match a given source catalog into the loaded reference catalog.

        Given array of set of points on the unit sphere and tolerances, we
        attempt to match a pinwheel like pattern between these input sources
        and the reference objects this class was created with.

        Parameters
        ----------
        source_catalog: float64_array
            An array of spherical x,y,z coordinates and a magnitude in units of
            objects having a lower value for sorting.
        n_check : int value
            Number of sources to attempt a match on. Not all may be checked if
            n_match criteria is met before hand.
        n_match : int value
            Number of objects to use in constructing a pattern.
        n_agree: int value
            Number of found patterns that must agree on their shift and
            rotation before exiting. Set this value to 1 to recover the
            expected behavior of Optimistic Pattern Matcher B.
        max_n_patters : int value
            Number of patterns to craete on from the input source objects to
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
            A pipebase struct containing the following outputs
                maches : a (N, 2) int array of matched pairs
                distances : a float array of pair distances
                pattern_idx : int index of matched pattern
                shift : a float value for the shift between the source and
                    reference objects.
        Returns None type if no match is found.
        """

        # Given our input source_catalog we sort on magnitude.
        sorted_catalog = source_catalog[source_catalog[:, -1].argsort(), :3]
        n_source = len(sorted_catalog)
        if n_source <= 0:
            self.log.warn("Source object array is empty. Unable to match. "
                          "Exiting matcher.")
            return None

        # To test if the shifts and rotations we find agree we first create
        # two test points situtated at the top and bottom of where the z axis
        # on the sphere bicects the source catalog.
        test_vect_list = self._compute_test_vectors(source_catalog[:, :3])

        # We now create an empty list of our resultant rotator vectors to
        # compare the different roations we find.
        rot_vect_list = []

        # Convert the tolerances to values we will use in the code.
        max_cos_shift = np.cos(np.radians(max_shift/3600.))
        max_cos_rot_sq = np.cos(np.radians(max_rotation)) ** 2
        max_dist_rad = np.radians(max_dist/3600.)

        # Loop through the sources from brightest to faintest grabbing a chucnk
        # of n_check each time.
        for pattern_idx in xrange(np.min((max_n_patterns,
                                          n_source - n_match))):

            # If this pattern is one that we matched on the past but we
            # now want to skip it, we can do so here.
            if pattern_skip_array is not None and \
               np.any(pattern_skip_array == pattern_idx):
                print("Skipping previously matched bad pattern %i..." %
                      pattern_idx)
                continue
            # Grab the sources to attempt to create this pattern.
            pattern = sorted_catalog[
                pattern_idx: np.min((pattern_idx + n_check, n_source)), :3]

            # Construct a pattern given the number of points we are using to
            # create it. We also convert our tolerances to radians at this
            # time.
            construct_return_struct = \
                self._construct_pattern_and_shift_rot_matrix(
                     pattern, n_match, max_cos_shift, max_cos_rot_sq,
                     max_dist_rad)
            if construct_return_struct is None:
                continue

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

            # Test for if we have enough rotations, which agree, or if we
            # are in optimistic mode.
            if self._test_rotation_agreement(rot_vect_list, max_dist_rad) < \
               n_agree - 1:
                continue

            # Perform final verify and check that we have enough matches.
            match_sources_struct = self._match_sources(
                    source_catalog[:, :3], shift_rot_matrix, max_dist_rad)
            shift = np.degrees(np.arccos(cos_shift)) * 3600.
            if len(match_sources_struct.matches) >= min_matches:
                self.log.debug("Succeeded after %i patterns." % pattern_idx)
                self.log.debug("\tShift %.4f arcsec" % shift)
                self.log.debug("\tRotation: %.4f deg" %
                               np.degrees(np.arcsin(sin_rot)))
                return pipeBase.Struct(
                    matches=match_sources_struct.matches,
                    distances=match_sources_struct.distances,
                    pattern_idx=pattern_idx,
                    shift=shift,)

        self.log.warn("Failed after %i patterns." % (pattern_idx + 1))
        return None

    def _compute_test_vectors(self, source_array):
        """Create a list of vectors we can use to test the
        different rotations against.

        Return a list of 3 vector arrays for testing rotation agreement.
        ------------------------------------------------------
        """

        # Get the center of source_array.
        if np.any(np.logical_not(np.isfinite(source_array))):
            self.log.warn("Input source objects contain non-finite values. "
                          "This could end badly.")
        center_vect = np.nanmean(source_array, axis=0)

        # So that our rotation test works overl the full sky we compute
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
        ------------------------------------------------------
        """

        # Create our place hold variables for the matched sources and
        # referces. The source list starts with the 0th and first indexed
        # objects as we are garunteed to use those.
        matched_src_list = [0, 1]
        matched_ref_list = []
        shift_rot_matrix = None

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

        # Our first test.
        ref_dist_index_array = self._find_candidate_reference_pairs(
            src_dist_array[0], self._dist_array, max_dist_rad)

        # Start our loop over the candidate reference objects.
        for ref_dist_idx in ref_dist_index_array:
            # We hav two candidates for which reference object correspondes
            # with the source at the center of our pattern. As such we loop
            # over and test both possiblities.
            tmp_ref_pair_list = self._id_array[ref_dist_idx]
            for pair_idx, ref_id in enumerate(tmp_ref_pair_list):
                matched_src_list = [0, 1]
                matched_ref_list = []
                # Test the angle between our candidate ref center and the
                # source center of our pattern.
                ref_center = self._reference_catalog[ref_id]
                cos_shift = np.dot(src_pattern_array[0], ref_center)
                if cos_shift < max_cos_theta_shift:
                    continue

                # We can now append this one as a candaite.
                matched_ref_list.append(ref_id)
                ref_delta = self._delta_array[ref_dist_idx]
                if pair_idx == 0:
                    matched_ref_list.append(tmp_ref_pair_list[1])
                else:
                    matched_ref_list.append(tmp_ref_pair_list[0])
                    ref_delta *= -1

                # For dense fields it will be faster to compute the absolute
                # rotation this pair suggests first rather than saving it
                # after all the spokes are found. We then compute the cos^2
                # of the rotation and first part of the rotation matrix from
                # source to reference frame.
                test_rot_struct = self._test_rotation(
                    src_pattern_array[0], ref_center, src_delta_array[0],
                    ref_delta, cos_shift, max_cos_rot_sq)
                if test_rot_struct is None:
                    continue

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
                    self._reference_catalog[ref_id], ref_id, ref_delta,
                    ref_dist, tmp_ref_delta_array, tmp_ref_dist_arary,
                    tmp_ref_id_array, max_dist_rad,
                    n_match)

                ref_spoke_list = pattern_spoke_struct.ref_spoke_list
                src_spoke_list = pattern_spoke_struct.src_spoke_list
                # If we don't find enough we can continue to the next reference
                # center pair.
                if len(ref_spoke_list) < n_match - 2 or \
                   len(src_spoke_list) < n_match - 2:
                    continue

                # If we have the right number of matched ids we store these.
                matched_ref_list.extend(ref_spoke_list)
                matched_src_list.extend(src_spoke_list)

                # We can now create our full rotation matrix for both the
                # shift and coordinate system rotation.
                shift_rot_struct = self._create_shift_rot_matrix(
                    cos_rot_sq, shift_matrix, src_delta_array[0],
                    self._reference_catalog[ref_id], ref_delta)
                # If we succeeded, return the pair_ids, the shift/rot matrix
                # and values of the found shift and rotation.
                if shift_rot_struct is None:
                    continue

                sin_rot = shift_rot_struct.sin_rot
                shift_rot_matrix = shift_rot_struct.shift_rot_matrix

                # Now that we have enough canidates we test to see if it passes
                # intermediate verify
                if self._intermediate_verify(src_pattern_array[
                                                 matched_src_list],
                                             self._reference_catalog[
                                                 matched_ref_list],
                                             shift_rot_matrix, max_dist_rad):
                    return pipeBase.Struct(
                        ref_candidates=matched_ref_list,
                        src_candidates=matched_src_list,
                        shift_rot_matrix=shift_rot_matrix,
                        cos_shift=cos_shift,
                        sin_rot=sin_rot)
        return None

    def _find_candidate_reference_pairs(self, src_dist, ref_dist_array,
                                        max_dist_rad):
        """Wrap numpy.searchsorted to find the range of reference spokes within
        a spoke length tolernace of our source spoke. Returns an array sorted
        from the smallest delta distance between source and reference spoke
        length. This sorting increases the speed for the pattern creation
        greatly.
        ------------------------------------------------------
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

        # Now we sort the indicies from smallest delta dist difference to
        # the largest and return the vector. This step greatly increases the
        # speed of the algorithm.
        tmp_diff_array = np.fabs(ref_dist_array[start_idx:end_idx] - src_dist)
        return tmp_diff_array.argsort() + start_idx

    def _test_rotation(self, src_center, ref_center, src_delta, ref_delta,
                       cos_shift, max_cos_rot_sq):
        """ Test if the rotation implied between the source
        pattern and reference pattern is within tolerance. To test this
        we need to create the first part of our spherical rotation matrix
        which we also return for use later.
        ------------------------------------------------------
        """

        # Make sure the sin is a real number.
        if cos_shift > 1.0:
            cos_shift = 1.
        elif cos_shift < -1.0:
            cos_shift = -1.
        sin_shift = np.sqrt(1 - cos_shift ** 2)

        # If our sign shift is zero we only need to use the idenity matrix
        # for the shift.
        if sin_shift > 0:
            rot_axis = np.cross(src_center, ref_center)
            rot_axis /= sin_shift

            rot_cross_matrix = np.array(
                [[0., -rot_axis[2], rot_axis[1]],
                 [rot_axis[2], 0., -rot_axis[0]],
                 [-rot_axis[1], rot_axis[0], 0.]], dtype=np.float64)
            shift_matrix = (cos_shift*np.identity(3) +
                            sin_shift*rot_cross_matrix +
                            (1. - cos_shift)*np.outer(rot_axis, rot_axis))
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
            return None
        return pipeBase.Struct(
            cos_rot_sq=cos_rot_sq,
            shift_matrix=shift_matrix,)

    def _create_pattern_spokes(self, src_ctr, src_delta_array, src_dist_array,
                               ref_ctr, ref_ctr_id, ref_delta, ref_dist,
                               ref_delta_array, ref_dist_array,
                               ref_id_array, max_dist_rad, n_match):
        """ Create the indiviual spokes that make up the pattern now tha the
        shift and rotation are within tolerance. If we can't create a valid
        pattern we exit early.
        ------------------------------------------------------
        """
        # Lists where we will be putting our results.
        src_list = []
        ref_list = []

        # Counter for number of spokes we failed to find a reference
        # candidate for. We break the loop if we haven't found enough.
        n_fail = 0
        for src_idx in xrange(1, len(src_dist_array)):
            if n_fail > len(src_dist_array) - (n_match - 1):
                break

            # Given our length tolerance we can use it to compute a tolernace
            # on the angle between our spoke.
            src_sin_tol = (max_dist_rad /
                           (src_dist_array[src_idx] + max_dist_rad))
            # Test if the small angle approximation will still hold. This is
            # defined as when sin(theta) ~ theta to within 0.1%. This also
            # implicitly sets a minimum spoke length that we can use.
            if src_sin_tol > 0.0447:
                n_fail += 1
                continue

            # Find the reference pairs that include our candiates pattern
            # center and sort them in increasing delta
            ref_dist_idx_array = self._find_candidate_reference_pairs(
                src_dist_array[src_idx], ref_dist_array, max_dist_rad)

            # Test the spokes and return the id of the reference object.
            # Return none if no match is found.
            ref_id = self._test_spoke(
                src_ctr, src_delta_array[src_idx], src_dist_array[src_idx],
                src_delta_array[0], src_dist_array[0], ref_ctr, ref_ctr_id,
                ref_delta, ref_dist, ref_dist_idx_array, ref_delta_array,
                ref_dist_array,
                ref_id_array, src_sin_tol)
            if ref_id is None:
                n_fail += 1
                continue

            # Append the successful indicies to our list. The src_idx needs
            # an extra iteration to skip the first and second source objects.
            ref_list.append(ref_id)
            src_list.append(src_idx + 1)
            # If we found enough reference objects we can return early.
            if len(ref_list) >= n_match - 2:
                break
        return pipeBase.Struct(
            ref_spoke_list=ref_list,
            src_spoke_list=src_list,)

    # TODO: Clean up test_spoke and make input simpler.

    def _test_spoke(self, src_ctr, src_delta, src_dist, src_ctr_delta,
                    src_ctr_dist, ref_ctr, ref_ctr_id, ref_delta, ref_dist,
                    ref_dist_idx_array, ref_delta_array, ref_dist_array,
                    ref_id_array, src_sin_tol):
        """Test the rotation of our source spoke against the candidate
        reference spoke. This method makes heavy use of the small angle
        approximation assumption to perform the comparison.
        ------------------------------------------------------
        """

        # Precompute all of the source only calculations so we don't have
        # to do it for each interation in the reference loop.
        cos_theta_src = (np.dot(src_delta, src_ctr_delta) /
                         (src_dist * src_ctr_dist))
        cross_src = (np.cross(src_delta, src_ctr_delta) /
                     (src_dist * src_ctr_dist))
        dot_cross_src = np.dot(cross_src, src_ctr)

        # Loop over our candidate reference objects.
        for ref_dist_idx in ref_dist_idx_array:
            # Check the direction of the delta vector.
            ref_sign = 1
            if ref_id_array[ref_dist_idx, 1] == ref_ctr_id:
                ref_sign = -1

            # Compute the cos between our "center" reference vector and the
            # current reference candidate.
            cos_theta_ref = ref_sign*(
                np.dot(ref_delta_array[ref_dist_idx], ref_delta) /
                (ref_dist_array[ref_dist_idx] * ref_dist))

            # Make sure we can safely can safely make the comparison in case
            # our "center" and candidate vectors are mostly aligned.
            cos_comp = cos_theta_ref ** 2 < 1 - src_sin_tol**2
            # Test the difference of the reference angle against the
            # source angle. Assums that the delta between the two is
            # small.
            if cos_comp:
                if not ((cos_theta_src - cos_theta_ref)**2 /
                        (1 - cos_theta_ref**2) < src_sin_tol**2):
                    continue
            else:
                if not ((cos_theta_src - cos_theta_ref)**2 /
                        src_sin_tol**2 < src_sin_tol**2):
                    continue

            # The cosine doesn't tests the manitude of the angle but not
            # it's direction. To do that we need to know the sine as well.
            # This cross product calculation does that.
            cross_ref = ref_sign*(
                np.cross(ref_delta_array[ref_dist_idx], ref_delta) /
                (ref_dist_array[ref_dist_idx]*ref_dist))
            dot_cross_ref = np.dot(cross_ref, ref_ctr)

            # Check the value of the cos again to make sure that it is not
            # near zero.
            sin_comp = -src_sin_tol < cos_theta_src < src_sin_tol
            if sin_comp:
                if not (-src_sin_tol <
                        (dot_cross_src - dot_cross_ref) /
                        src_sin_tol < src_sin_tol):
                    continue
            else:
                if not (-src_sin_tol <
                        (dot_cross_src - dot_cross_ref) /
                        cos_theta_ref < src_sin_tol):
                    continue

            # Return the correct id depending on the direction of the vector.
            if ref_sign == 1:
                return ref_id_array[ref_dist_idx, 1]
            else:
                return ref_id_array[ref_dist_idx, 0]
        return None

    def _create_shift_rot_matrix(self, cos_rot_sq, shift_matrix, src_delta,
                                 ref_ctr, ref_delta):
        """ Create the final part of our spherical rotation matrix.
        ------------------------------------------------------
        """
        cos_rot = np.sqrt(cos_rot_sq)
        rot_src_delta = np.dot(shift_matrix, src_delta)
        delta_dot_cross = np.dot(np.cross(rot_src_delta, ref_delta), ref_ctr)

        sin_rot = np.sign(delta_dot_cross)*np.sqrt(1 - cos_rot_sq)
        ref_cross_matrix = np.array(
            [[0., -ref_ctr[2], ref_ctr[1]],
             [ref_ctr[2], 0., -ref_ctr[0]],
             [-ref_ctr[1], ref_ctr[0], 0.]], dtype=np.float64)
        rot_matrix = (
            cos_rot*np.identity(3) +
            sin_rot*ref_cross_matrix +
            (1. - cos_rot)*np.outer(ref_ctr, ref_ctr))

        shift_rot_matrix = np.dot(rot_matrix, shift_matrix)

        return pipeBase.Struct(
            sin_rot=sin_rot,
            shift_rot_matrix=shift_rot_matrix,)

    def _intermediate_verify(self, src_pattern, ref_pattern, shift_rot_matrix,
                             max_dist_rad):
        """ Perform an intermediate verify step. Rotate the matches references
        into the source frame and test their distances against tolerance. Only
        return true if all points are within tolerance.
        ------------------------------------------------------
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
        ------------------------------------------------------
        """

        print("Comparing pattern %i to previous %i rotations..." %
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

    def _match_sources(self, source_catalog, shift_rot_matrix, max_dist_rad):
        """ Shift both the reference and source catalog to their other frames
        finds their nearest neighbor using a kdTree. Removes all matches
        who do not agree when either the refernce or source catalog is rotated
        and removes all matches greated than the requested distance.
        ----------------------------------------------------------------------
        """
        shifted_references = np.dot(
            shift_rot_matrix.transpose(),
            self._reference_catalog.transpose()).transpose()
        shifted_sources = np.dot(
            shift_rot_matrix,
            source_catalog.transpose()).transpose()

        ref_matches = np.empty((len(shifted_references), 2),
                               dtype=np.int32)
        src_matches = np.empty((len(shifted_sources), 2),
                               dtype=np.int32)

        ref_matches[:, 1] = np.arange(len(shifted_references),
                                      dtype=np.int32)
        src_matches[:, 0] = np.arange(len(shifted_sources),
                                      dtype=np.int32)

        ref_kdtree = cKDTree(self._reference_catalog)
        src_kdtree = cKDTree(source_catalog)

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
        and reference objects agree they they are each other's
        nearist neighbor.
        ------------------------------------------------------
        """

        handshake_mask_array = np.zeros(len(matches_ref))
        for ref_match_idx, match in enumerate(matches_ref):
            src_match_idx = np.searchsorted(matches_src[:, 0], match[0])
            if match[1] == matches_src[src_match_idx, 1]:
                handshake_mask_array[ref_match_idx] = True
        return handshake_mask_array
