
from __future__ import division, print_function, absolute_import

from copy import copy

import numpy as np
from scipy.spatial import cKDTree

__deg_to_rad__ = np.pi/180.0


class OptimisticPatternMatcherB(object):
    """Class implimenting Optimistic Pattern Matcher B from Tabur 2007. The
    class loads and stores the reference object in a convienent data structure
    for matching any set of source obejcts that are assumed to contain each
    other.
    ----------------------------------------------------------------------------
    Attributes:
        reference_catalog: Input array of x, y, mag of reference objects.
        max_rotation_theta: The max "shift" distance allowed in degrees.
        max_rotation_phi: The max rotation allowed in degrees
        dist_tol: float epislon distance to consider a pair of objects the
        same. Units are in degrees.
        max_dist_cand: int max number of candidates to test For dense fields
        even a small dist_tol can lead to looping over a large number of
        candidates. This cuts off that loop early so we don't slow the code by
        looking over too many candidates.
        ang_tol: float max tolerance to consider the angles between two
        spokes in the pattern the same. Units are degrees.
        max_match_dist: float Maximum distance after shift and rotation are
        applied to consider two objects a match in the KDTree.
        min_matches: int minimum number of objects required to be matched to
        consider the match valid.
        max_n_patterns: int Number of patterns to attempt to create before
        exiting with a failure.
    """

    def __init__(self, reference_catalog, max_rotation_theta,
                 max_rotation_phi, dist_tol, max_dist_cand, ang_tol,
                 max_match_dist, min_matches, max_n_patterns):
        self._reference_catalog = copy(reference_catalog[:, :3])
        self._n_reference = len(self._reference_catalog)

        self._max_cos_theta = np.cos(max_rotation_theta*__deg_to_rad__)
        self._max_cos_phi_sq = np.cos(max_rotation_phi*__deg_to_rad__)**2

        self._dist_tol = dist_tol*__deg_to_rad__
        self._max_dist_cand = max_dist_cand
        self._ang_tol = ang_tol*__deg_to_rad__

        self._max_match_dist = max_match_dist*__deg_to_rad__
        self._min_matches = min_matches
        self._max_n_patterns = max_n_patterns
        # These two tests set limits the values of cosine theta in the spoken
        # pattern part of the pattern matcher. These avoid divide by zero
        # and also enforce a small angle approximation for the spokes at
        # 0, 180 degrees of the first pair and 90, 270 respectively.
        self._cos_limit = np.cos(np.max((1.4e-8, self._ang_tol)))
        self._sin_limit = np.sin(np.max((1.4e-8, self._ang_tol)))

        self._is_valid_rotation = False

        self._build_distances_and_angles()

    def _build_distances_and_angles(self):
        """Internal function for constructing for searchable distances and
        angles between pairs of objects in the reference catalog.
        """
        # Build the empty arrays we will store the distances, delta vectors,
        # and pair_ids into.
        self._id_array = np.empty(
            (int(self._n_reference*(self._n_reference - 1)/2), 2),
            dtype=np.int)
        self._dist_array = np.empty(
            int(self._n_reference*(self._n_reference - 1)/2),
            dtype=np.float64)
        self._delta_array = np.empty(
            (int(self._n_reference*(self._n_reference - 1)/2), 3),
            dtype=np.float64)
        pair_idx_array = np.empty(
            (self._n_reference, self._n_reference - 1),
            dtype=np.int)
        # Start the loop over the n choose 2 pairs.
        start_idx = 0
        start_idx_list = []
        for ref_idx, ref_obj in enumerate(self._reference_catalog):
            end_idx = self._n_reference - 1 - ref_idx
            self._id_array[start_idx: start_idx + end_idx, 0] = ref_idx
            self._id_array[start_idx: start_idx + end_idx, 1] = np.arange(
                ref_idx + 1, self._n_reference, dtype=np.int_)
            self._delta_array[start_idx: start_idx + end_idx, 0] = (
                self._reference_catalog[ref_idx + 1:, 0] - ref_obj[0])
            self._delta_array[start_idx: start_idx + end_idx, 1] = (
                self._reference_catalog[ref_idx + 1:, 1] - ref_obj[1])
            self._delta_array[start_idx: start_idx + end_idx, 2] = (
                self._reference_catalog[ref_idx + 1:, 2] - ref_obj[2])
            self._dist_array[start_idx: start_idx + end_idx] = (
                self._delta_array[start_idx: start_idx + end_idx, 0]**2 +
                self._delta_array[start_idx: start_idx + end_idx, 1]**2 +
                self._delta_array[start_idx: start_idx + end_idx, 2]**2)
            tmp_pair_idx_second = []
            for shift_idx, tmp_start_idx in enumerate(start_idx_list):
                tmp_pair_idx_second.append(ref_idx - 1 - shift_idx +
                                           tmp_start_idx)
            tmp_pair_idx_first = np.arange(start_idx, start_idx + end_idx,
                                           dtype=np.int)
            pair_idx_array[ref_idx, :] = np.concatenate(
                (np.array(tmp_pair_idx_second), tmp_pair_idx_first))
            start_idx_list.append(start_idx)
            start_idx += end_idx
        # Sort our arrays by the distance of the pair.
        self._dist_array = np.sqrt(self._dist_array)
        tmp_dist_map_array = self._dist_array[pair_idx_array]
        self._sorted_args = self._dist_array.argsort()
        self._dist_array = self._dist_array[self._sorted_args]
        self._pair_idx_array = np.empty_like(pair_idx_array)
        for ref_idx, tmp_dist_array in enumerate(tmp_dist_map_array):
            tmp_sort_dist_array = np.sort(tmp_dist_array)
            self._pair_idx_array[ref_idx, :] = np.searchsorted(
                self._dist_array, tmp_sort_dist_array)
        self._id_array = self._id_array[self._sorted_args]
        self._delta_array = self._delta_array[self._sorted_args]

        return None

    def _candidate_sort(self, dist_array, cand_dist, start_idx):
        """Internal function for sorting an array by distance relative to the
        candidate distance out.
        """
        return np.argsort(np.fabs(dist_array - cand_dist)) + start_idx

    def _construct_and_match_pattern(self, source_candidates, n_match):
        """Given a list of source canidates we check the pinwheel pattern they
        create against the reference catalog by checking distances an angles.
        We keep the calculations as simple and fast as posible by using mostly
        # cosine and sine relations.
        """
        matched_references = []
        # Create our vector and distances for the source object pinwheel.
        source_delta = np.empty((len(source_candidates) - 1, 3))
        source_delta[:, 0] = source_candidates[1:, 0] - source_candidates[0, 0]
        source_delta[:, 1] = source_candidates[1:, 1] - source_candidates[0, 1]
        source_delta[:, 2] = source_candidates[1:, 2] - source_candidates[0, 2]
        source_dist_array = np.sqrt(source_delta[:, 0]**2 +
                                    source_delta[:, 1]**2 +
                                    source_delta[:, 2]**2)
        # We first test if the distance of the first (AB) spoke of our source
        # pinwheel can be found in the array of reference catalog pairs.
        start_idx = np.searchsorted(
            self._dist_array,
            source_dist_array[0] - self._dist_tol)
        end_idx = np.searchsorted(
            self._dist_array,
            source_dist_array[0] + self._dist_tol,
            side='right')
        # If we couldn't find any candidate references distances we exit. We
        # also test if the edges to make sure we are not running over the array
        # size.
        if start_idx == end_idx:
            return ([], None, None)
        if start_idx < 0:
            start_idx = 0
        if end_idx > self._dist_array.shape[0]:
            end_idx = self._dist_array.shape[0]
        # Now that we have candiates reference distances for the first spoke of
        # the pinwheel we loop over them and attempt to construct the rest of
        # the pinwheel. We loop over them from the smallest difference in
        # distane between the reference and source to the largest.
        cand_idx_array = self._candidate_sort(
            self._dist_array[start_idx:end_idx], source_dist_array[0],
            start_idx)
        # if cand_idx_array.shape[0] > self._max_dist_cand:
            # print("Pattern cand array greater than max_dist_cand...")
            # print("\tMax dist diff will be %.4f arcsec" %
            #       ((source_dist_array[0] -
            #         self._dist_array[
            #             cand_idx_array[self._max_dist_cand - 1]]) /
            #        __deg_to_rad__*3600))
        for dist_idx in cand_idx_array[:self._max_dist_cand]:
            # Compute the value of cos_theta between our source candidate and
            # both reference objects. We will pick the one with the smaller
            # difference in angle, larger cosine to our candicate source.
            tmp_cos_theta_tuple = (
                np.dot(source_candidates[0],
                       self._reference_catalog[self._id_array[dist_idx][0]]),
                np.dot(source_candidates[0],
                       self._reference_catalog[self._id_array[dist_idx][1]]))
            for cos_idx, cos_theta in enumerate(tmp_cos_theta_tuple):
                matched_references = []
                # Now we can test the displacement we have on the sky between
                # the centers of our source and reference pinwheels, exiting if
                # it is to distant.
                if cos_theta < self._max_cos_theta:
                    continue
                # Store the candidate reference center that matches our source
                # in the correct order.
                if cos_idx == 0:
                    matched_references.append(self._id_array[dist_idx][0])
                    matched_references.append(self._id_array[dist_idx][1])
                    ref_delta = self._delta_array[dist_idx]
                else:
                    matched_references.append(self._id_array[dist_idx][1])
                    matched_references.append(self._id_array[dist_idx][0])
                    ref_delta = -self._delta_array[dist_idx]
                # Store the distance of this reference delta vector to speed up
                # later computation.
                ref_delta_dist = self._dist_array[dist_idx]
                # Since we already have the first two reference candidate ids
                # we can narrow our search to only those pairs that contain our
                # pinwheel reference and exclude the reference we have already
                # used to match the first spoke.
                # id_mask = np.logical_or(
                #     np.logical_and(
                #         self._id_array[:, 0] == matched_references[0],
                #         self._id_array[:, 1] != matched_references[1]),
                #     np.logical_and(
                #         self._id_array[:, 1] == matched_references[0],
                #         self._id_array[:, 0] != matched_references[1]))
                # tmp_ref_dist_arary = self._dist_array[id_mask]
                # tmp_ref_delta_array = self._delta_array[id_mask]
                # tmp_ref_id_array = self._id_array[id_mask]
                tmp_ref_dist_arary = self._dist_array[
                    self._pair_idx_array[matched_references[0]]]
                tmp_ref_delta_array = self._delta_array[
                    self._pair_idx_array[matched_references[0]]]
                tmp_ref_id_array = self._id_array[
                    self._pair_idx_array[matched_references[0]]]
                # Now we can start our loop to look for the remaining candidate
                # spokes of our pinwheel.
                n_failed = 0
                for cand_idx in xrange(1, len(source_dist_array)):
                    # Our test on the individual spokes.
                    match = self._pattern_spoke_test(
                        source_dist_array[cand_idx], source_delta[cand_idx],
                        source_candidates[0], source_delta[0],
                        source_dist_array[0], matched_references[0], ref_delta,
                        ref_delta_dist, tmp_ref_dist_arary, tmp_ref_delta_array,
                        tmp_ref_id_array)
                    # If we don't find a mach for this spoke we can exit early.
                    if match is None:
                        n_failed += 1
                        if n_failed >= len(source_candidates) - n_match:
                            break
                        continue
                    matched_references.append(match)
                    # If we've found enough spokes that agree with the referce
                    # pattern we can exit early.
                    if len(matched_references) >= n_match:
                        break
                # If if we've found a match for each spoke we can exit early
                # and then return the matches. We can also send off the
                # rotations we have already computed.
                if len(matched_references) >= n_match:
                    self._construct_rotation_matricies(
                       source_candidates[0],
                       self._reference_catalog[matched_references[0]],
                       source_delta[0], ref_delta, cos_theta)
                    if self._is_valid_rotation:
                        break
            if self._is_valid_rotation:
                break
        # Return the matches. If found.
        if len(matched_references) >= n_match:
            return matched_references
        return []

    def _pattern_spoke_test(self, cand_dist, cand_delta, source_center,
                            source_delta, source_delta_dist, ref_center_id,
                            ref_delta, ref_delta_dist, ref_dist_array,
                            ref_delta_array, ref_id_array):
        """Internal function finding matches for the remaining spokes of our
        candidate pinwheel.
        """
        # As before we first check references with matching distances, exiting
        # early if we find none.
        start_idx = np.searchsorted(
            ref_dist_array, cand_dist - self._dist_tol)
        end_idx = np.searchsorted(
            ref_dist_array, cand_dist + self._dist_tol,
            side='right')
        if start_idx == end_idx:
            return None
        if start_idx < 0:
            start_idx = 0
        if end_idx > ref_dist_array.shape[0]:
            end_idx = ref_dist_array.shape[0]
        # Loop over the posible matches and test them for quality.
        hold_id = -99
        cand_idx_array = self._candidate_sort(
            ref_dist_array[start_idx:end_idx], cand_dist, start_idx)
        # if cand_idx_array.shape[0] > self._max_dist_cand:
        #     print("Spoke cand array greater than max_dist_cand...")
        #     print("\tMax dist diff will be %.4f arcsec" %
        #           ((cand_dist -
        #             ref_dist_array[
        #                 cand_idx_array[self._max_dist_cand - 1]]) /
        #            __deg_to_rad__*3600))
        for dist_idx in cand_idx_array[:self._max_dist_cand]:
            # First we compute the dot product between our delta
            # vectors in each of the source and reference pinwheels
            # and test that they are the same within tolerance. Since
            # we know the distances of these deltas already we can
            # normalize the vectors and compare the values of cos_theta
            # between the two legs.
            ref_sign = 1
            if ref_id_array[dist_idx, 1] == ref_center_id:
                ref_sign = -1
            cos_theta_source = (np.dot(cand_delta, source_delta) /
                                (cand_dist*source_delta_dist))
            cos_theta_ref = ref_sign*(
                np.dot(ref_delta_array[dist_idx], ref_delta) /
                (ref_dist_array[dist_idx]*ref_delta_dist))
            # We need to test that the vectors are not completely aligned.
            # If they are our first test will be invalid thanks to
            # 1 - cos_theta**2 equaling zero.
            # Using a few trig relations and taylor expantions around
            # _ang_tol we compare the opening angles of our pinwheel
            # legs to see if they are within tolerance.
            cos_comparison = -self._cos_limit < cos_theta_ref < self._cos_limit
            if (cos_comparison and
                not ((cos_theta_source - cos_theta_ref)**2 /
                     (1 - cos_theta_ref**2) < self._ang_tol**2)):
                continue
            elif (not cos_comparison and
                  not ((cos_theta_source - cos_theta_ref)**2 /
                       self._sin_limit**2 < self._ang_tol**2)):
                continue
            # Now we compute the cross product between the first
            # rungs of our spokes and our candidate rungs. We then
            # dot these into our center vector to make sure the
            # rotation direction and amount of rotation are correct.
            # If they are not we move on to the next candidate.
            cross_source = (np.cross(cand_delta, source_delta) /
                            (cand_dist*source_delta_dist))
            cross_ref = ref_sign*(
                np.cross(ref_delta_array[dist_idx], ref_delta) /
                (ref_dist_array[dist_idx]*ref_delta_dist))
            dot_cross_source = np.dot(cross_source, source_center)
            dot_cross_ref = np.dot(cross_ref,
                                   self._reference_catalog[ref_center_id])
            # The test below tests for both the corality of the rotation
            # and its amount. We assume that the cross product is aligned
            # with the centeral vectors. Again using trig relations and
            # small angle aproximation on _ang_tol we arrive at the
            # folloing relation.
            sin_comparison = -self._sin_limit < cos_theta_ref < self._sin_limit
            if (not sin_comparison and
                not (-self._ang_tol <
                     (dot_cross_source - dot_cross_ref)/cos_theta_ref <
                     self._ang_tol)):
                continue
            elif (sin_comparison and
                  not (-self._ang_tol <
                       (dot_cross_source - dot_cross_ref)/self._sin_limit <
                       self._ang_tol)):
                continue
            # Check to see which id we should return.
            if ref_sign == 1:
                hold_id = ref_id_array[dist_idx, 1]
            else:
                hold_id = ref_id_array[dist_idx, 0]
            break
        # Return the id of our matched object that makes up this spoke if we
        # found it.
        if hold_id >= 0:
            return hold_id
        return None

    def _construct_rotation_matricies(self, source_candidate, ref_candidate,
                                      source_delta, ref_delta, cos_theta):
        self._is_valid_rotation = False
        # First we compute the unit vector for the axis of rotation between
        # our two candidate centers. This gives us the overal shift.
        if cos_theta > 1.0:
            cos_theta = 1.
        elif cos_theta < -1.0:
            cos_theta = -1.
        sin_theta = np.sqrt(1. - cos_theta**2)
        # We need to test that we actually have to do this rotation. If the
        # vectors are already aligned we can skip the first rotation and just
        # store the identidity.
        if sin_theta != 0:
            rot_axis = np.cross(source_candidate, ref_candidate)
            rot_axis /= sin_theta
            # Now that we have our axis and cos_theta from before we can rotate
            # about it to align the source and candidate vectors. This is our
            # first rotation matrix.
            rot_cross_matrix = np.array(
                [[0., -rot_axis[2], rot_axis[1]],
                 [rot_axis[2], 0., -rot_axis[0]],
                 [-rot_axis[1], rot_axis[0], 0.]], dtype=np.float64)
            self.theta_rot_matrix = (
                cos_theta*np.identity(3) +
                sin_theta*rot_cross_matrix +
                (1. - cos_theta)*np.outer(rot_axis, rot_axis))
        else:
            self.theta_rot_matrix = np.identity(3)
        # Now we rotate our source delta to the frame of the reference.
        rot_source_delta = np.dot(self.theta_rot_matrix, source_delta)
        # Now that we've rotated our source_delta vector we can dot it
        # into the reference delta vector and test that the amount of
        # rotation that would align the two is within our rotation
        # tolerance.
        cos_phi_sq = (np.dot(rot_source_delta, ref_delta)**2 /
                      (np.dot(rot_source_delta, rot_source_delta) *
                       np.dot(ref_delta, ref_delta)))
        if cos_phi_sq < self._max_cos_phi_sq:
            self._is_valid_rotation = False
            return None
        # Now that we know the rotation is valid we can compute the cos_phi
        # and sin_phi that we need to compute our rotation matrix.
        cos_phi = np.sqrt(cos_phi_sq)
        delta_dot_cross = np.dot(np.cross(rot_source_delta, ref_delta),
                                 ref_candidate)
        sin_phi = np.sign(delta_dot_cross)*np.sqrt(1 - cos_phi_sq)
        # Construct our phi rotation matrix.
        ref_cross_matrix = np.array(
            [[0., -ref_candidate[2], ref_candidate[1]],
             [ref_candidate[2], 0., -ref_candidate[0]],
             [-ref_candidate[1], ref_candidate[0], 0.]], dtype=np.float64)
        self.phi_rot_matrix = (
            cos_phi*np.identity(3) +
            sin_phi*ref_cross_matrix +
            (1. - cos_phi)*np.outer(ref_candidate, ref_candidate))
        # Now that we have our two rotations around theta and phi we can dot
        # the two matrices together to create our full rotation matrix.
        self.rot_matrix = np.dot(self.phi_rot_matrix, self.theta_rot_matrix)
        # Store cos_theta and sin_phi for our rotation matrix for later use.
        self._cos_theta = cos_theta
        self._sin_phi = sin_phi
        # This is a valid rotation.
        self._is_valid_rotation = True
        return None

    def _compute_shift_and_match_sources(self, source_catalog):
        """Given an input source catalog, pinwheel centers in the source and
        reference catalog, and a cosine and sine rotation return a shifted
        catalog for matching.
        """
        # Since the sources are the kd-tree objects we rotate the
        # references.
        shifted_references = np.dot(
            self.rot_matrix.transpose(),
            self._reference_catalog[:, :3].transpose()).transpose()
        # Empty arrays for output.
        output_matches = np.empty((len(shifted_references), 2),
                                  dtype=np.int_)
        # Store the "id" of the references
        output_matches[:, 1] = np.arange(len(shifted_references),
                                         dtype=np.int_)
        # Find the matches.
        tmp_src_dist, tmp_src_idx = self._kdtree.query(
            shifted_references[:, :3])
        output_matches[:, 0] = tmp_src_idx
        tmp_src_dist = tmp_src_dist
        unique_mask = self._test_unique_matches(output_matches[:, 0],
                                                tmp_src_dist)
        # Mask on the max distance and return the masked arrays.
        dist_mask = np.logical_and(unique_mask,
                                   tmp_src_dist < self._max_match_dist)
        return output_matches[dist_mask], tmp_src_dist[dist_mask]

    def _test_unique_matches(self, idx_array, dist_array):
        """Internal function for creating a mask of matches with unique indices.
        """
        # TODO:
        #     Only returns matches that are truly unique. Could later add test
        # on distances to return the closest match only.
        unique_array, unique_idx_array, inverse_array = np.unique(
            idx_array, return_index=True, return_inverse=True)
        unique_mask = np.zeros(len(idx_array), dtype=np.bool)
        unique_mask[unique_idx_array] = True
        used_non_unique = np.array([], np.int64)
        for non_unique_idx in idx_array[inverse_array[unique_array.shape[0]:]]:
            if np.any(non_unique_idx == used_non_unique):
                continue
            used_non_unique = np.concatenate((used_non_unique,
                                              [non_unique_idx]))
            dub_idx_array = np.argwhere(idx_array == non_unique_idx)
            unique_mask[
                dub_idx_array[np.argmin(dist_array[dub_idx_array])]] = True
        return unique_mask

    def match(self, source_catalog, n_check, n_match):
        """Function for matching a given source catalog into the loaded
        reference catalog.
        ----------------------------------------------------------------------
        Args:
            source_catalog: float array of spherical x,y,z coordinates and a
                magnitude.
            n_check: int value specifying the number of sources to attempt a
                match on. Not all may be checked if n_match criteria is met
                before hand. n_check should be greater than n_match by 1-3
                objects.
            n_match: Number of objects to use in constructing a pattern.
        Returns:
            tuple (2D int array of matched pairs,
                   float array of pair distances)
        """
        # Given our input source_catalog we sort on magnitude.
        sorted_catalog = source_catalog[source_catalog[:, -1].argsort()]
        n_source = len(sorted_catalog)
        # If there are more sources we store them in the kd-tree. Opposite
        # if there are more references. This we we will always have unique
        # matches.
        self._kdtree = cKDTree(source_catalog[:, :3])
        # Loop through the sources from brightest to faintest grabbing a chucnk
        # of n_check each time.
        for pattern_idx in xrange(np.min((self._max_n_patterns,
                                          n_source - n_check))):
            matches = []
            distances = []
            # Grab the sources to attempt to create this pattern.
            pattern = sorted_catalog[pattern_idx: pattern_idx + n_check, :3]
            # Construct a pattern given the number of points we are using to
            # create it.
            ref_candidates = self._construct_and_match_pattern(pattern,
                                                               n_match)
            # If we have enough candidates we can shift and attempt to match
            # the two catalogs.
            if len(ref_candidates) >= n_match:
                print('Matching...')
                matches, distances = self._compute_shift_and_match_sources(
                    source_catalog)
                print('Matches:', len(matches))
                # If the number of matched objects satifies our criteria we
                # can print summary statistics and exit. If not we start the
                # loop over with the next pattern.
                if len(matches) > self._min_matches:
                    print("Succeeded after %i patterns." % pattern_idx)
                    print("\tShift %.4f arcsec" %
                          (np.arccos(self._cos_theta)*3600/__deg_to_rad__))
                    print("\tRotation: %.4f deg" %
                          (np.arcsin(self._sin_phi)/__deg_to_rad__))
                    break
        if len(matches) < n_match:
            print("Failed after %i patterns." % pattern_idx)
            return ([], [])
        return (matches, distances)
