// -*- LSST-C++ -*-

/*
 * This file is part of meas_astrom.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <ndarray.h>

namespace lsst {
namespace meas {
namespace astrom {
/**
 * Find the range of reference spokes within a spoke distance tolerance of our
 * source spoke.
 *
 * Returns an the min and max index spanning src_dist +/- max_dist_rad.
 *
 * @param[in] src_dist float value of the distance we would like to search for
 *     in the reference array in radians.
 * @param[in] ref_dist_array sorted array of distances in radians.
 * @param[in] max_dist_rad maximum plus/minus search to find in the reference
 *     array in radians.
 * @return pair of indices for the min and max range of indices to search.
 */
std::pair<size_t, size_t> find_candidate_reference_pair_range(
        float src_dist, ndarray::Array<float, 1, 1> const& ref_dist_array, float max_dist_rad);
/**
 * Create the individual spokes that make up the pattern now that the
 * shift and rotation are within tolerance.
 *
 * If we can't create a valid pattern we exit early with a partial result.
 *
 * @param[in] src_ctr 3 vector of the source pinwheel center
 * @param[in] src_delta_array Array of 3 vector deltas between the source
 *     center and the pairs that make up the remaining spokes of the pinwheel.
 * @param[in] src_dist_array Array of the distances of each src_delta in the
 *     pinwheel.
 * @param[in] ref_ctr 3 vector of the candidate reference center.
 * @param[in] proj_ref_ctr_delta Plane projected 3 vector formed from the
 *     center point of the candidate pin-wheel and the second point in the
 *     pattern to create the first spoke pair. This is the candidate pair that
 *     was matched in the main _construct_pattern_and_shift_rot_matrix loop.
 * @param[in] ref_dist_array Array of vector distances for each of the
 *     reference pairs
 * @param[in] ref_id_array Array of id lookups into the master reference array
 *     that our center id object is paired with.
 * @param[in] reference_array Full 3 vector data for the reference catalog.
 * @param[in] max_dist_rad Maximum search radius for distances.
 * @param[in] n_match Number of source deltas that must be matched into the
 *     reference deltas in order to consider this a successful pattern match.
 * @return Return pairs of reference ids and their matched src ids.
 */
std::vector<std::pair<size_t, size_t>> create_pattern_spokes(
        ndarray::Array<double, 1, 1> const& src_ctr, ndarray::Array<double, 2, 1> const& src_delta_array,
        ndarray::Array<double, 1, 1> const& src_dist_array, ndarray::Array<double, 1, 1> const& ref_ctr,
        ndarray::Array<double, 1, 1> const& proj_ref_ctr_delta,
        ndarray::Array<float, 1, 1> const& ref_dist_array, ndarray::Array<uint16_t, 1, 1> const& ref_id_array,
        ndarray::Array<double, 2, 1> const& reference_array, double max_dist_rad, size_t n_match);

/**
 * Check the opening angle between the first spoke of our pattern for the
 * source object against the reference object.
 *
 * This method makes heavy use of the small angle approximation to perform
 * the comparison.
 *
 * @param[in] cos_theta_src Cosine of the angle between the current candidate
 *     source spoke and the first spoke.
 * @param[in] sin_theta_src Sine of the angle between the current candidate
 *     source spoke and the first spoke.
 * @param[in] ref_ctr 3 vector of the candidate reference center
 * @param[in] proj_ref_ctr_delta Plane projected first spoke in the reference
 *     pattern using the pattern center as normal.
 * @param[in] proj_ref_ctr_dist_sq Squared length of the projected vector.
 * @param[in] candidate_range Min and max index locations in ref_id_array that
 *     have pair lengths within the tolerance range.
 * @param[in] ref_id_array Array of id lookups into the master reference array
 *     that our center id object is paired with. uint16 type is to limit the
 *     amount of memory used and is set in the pessimistic_pattern_matcher_3d
 *     python class with reference catalogs trimmed by the matchPessimisticB
 *     runner class.
 * @param[in] reference_array Array of three vectors representing the locations
 *     of all reference objects.
 * @param[in] src_sin_tol Sine of tolerance allowed between source and
 *     reference spoke opening angles.
 * @return ID of the candidate reference object successfully matched or -1 if
 *     no match is found.
 */
int check_spoke(double cos_theta_src, double sin_theta_src, ndarray::Array<double, 1, 1> const& ref_ctr,
                ndarray::Array<double, 1, 1> const& proj_ref_ctr_delta, double proj_ref_ctr_dist_sq,
                std::pair<size_t, size_t> const& candidate_range,
                ndarray::Array<uint16_t, 1, 1> const& ref_id_array,
                ndarray::Array<double, 2, 1> const& reference_array, double src_sin_tol);

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
