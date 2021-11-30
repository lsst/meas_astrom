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

#include <Eigen/Dense>
#include <ndarray.h>

namespace lsst {
namespace meas {
namespace astrom {

struct PatternResult {
    std::vector<std::pair<uint16_t, uint16_t>> candidate_pairs;
    std::vector<double> shift_rot_matrix;
    double cos_shift;
    double sin_rot;
    bool success = false;
};

struct RotationTestResult {
    double cos_rot_sq;
    Eigen::Vector3d proj_ref_ctr_delta;
    Eigen::Matrix<double, 3, 3> shift_matrix;
    bool success = false;
};

struct SortedArrayResult {
    std::vector<float> dists;
    std::vector<uint16_t> ids;
};

struct ShiftRotMatrixResult {
    double sin_rot;
    Eigen::Matrix<double, 3, 3> shift_rot_matrix;
};

/**
 * Test an input source pattern against the reference catalog.
 *
 * @param[in] src_pattern_array Sub selection of source 3 vectors to create a
 *     pattern from.
 * @param[in] src_delta_array Deltas of pairs from the source pattern center.
 *     (idx=0).
 * @param[in] src_dist_array Distances of the pairs in src_delta_array.
 * @param[in] dist_array Set of all distances between pairs of reference
 *     objects.
 * @param[in] id_array Pairs of indices creating the reference pairs. Sorted
 *     the same as dist_array.
 * @param[in] reference_array Set of all reference points.
 * @param[in] n_match Number of points to attempt to create a pattern from.
 *     Must be <= len(src_pattern_array)
 * @param[in] max_cos_theta_shift Maximum shift allowed between two patterns'
 *     centers. Equivalent to the shift on the sky.
 * @param[in] max_cos_rot_sq Maximum rotation between two patterns that have
 *     been shifted to have their centers on top of each other.
 * @param[in] max_dist_rad Maximum difference allowed between the source
 *     and reference pair distances to consider the reference pair a candidate
 *     for the source pair. Also sets the tolerance between the opening angles
 *     of the spokes when compared to the reference.
 * @return Struct containing: candidate pairs for this matched pattern; the
 *     rotation matrix to match the patterns, the implied rotations to create
 *     the matrix, a boolean representing success of the algorithm.
 */
PatternResult construct_pattern_and_shift_rot_matrix(
        ndarray::Array<double, 2, 1> src_pattern_array, ndarray::Array<double, 2, 1> src_delta_array,
        ndarray::Array<double, 1, 1> src_dist_array, ndarray::Array<float, 1, 1> dist_array,
        ndarray::Array<uint16_t, 2, 1> id_array, ndarray::Array<double, 2, 1> reference_array, int n_match,
        double max_cos_theta_shift, double max_cos_rot_sq, double max_dist_rad);

/**
 * Find the range of reference pairs within the distance tolerance of our
 * source pair spoke.
 *
 * Returns the min and max index spanning src_dist +/- max_dist_rad.
 *
 * @param[in] src_dist Value of the distance we would like to search for in
 *     the reference array in radians.
 * @param[in] ref_dist_array sorted array of distances in radians.
 * @param[in] max_dist_rad maximum plus/minus search to find in the reference
 *     array in radians.
 * @return pair of indices for the min and max range of indices to search.
 */
std::pair<size_t, size_t> find_candidate_reference_pair_range(
        float src_dist, ndarray::Array<float, 1, 1> const& ref_dist_array, double max_dist_rad);
/**
 * Find the range of reference pairs within the distance tolerance of our
 * source pair spoke.
 *
 * Returns the min and max index spanning src_dist +/- max_dist_rad.
 *
 * @param[in] src_dist float value of the distance we would like to search for
 *     in the reference array in radians.
 * @param[in] ref_dist_array sorted array of distances in radians.
 * @param[in] max_dist_rad maximum plus/minus search to find in the reference
 *     array in radians.
 * @return pair of indices for the min and max range of indices spanning +/-
 *     max_dist_rad.
 */
std::pair<size_t, size_t> find_candidate_reference_pair_range(float src_dist,
                                                              std::vector<float> const& ref_dist_array,
                                                              double max_dist_rad);

/**
 * Test if the rotation implied between the source pattern and reference
 * pattern is within tolerance. To test this we need to create the first part
 * of our spherical rotation matrix which we also return for use later.
 *
 * @param[in] src_center 3 vector defining the center of the candidate source
 *     pinwheel pattern.
 * @param[in] ref_center 3 vector defining the center of the candidate
 *     reference pinwheel pattern.
 * @param[in] src_delta 3 vector delta between the source pattern center and
 *     the end of the pinwheel spoke.
 * @param[in] ref_delta 3 vector delta of the candidate matched reference pair
 * @param[in] cos_shift Cosine of the angle between the source and reference
 *     candidate centers.
 * @param[in] max_cos_rot_sq Maximum allowed rotation of the pinwheel pattern.
 * @return Struct containing the magnitude of the rotation needed to align the
        two patterns after their centers are shifted on top of each and the
        rotation matrix describing the shift needed to align the source and
        candidate reference center.
 */
RotationTestResult test_rotation(ndarray::Array<double, 1, 1> const& src_center,
                                 ndarray::Array<double, 1, 1> const& ref_center,
                                 ndarray::Array<double, 1, 1> const& src_delta,
                                 ndarray::Array<double, 1, 1> const& ref_delta, double cos_shift,
                                 double max_cos_rot_sq);

/**
 * Create arrays sorted on the distance between the center of this
 * candidate reference object and the all reference objects.
 *
 * @param[in] ref_center Center point of the candidate pattern.
 * @param[in] reference_array Array of all reference object points.
 * @return Struct containing to vectors. First is the sorted distances and
 * second is the index locations of the pair that creates those distances in
 * reference_array.
 */
SortedArrayResult create_sorted_arrays(ndarray::Array<double, 1, 1> const& ref_center,
                                       ndarray::Array<double, 2, 1> const& reference_array);

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
        Eigen::Vector3d const& proj_ref_ctr_delta, std::vector<float> const& ref_dist_array,
        std::vector<uint16_t> const& ref_id_array, ndarray::Array<double, 2, 1> const& reference_array,
        double max_dist_rad, size_t n_match);

/**
 * Create the final part of our spherical rotation matrix.
 *
 * @param[in] cos_rot_sq cosine of the rotation needed to align our source and
 *     reference candidate patterns.
 * @param[in] shift_matrix 3x3 rotation matrix for shifting the source pattern
 *     center on top of the candidate reference pattern center.
 * @param[in] src_delta 3 vector delta representing the first spoke of the
 *     source pattern
 * @param[in] ref_ctr 3 vector on the unit-sphere representing the center of
 *     our reference pattern.
 * @param[in] ref_delta 3 vector delta made by the first pair of the reference
 *     pattern.
 * @return Struct containing constructed matrix and implied rotation.
 */
ShiftRotMatrixResult create_shift_rot_matrix(double cos_rot_sq,
                                             Eigen::Matrix<double, 3, 3> const& shift_matrix,
                                             ndarray::Array<double, 1, 1> const& src_delta,
                                             ndarray::Array<double, 1, 1> const& ref_ctr,
                                             ndarray::Array<double, 1, 1> const& ref_delta);

/**
 * Construct a generalized 3D rotation matrix about a given axis.
 *
 * param[in] rot_axis 3 vector defining the axis to rotate about.
 * param[in] cos_rotation cosine of the rotation angle.
 * param[in] sin_rotion sine of the rotation angle.
 * @return 3x3 spherical, rotation matrix.
 */
Eigen::Matrix<double, 3, 3> create_spherical_rotation_matrix(Eigen::Vector3d const& rot_axis,
                                                             double cos_rotation, double sin_rotion);

/**
 * Test the input rotation matrix against one input pattern and a second one.
 *
 * If every point in the pattern after rotation is within a distance of
 * max_dist_rad to its candidate point in the other pattern, we return True.
 *
 * @param[in] src_pattern Vector of 3 vectors representing the points that make
 *     up our source pinwheel pattern.
 * @param[in] ref_pattern Vector of 3 vectors representing our candidate
 *     reference pinwheel pattern.
 * @param[in] shift_rot_matrix 3x3 rotation matrix that takes the source
 *     objects and rotates them onto the frame of the reference objects.
 * @param[in] max_dist_rad Maximum distance allowed to consider two objects the
 *     same.
 * @param[in] n_match Number of points making up a pattern.
 * @return Comparison passes.
 */
bool intermediate_verify_comparison(std::vector<Eigen::Vector3d> const& src_pattern,
                                    std::vector<Eigen::Vector3d> const& ref_pattern,
                                    Eigen::Matrix<double, 3, 3> const& shift_rot_matrix, double max_dist_rad,
                                    int n_match);

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
                Eigen::Vector3d const& proj_ref_ctr_delta, double proj_ref_ctr_dist_sq,
                std::pair<size_t, size_t> const& candidate_range, std::vector<uint16_t> const& ref_id_array,
                ndarray::Array<double, 2, 1> const& reference_array, double src_sin_tol);

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
