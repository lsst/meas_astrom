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

#include "stdio.h"
#include <cmath>
#include "ndarray/eigen.h"
#include "lsst/meas/astrom/pessimisticPatternMatcherUtils.h"

namespace {

/// Return -1, 0, or 1, depending on whether val is negative, zero, or positive.
int sgn(double val) { return (0.0 < val) - (val < 0.0); }
}  // namespace

namespace lsst {
namespace meas {
namespace astrom {

PatternResult construct_pattern_and_shift_rot_matrix(
        ndarray::Array<double, 2, 1> src_pattern_array, ndarray::Array<double, 2, 1> src_delta_array,
        ndarray::Array<double, 1, 1> src_dist_array, ndarray::Array<float, 1, 1> dist_array,
        ndarray::Array<uint16_t, 2, 1> id_array, ndarray::Array<double, 2, 1> reference_array, size_t n_match,
        double max_cos_theta_shift, double max_cos_rot_sq, double max_dist_rad) {
    // Our first test. We search the reference dataset for pairs that have the same length as our first source
    // pairs to within plus/minus the max_dist tolerance.
    std::pair<size_t, size_t> candidate_range =
            find_candidate_reference_pair_range(src_dist_array[0], dist_array, max_dist_rad);
    size_t ref_dist_idx = (candidate_range.first + candidate_range.second) / 2;

    // Start our loop over the candidate reference objects. Looping from the inside (minimum difference to our
    // source dist) to the outside.
    for (size_t idx = 0; idx < candidate_range.second - candidate_range.first; idx++) {
        // TODO DM-33514: cleanup this loop to use an iterator that handles the "inside-out" iteration.
        if (idx % 2 == 0) {
            ref_dist_idx = ref_dist_idx + idx;
        } else {
            ref_dist_idx = ref_dist_idx - idx;
        }
        // We have two candidates for which reference object corresponds with the source at the center of our
        // pattern. As such we loop over and test both possibilities.
        ndarray::Array<uint16_t, 1, 1> tmp_ref_pair_list = id_array[ref_dist_idx];
        for (uint16_t ref_pair_idx = 0; ref_pair_idx < 2; ref_pair_idx++) {
            uint16_t ref_id = id_array[ref_dist_idx][ref_pair_idx];
            std::vector<std::pair<uint16_t, uint16_t>> candidate_pairs;
            // Test the angle between our candidate ref center and the source center of our pattern. This
            // angular distance also defines the shift we will later use.
            ndarray::Array<double, 1, 1> ref_center = reference_array[ref_id];
            double cos_shift =
                    ndarray::asEigenMatrix(src_pattern_array[0]).dot(ndarray::asEigenMatrix(ref_center));
            fprintf(stdout, "ref_center = %.20f, %.20f\n", ref_center[0], ref_center[1]);
            fprintf(stdout, "src_pattern = %.20f, %.20f\n", src_pattern_array[0][0], src_pattern_array[0][1]);
            fprintf(stdout, "cos_shift = %.20f, %.20f\n", cos_shift, max_cos_theta_shift);
            if (cos_shift < max_cos_theta_shift) {
                continue;
            }
            // We can now append this one as a candidate.
            candidate_pairs.push_back(std::make_pair(ref_id, 0));
            ndarray::Array<double, 1, 1> ref_delta;
            // Test to see which reference object to use in the pair.
            if (ref_id == *tmp_ref_pair_list.begin()) {
                candidate_pairs.push_back(std::make_pair(tmp_ref_pair_list[1], 1));
                ref_delta = copy(reference_array[tmp_ref_pair_list[1]] - ref_center);
            } else {
                candidate_pairs.push_back(std::make_pair(tmp_ref_pair_list[0], 1));
                ref_delta = copy(reference_array[tmp_ref_pair_list[0]] - ref_center);
            }
            // For dense fields it will be faster to compute the absolute rotation this pair suggests first
            // rather than saving it after all the spokes are found. We then compute the cos^2 of the rotation
            // and first part of the rotation matrix from source to reference frame.
            RotationTestResult test_rot_result =
                    test_rotation(src_pattern_array[0], ref_center, src_delta_array[0], ref_delta, cos_shift,
                                  max_cos_rot_sq);
            if (!test_rot_result.success) {
                continue;
            }
            // Now that we have a candidate first spoke and reference pattern center, we mask our future
            // search to only those pairs that contain our candidate reference center.
            SortedArrayResult sorted_array_struct = create_sorted_arrays(ref_center, reference_array);
            // Now we feed this sub data to match the spokes of our pattern.
            std::vector<std::pair<size_t, size_t>> pattern_spokes =
                    create_pattern_spokes(src_pattern_array[0], src_delta_array, src_dist_array, ref_center,
                                          test_rot_result.proj_ref_ctr_delta, sorted_array_struct.dists,
                                          sorted_array_struct.ids, reference_array, max_dist_rad, n_match);
            // If we don't find enough candidates we can continue to the next reference center pair.
            if (pattern_spokes.size() < n_match - 2) {
                continue;
            }
            // If we have the right number of matched ids we store these.
            candidate_pairs.reserve(candidate_pairs.size() + pattern_spokes.size());
            candidate_pairs.insert(candidate_pairs.end(), pattern_spokes.begin(), pattern_spokes.end());
            // We can now create our full matrix for both the shift and rotation. The shift aligns
            // the pattern centers, while the rotation rotates the spokes on top of each other.
            ShiftRotMatrixResult shift_rot_result =
                    create_shift_rot_matrix(test_rot_result.cos_rot_sq, test_rot_result.shift_matrix,
                                            src_delta_array[0], ref_center, ref_delta);

            // concatenate our final patterns.
            Eigen::Matrix3d shift_rot_matrix = shift_rot_result.shift_rot_matrix;
            std::vector<Eigen::Vector3d> ref_pattern;
            std::vector<Eigen::Vector3d> src_pattern;
            ref_pattern.reserve(n_match);
            src_pattern.reserve(n_match);
            for (size_t idx = 0; idx < n_match; idx++) {
                ref_pattern[idx] = ndarray::asEigenMatrix(reference_array[candidate_pairs[idx].first]);
                src_pattern[idx] = ndarray::asEigenMatrix(src_pattern_array[candidate_pairs[idx].second]);
            }
            // TODO: DM-32985 Implement least squares linear fitting here. Look at python example linked in
            // ticket for how to implement.
            // Test that the all points in each pattern are within tolerance after shift/rotation.
            bool passed =
                    intermediate_verify_comparison(src_pattern, ref_pattern, shift_rot_matrix, max_dist_rad);
            if (passed) {
                return PatternResult(candidate_pairs, shift_rot_matrix, cos_shift, shift_rot_result.sin_rot);
            }
        }
    }
    // failed fit
    return PatternResult();
}

SortedArrayResult create_sorted_arrays(ndarray::Array<double, 1, 1> const& ref_center,
                                       ndarray::Array<double, 2, 1> const& reference_array) {
    SortedArrayResult result;
    // NOTE: this algorithm is quadratic in the length of reference_array. It might be worth using std::sort
    // instead of this approach, if reference_array is long, but the length at which that matters should
    // be tested because it would require caching of distances.
    for (uint16_t idx = 0; idx < reference_array.getShape()[0]; idx++) {
        Eigen::Vector3d diff = ndarray::asEigenMatrix(copy(reference_array[idx] - ref_center));
        double dist = sqrt(diff.dot(diff));
        auto dists_itr = std::lower_bound(result.dists.begin(), result.dists.end(), dist);
        auto ids_itr = result.ids.begin() + (dists_itr - result.dists.begin());
        result.ids.insert(ids_itr, idx);
        result.dists.insert(dists_itr, dist);
    }
    return result;
}

std::pair<size_t, size_t> find_candidate_reference_pair_range(
        float src_dist, ndarray::Array<float, 1, 1> const& ref_dist_array, double max_dist_rad) {
    auto itr =
            std::lower_bound(ref_dist_array.begin(), ref_dist_array.end(), src_dist - max_dist_rad - 1e-16);
    auto itrEnd =
            std::upper_bound(ref_dist_array.begin(), ref_dist_array.end(), src_dist + max_dist_rad + 1e-16);

    size_t startIdx = itr - ref_dist_array.begin();
    size_t endIdx = itrEnd - ref_dist_array.begin();
    return std::make_pair(startIdx, endIdx);
}

std::pair<size_t, size_t> find_candidate_reference_pair_range(float src_dist,
                                                              std::vector<float> const& ref_dist_array,
                                                              double max_dist_rad) {
    auto itr =
            std::lower_bound(ref_dist_array.begin(), ref_dist_array.end(), src_dist - max_dist_rad - 1e-16);
    auto itrEnd =
            std::upper_bound(ref_dist_array.begin(), ref_dist_array.end(), src_dist + max_dist_rad + 1e-16);

    size_t startIdx = itr - ref_dist_array.begin();
    size_t endIdx = itrEnd - ref_dist_array.begin();
    return std::make_pair(startIdx, endIdx);
}

RotationTestResult test_rotation(ndarray::Array<double, 1, 1> const& src_center,
                                 ndarray::Array<double, 1, 1> const& ref_center,
                                 ndarray::Array<double, 1, 1> const& src_delta,
                                 ndarray::Array<double, 1, 1> const& ref_delta, double cos_shift,
                                 double max_cos_rot_sq) {
    // Make sure the sine is a real number.
    if (cos_shift > 1.0) {
        cos_shift = 1;
    } else if (cos_shift < -1.0) {
        cos_shift = -1;
    }
    double sin_shift = sqrt(1 - cos_shift * cos_shift);

    // If the sine of our shift is zero we only need to use the identity matrix for the shift. Else we
    // construct the rotation matrix for shift.
    auto ref_center_eigen = ndarray::asEigenMatrix(ref_center).head<3>();
    Eigen::Matrix3d shift_matrix;
    if (sin_shift > 0) {
        Eigen::Vector3d rot_axis = ndarray::asEigenMatrix(src_center).head<3>().cross(ref_center_eigen);
        rot_axis /= sin_shift;
        shift_matrix = create_spherical_rotation_matrix(rot_axis, cos_shift, sin_shift);
    } else {
        shift_matrix = Eigen::Matrix3d::Identity();
    }

    // Now that we have our shift we apply it to the src delta vector and check the rotation.
    Eigen::Vector3d rot_src_delta = shift_matrix * ndarray::asEigenMatrix(src_delta);
    Eigen::Vector3d proj_src_delta = rot_src_delta - rot_src_delta.dot(ref_center_eigen) * ref_center_eigen;
    Eigen::Vector3d ref_delta_eigen = ndarray::asEigenMatrix(ref_delta);

    Eigen::Vector3d proj_ref_delta =
            ref_delta_eigen - ref_delta_eigen.dot(ref_center_eigen) * ref_center_eigen;
    double proj_src_delta_sq = proj_src_delta.dot(proj_ref_delta);
    double cos_rot_sq = proj_src_delta_sq * proj_src_delta_sq /
                        (proj_src_delta.dot(proj_src_delta) * proj_ref_delta.dot(proj_ref_delta));
    if (cos_rot_sq < max_cos_rot_sq) {
        // Return failure if the rotation isn't within tolerance.
        return RotationTestResult();
    }
    return RotationTestResult(cos_rot_sq, proj_ref_delta, shift_matrix);
}

Eigen::Matrix3d create_spherical_rotation_matrix(Eigen::Vector3d const& rot_axis, double cos_rotation,
                                                 double sin_rotation) {
    Eigen::Matrix3d rot_cross_matrix;
    // clang-format off
    rot_cross_matrix << 0., -rot_axis[2], rot_axis[1],
                        rot_axis[2], 0., -rot_axis[0],
                        -rot_axis[1], rot_axis[0], 0.;
    // clang-format on
    Eigen::Matrix3d shift_matrix = cos_rotation * Eigen::Matrix3d::Identity() +
                                   sin_rotation * rot_cross_matrix +
                                   (1 - cos_rotation) * rot_axis * rot_axis.transpose();
    return shift_matrix;
}

ShiftRotMatrixResult create_shift_rot_matrix(double cos_rot_sq, Eigen::Matrix3d const& shift_matrix,
                                             ndarray::Array<double, 1, 1> const& src_delta,
                                             ndarray::Array<double, 1, 1> const& ref_ctr,
                                             ndarray::Array<double, 1, 1> const& ref_delta) {
    double cos_rot = sqrt(cos_rot_sq);
    Eigen::Vector3d src_delta_eigen = ndarray::asEigenMatrix(src_delta);
    Eigen::Vector3d rot_src_delta = shift_matrix * src_delta_eigen;
    Eigen::Vector3d tmp_cross = rot_src_delta.cross(ndarray::asEigenMatrix(ref_delta).head<3>());
    double delta_dot_cross = tmp_cross.dot(ndarray::asEigenMatrix(ref_ctr));

    double sin_rot = sgn(delta_dot_cross) * sqrt(1 - cos_rot_sq);
    Eigen::Matrix3d rot_matrix =
            create_spherical_rotation_matrix(ndarray::asEigenMatrix(ref_ctr), cos_rot, sin_rot);

    Eigen::Matrix3d shift_rot_matrix = rot_matrix * shift_matrix;

    ShiftRotMatrixResult result;
    result.sin_rot = sin_rot;
    result.shift_rot_matrix = shift_rot_matrix;
    return result;
}

bool intermediate_verify_comparison(std::vector<Eigen::Vector3d> const& src_pattern,
                                    std::vector<Eigen::Vector3d> const& ref_pattern,
                                    Eigen::Matrix3d const& shift_rot_matrix, double max_dist_rad) {
    double max_dist_sq = max_dist_rad * max_dist_rad;
    bool passed = true;
    auto iSrc = src_pattern.begin();
    auto iRef = ref_pattern.begin();
    for (auto iSrc = src_pattern.begin(), iRef = ref_pattern.begin();
         iSrc != src_pattern.end() && iRef != ref_pattern.end(); iSrc++, iRef++) {
        Eigen::Vector3d rot_src_vect = shift_rot_matrix * *iSrc;
        Eigen::Vector3d diff_vect = rot_src_vect - *iRef;
        if (max_dist_sq < diff_vect.dot(diff_vect)) {
            passed = false;
            break;
        }
    }
    return passed;
}

std::vector<std::pair<size_t, size_t>> create_pattern_spokes(
        ndarray::Array<double, 1, 1> const& src_ctr, ndarray::Array<double, 2, 1> const& src_delta_array,
        ndarray::Array<double, 1, 1> const& src_dist_array, ndarray::Array<double, 1, 1> const& ref_ctr,
        Eigen::Vector3d const& proj_ref_ctr_delta, std::vector<float> const& ref_dist_array,
        std::vector<uint16_t> const& ref_id_array, ndarray::Array<double, 2, 1> const& reference_array,
        double max_dist_rad, size_t n_match) {
    // Struct where we will be putting our results.
    std::vector<std::pair<size_t, size_t>> output_spokes;

    // Counter for number of spokes we failed to find a reference
    // candidate for. We break the loop if we haven't found enough.
    size_t n_fail = 0;

    // Plane project the center/first spoke of the source pattern using
    // the center vector of the pattern as normal.
    auto src_ctr_eigen = ndarray::asEigenMatrix(src_ctr);
    double src_delta_ctr_dot = ndarray::asEigenMatrix(src_delta_array[0]).dot(src_ctr_eigen);

    ndarray::Array<double, 1, 1> proj_src_ctr_delta = copy(src_delta_array[0] - src_delta_ctr_dot * src_ctr);
    auto proj_src_ctr_delta_eigen = ndarray::asEigenMatrix(proj_src_ctr_delta);
    double proj_src_ctr_dist_sq = proj_src_ctr_delta_eigen.dot(proj_src_ctr_delta_eigen);

    // Pre - compute the squared length of the projected reference vector.
    double proj_ref_ctr_dist_sq = proj_ref_ctr_delta.dot(proj_ref_ctr_delta);

    // Value of sin where sin(theta) ~= theta to within 0.1%. Used to make
    // sure we are still in small angle approximation.
    double max_sin_tol = 0.0447;
    // Loop over the source pairs.
    for (size_t src_idx = 1; src_idx < src_dist_array.size(); src_idx++) {
        if (n_fail > src_dist_array.size() - (n_match - 1)) {
            break;
        }
        // Find the reference pairs that include our candidate pattern center
        // and sort them in increasing delta. Check this first so we don't
        // compute anything else if no candidates exist.
        std::pair<size_t, size_t> candidate_range =
                find_candidate_reference_pair_range(src_dist_array[src_idx], ref_dist_array, max_dist_rad);
        if (candidate_range.first == candidate_range.second) {
            n_fail++;
            continue;
        }

        // Given our length tolerance we can use it to compute a tolerance
        // on the angle between our spoke.
        double src_sin_tol = max_dist_rad / (src_dist_array[src_idx] + max_dist_rad);

        // Test if the small angle approximation will still hold. This is
        // defined as when sin(theta) ~= theta to within 0.1% of each
        // other. If the implied opening angle is too large we set it to
        // the 0.1% threshold.
        if (src_sin_tol > max_sin_tol) {
            src_sin_tol = max_sin_tol;
        }

        // Plane project the candidate source spoke and compute the cosine
        // and sine of the opening angle.
        double proj_src_delta_dot = ndarray::asEigenMatrix(src_delta_array[src_idx]).dot(src_ctr_eigen);

        ndarray::Array<double, 1, 1> proj_src_delta =
                copy(src_delta_array[src_idx] - proj_src_delta_dot * src_ctr);
        auto proj_src_delta_eigen = ndarray::asEigenMatrix(proj_src_delta);
        double geom_dist_src = sqrt(proj_src_delta_eigen.dot(proj_src_delta_eigen) * proj_src_ctr_dist_sq);

        // Compute cosine and sine of the delta vector opening angle.
        double cos_theta_src = proj_src_delta_eigen.dot(proj_src_ctr_delta_eigen) / geom_dist_src;
        Eigen::Vector3d cross_src =
                proj_src_delta_eigen.head<3>().cross(proj_src_ctr_delta_eigen.head<3>()) / geom_dist_src;
        double sin_theta_src = cross_src.dot(src_ctr_eigen);

        // Test the spokes and return the id of the reference object.
        // Return -1 if no match is found.
        int ref_id =
                check_spoke(cos_theta_src, sin_theta_src, ref_ctr, proj_ref_ctr_delta, proj_ref_ctr_dist_sq,
                            candidate_range, ref_id_array, reference_array, src_sin_tol);
        if (ref_id < 0) {
            n_fail++;
            continue;
        }

        /// Append the successful indices to our list. The src_idx needs
        /// an extra iteration to skip the first and second source objects.
        output_spokes.push_back(std::make_pair(static_cast<size_t>(ref_id), src_idx + 1));
        if (output_spokes.size() >= n_match - 2) {
            break;
        }
    }
    return output_spokes;
}

int check_spoke(double cos_theta_src, double sin_theta_src, ndarray::Array<double, 1, 1> const& ref_ctr,
                Eigen::Vector3d const& proj_ref_ctr_delta, double proj_ref_ctr_dist_sq,
                std::pair<size_t, size_t> const& candidate_range, std::vector<uint16_t> const& ref_id_array,
                ndarray::Array<double, 2, 1> const& reference_array, double src_sin_tol) {
    // Loop over our candidate reference objects. candidate_range is the min
    // and max of for pair candidates and are view into ref_id_array. Here we
    // start from the midpoint of min and max values and step outward.
    size_t midpoint = (candidate_range.first + candidate_range.second) / 2;
    for (size_t idx = 0; idx < candidate_range.second - candidate_range.first; idx++) {
        // TODO DM-33514: cleanup this loop to use an iterator that handles the "inside-out" iteration.
        if (idx % 2 == 0) {
            midpoint = midpoint + idx;
        } else {
            midpoint = midpoint - idx;
        }
        // Compute the delta vector from the pattern center.
        uint16_t ref_id = ref_id_array[midpoint];
        ndarray::Array<double, 1, 1> ref_delta = copy(reference_array[ref_id] - ref_ctr);

        double ref_dot = ndarray::asEigenMatrix(ref_delta).dot(ndarray::asEigenMatrix(ref_ctr));
        ndarray::Array<double, 1, 1> proj_ref_delta = copy(ref_delta - ref_dot * ref_ctr);
        // Compute the cos between our "center" reference vector and the
        // current reference candidate.
        auto proj_ref_delta_eigen = ndarray::asEigenMatrix(proj_ref_delta);
        double proj_delta_dist_sq = proj_ref_delta_eigen.dot(proj_ref_delta_eigen);
        double geom_dist_ref = sqrt(proj_ref_ctr_dist_sq * proj_delta_dist_sq);
        double cos_theta_ref = proj_ref_delta_eigen.dot(proj_ref_ctr_delta) / geom_dist_ref;

        // Make sure we can safely make the comparison in case
        // our "center" and candidate vectors are mostly aligned.
        double cos_sq_comparison;
        if (cos_theta_ref * cos_theta_ref < (1 - src_sin_tol * src_sin_tol)) {
            cos_sq_comparison = (cos_theta_src - cos_theta_ref) * (cos_theta_src - cos_theta_ref) /
                                (1 - cos_theta_ref * cos_theta_ref);
        } else {
            cos_sq_comparison = (cos_theta_src - cos_theta_ref) * (cos_theta_src - cos_theta_ref) /
                                (src_sin_tol * src_sin_tol);
        }
        // Test the difference of the cosine of the reference angle against
        // the source angle. Assumes that the delta between the two is
        // small.
        if (cos_sq_comparison > src_sin_tol * src_sin_tol) {
            continue;
        }
        // The cosine tests the magnitude of the angle but not
        // its direction. To do that we need to know the sine as well.
        // This cross product calculation does that.
        Eigen::Vector3d cross_ref = proj_ref_delta_eigen.head<3>().cross(proj_ref_ctr_delta) / geom_dist_ref;
        double sin_theta_ref = cross_ref.dot(ndarray::asEigenMatrix(ref_ctr));
        // Check the value of the cos again to make sure that it is not
        // near zero.
        double sin_comparison;
        if (abs(cos_theta_src) < src_sin_tol) {
            sin_comparison = (sin_theta_src - sin_theta_ref) / src_sin_tol;
        } else {
            sin_comparison = (sin_theta_src - sin_theta_ref) / cos_theta_ref;
        }

        // Return the correct id of the candidate we found.
        if (abs(sin_comparison) < src_sin_tol) {
            return ref_id;
        }
    }
    return -1;
}

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
