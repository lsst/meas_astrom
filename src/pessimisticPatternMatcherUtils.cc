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

#include <cmath>
#include <Eigen/Dense>
#include "ndarray/eigen.h"
#include "lsst/meas/astrom/pessimisticPatternMatcherUtils.h"

namespace lsst {
namespace meas {
namespace astrom {

std::pair<size_t, size_t> find_candidate_reference_pair_range(
        float src_dist, ndarray::Array<float, 1, 1> const& ref_dist_array, float max_dist_rad) {
    auto itr =
            std::upper_bound(ref_dist_array.begin(), ref_dist_array.end(), src_dist - max_dist_rad - 1e-16);
    auto itrEnd =
            std::upper_bound(ref_dist_array.begin(), ref_dist_array.end(), src_dist + max_dist_rad + 1e-16);

    size_t startIdx = itr - ref_dist_array.begin();
    size_t endIdx = itrEnd - ref_dist_array.begin();
    return std::make_pair(startIdx, endIdx);
}

std::vector<std::pair<size_t, size_t>> create_pattern_spokes(
        ndarray::Array<double, 1, 1> const& src_ctr, ndarray::Array<double, 2, 1> const& src_delta_array,
        ndarray::Array<double, 1, 1> const& src_dist_array, ndarray::Array<double, 1, 1> const& ref_ctr,
        ndarray::Array<double, 1, 1> const& proj_ref_ctr_delta,
        ndarray::Array<float, 1, 1> const& ref_dist_array, ndarray::Array<uint16_t, 1, 1> const& ref_id_array,
        ndarray::Array<double, 2, 1> const& reference_array, double max_dist_rad, size_t n_match) {
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
    auto proj_ref_ctr_delta_eigen = ndarray::asEigenMatrix(proj_ref_ctr_delta);
    double proj_ref_ctr_dist_sq = proj_ref_ctr_delta_eigen.dot(proj_ref_ctr_delta_eigen);

    // Loop over the source pairs.
    // Value of sin where sin(theta) ~= theta to within 0.1%.
    double max_sin_tol = 0.0447;
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
                ndarray::Array<double, 1, 1> const& proj_ref_ctr_delta, double proj_ref_ctr_dist_sq,
                std::pair<size_t, size_t> const& candidate_range,
                ndarray::Array<uint16_t, 1, 1> const& ref_id_array,
                ndarray::Array<double, 2, 1> const& reference_array, double src_sin_tol) {
    // Loop over our candidate reference objects. candidate_range is the min
    // and max of for pair candidates and are view into ref_id_array. Here we
    // start from the midpoint of min and max values and step outward.
    size_t midpoint = (candidate_range.first + candidate_range.second) / 2;
    for (size_t idx = 0; idx < candidate_range.second - candidate_range.first; idx++) {
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
        auto proj_ref_ctr_delta_eigen = ndarray::asEigenMatrix(proj_ref_ctr_delta);
        double proj_delta_dist_sq = proj_ref_delta_eigen.dot(proj_ref_delta_eigen);
        double geom_dist_ref = sqrt(proj_ref_ctr_dist_sq * proj_delta_dist_sq);
        double cos_theta_ref = proj_ref_delta_eigen.dot(proj_ref_ctr_delta_eigen) / geom_dist_ref;

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
        Eigen::Vector3d cross_ref =
                proj_ref_delta_eigen.head<3>().cross(proj_ref_ctr_delta_eigen.head<3>()) / geom_dist_ref;
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
