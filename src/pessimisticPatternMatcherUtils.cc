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

int check_spoke(double cos_theta_src, double sin_theta_src, ndarray::Array<double, 1, 1> const& ref_ctr,
                ndarray::Array<double, 1, 1> const& proj_ref_ctr_delta, double proj_ref_ctr_dist_sq,
                ndarray::Array<long int, 1, 1> const& ref_dist_idx_array,
                ndarray::Array<uint16_t, 1, 1> const& ref_id_array,
                ndarray::Array<double, 2, 1> const& reference_array, double src_sin_tol) {
    // Loop over our candidate reference objects. ref_dist_idx_array is a view into ref_id_array
    // and is not the same length as ref_id_array.
    for (auto idx = ref_dist_idx_array.begin(); idx != ref_dist_idx_array.end(); idx++) {
        // Compute the delta vector from the pattern center.
        unsigned int ref_id = ref_id_array[*idx];
        ndarray::Array<double, 1, 1> ref_delta = copy(reference_array[ref_id] - ref_ctr);

        double ref_dot = ndarray::asEigenMatrix(ref_delta).dot(ndarray::asEigenMatrix(ref_ctr));
        ndarray::Array<double, 1, 1> proj_ref_delta = copy(ref_delta - ref_dot * ref_ctr);
        // Compute the cos between our "center" reference vector and the
        // current reference candidate.
        double proj_delta_dist_sq =
                ndarray::asEigenMatrix(proj_ref_delta).dot(ndarray::asEigenMatrix(proj_ref_delta));
        double geom_dist_ref = sqrt(proj_ref_ctr_dist_sq * proj_delta_dist_sq);
        double cos_theta_ref =
                ndarray::asEigenMatrix(proj_ref_delta).dot(ndarray::asEigenMatrix(proj_ref_ctr_delta)) /
                geom_dist_ref;

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
        Eigen::Vector3d cross_ref = ndarray::asEigenMatrix(proj_ref_delta)
                                            .head<3>()
                                            .cross(ndarray::asEigenMatrix(proj_ref_ctr_delta).head<3>()) /
                                    geom_dist_ref;
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
