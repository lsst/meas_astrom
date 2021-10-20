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
 * @param[in] ref_ctr_id id lookup of the ref_ctr into the master reference
 *     array. uint16 type is to limit the amount of memory used and is set
 *     in the pessimistic_pattern_matcher_3d python class with reference
 *     catalogs trimmed by the matchPessimisticB runner class.
 * @param[in] proj_ref_ctr_delta Plane projected first spoke in the reference
 *     pattern using the pattern center as normal.
 * @param[in] proj_ref_ctr_dist_sq Squared length of the projected vector.
 * @param[in] ref_dist_idx_array Indices sorted by the delta distance between
 *     the source spoke we are trying to test and the candidate reference
 *     spokes.
 * @param[in] ref_id_array Array of id lookups into the master reference array
 *     that our center id object is paired with.
 * @param[in] reference_array Array of three vectors representing the locations
 *     of all reference objects.
 * @param[in] src_sin_tol Sine of tolerance allowed between source and
 *     reference spoke opening angles.
 * @return ID of the candidate reference object successfully matched or -1 if
 *     no match is found.
 */
int check_spoke(double cos_theta_src, double sin_theta_src, ndarray::Array<double, 1, 1> const& ref_ctr,
                ndarray::Array<double, 1, 1> const& proj_ref_ctr_delta, double proj_ref_ctr_dist_sq,
                ndarray::Array<long int, 1, 1> const& ref_dist_idx_array,
                ndarray::Array<uint16_t, 1, 1> const& ref_id_array,
                ndarray::Array<double, 2, 1> const& reference_array, double src_sin_tol);

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
