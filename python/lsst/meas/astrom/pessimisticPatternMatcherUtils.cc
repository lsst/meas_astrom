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

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/meas/astrom/pessimisticPatternMatcherUtils.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {

PYBIND11_MODULE(pessimisticPatternMatcherUtils, mod) {
    mod.def("find_candidate_reference_pair_range", &find_candidate_reference_pair_range, "src_dist"_a,
            "ref_dist_array"_a, "max_dist_rad"_a);
    mod.def("create_pattern_spokes", &create_pattern_spokes, "src_ctr"_a, "src_delta_array"_a,
            "src_dist_array"_a, "ref_ctr"_a, "proj_ref_ctr_delta"_a, "ref_dist_array"_a, "ref_id_array"_a,
            "reference_array"_a, "max_dist_rad"_a, "n_match"_a);
    mod.def("check_spoke", &check_spoke, "cos_theta_src"_a, "sin_theta_src"_a, "ref_ctr"_a,
            "proj_ref_ctr_delta"_a, "proj_ref_ctr_dist_sq"_a, "candidate_range"_a, "ref_id_array"_a,
            "reference_array"_a, "src_sin_tol"_a);
}

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
