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
#include "pybind11/eigen.h"
#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/meas/astrom/pessimisticPatternMatcherUtils.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {

PYBIND11_MODULE(pessimisticPatternMatcherUtils, mod) {
    py::class_<PatternResult>(mod, "PatternResult")
            .def_readonly("candidate_pairs", &PatternResult::candidate_pairs)
            .def_readonly("shift_rot_matrix", &PatternResult::shift_rot_matrix)
            .def_readonly("cos_shift", &PatternResult::cos_shift)
            .def_readonly("sin_rot", &PatternResult::sin_rot)
            .def_readonly("success", &PatternResult::success);
    mod.def("construct_pattern_and_shift_rot_matrix", &construct_pattern_and_shift_rot_matrix,
            "src_pattern_array"_a, "src_delta_array"_a, "src_dist_array"_a, "dist_array"_a, "id_array"_a,
            "reference_array"_a, "n_match"_a, "max_cos_theta_shift"_a, "max_cos_rot_sq"_a, "max_dist_rad"_a);
}

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
