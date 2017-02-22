/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lsst/pex/config/python.h" // defines LSST_DECLARE_CONTROL_FIELD
#include "lsst/afw/table/Match.h"
#include "lsst/meas/astrom/ScaledPolynomialTransformFitter.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {

namespace {

void declareOutlierRejectionControl(py::module & mod) {
    py::class_<OutlierRejectionControl> cls(mod, "OutlierRejectionControl");

    cls.def(py::init<>());

    LSST_DECLARE_CONTROL_FIELD(cls, OutlierRejectionControl, nSigma);
    LSST_DECLARE_CONTROL_FIELD(cls, OutlierRejectionControl, nClipMin);
    LSST_DECLARE_CONTROL_FIELD(cls, OutlierRejectionControl, nClipMax);
}

void declareScaledPolynomialTransformFitter(py::module & mod) {
    py::class_<ScaledPolynomialTransformFitter> cls(mod, "ScaledPolynomialTransformFitter");

    cls.def_static("fromMatches", &ScaledPolynomialTransformFitter::fromMatches);
    cls.def_static("fromGrid", &ScaledPolynomialTransformFitter::fromGrid);

    cls.def("fit", &ScaledPolynomialTransformFitter::fit, "order"_a=-1);
    cls.def("updateModel", &ScaledPolynomialTransformFitter::updateModel);
    cls.def("updateIntrinsicScatter", &ScaledPolynomialTransformFitter::updateIntrinsicScatter);
    cls.def("getIntrinsicScatter", &ScaledPolynomialTransformFitter::getIntrinsicScatter);
    cls.def("rejectOutliers", &ScaledPolynomialTransformFitter::rejectOutliers, "ctrl"_a);
    cls.def("getData", &ScaledPolynomialTransformFitter::getData,
            py::return_value_policy::reference_internal);
    cls.def("getTransform", &ScaledPolynomialTransformFitter::getTransform, py::return_value_policy::copy);
    cls.def("getPoly", &ScaledPolynomialTransformFitter::getPoly, py::return_value_policy::copy);
    cls.def("getInputScaling", &ScaledPolynomialTransformFitter::getInputScaling,
            py::return_value_policy::copy);
    cls.def("getOutputScaling", &ScaledPolynomialTransformFitter::getOutputScaling,
            py::return_value_policy::copy);
}

}  // namespace lsst::meas::astrom::<anonymous>

PYBIND11_PLUGIN(scaledPolynomialTransformFitter) {
    py::module mod("scaledPolynomialTransformFitter");

    declareOutlierRejectionControl(mod);
    declareScaledPolynomialTransformFitter(mod);

    return mod.ptr();
}

}}}  // namespace lsst::meas::astrom
