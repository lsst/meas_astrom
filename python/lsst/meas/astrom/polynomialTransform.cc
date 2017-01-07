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
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"
#include "ndarray/converter.h"

#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/meas/astrom/SipTransform.h"
#include "lsst/meas/astrom/PolynomialTransform.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {

namespace {

void declarePolynomialTransform(py::module & mod) {
    py::class_<PolynomialTransform, std::shared_ptr<PolynomialTransform>> cls(mod, "PolynomialTransform");

    cls.def(py::init<ndarray::Array<double const,2,2> const &, ndarray::Array<double const,2,2> const &>(),
            "xCoeffs"_a, "yCoeffs"_a);
    cls.def(py::init<PolynomialTransform const &>(), "other"_a);

    cls.def_static("convert",
                   (PolynomialTransform (*)(ScaledPolynomialTransform const &)) &PolynomialTransform::convert,
                   "other"_a);
    cls.def_static("convert",
                   (PolynomialTransform (*)(SipForwardTransform const &)) &PolynomialTransform::convert,
                   "other"_a);
    cls.def_static("convert",
                   (PolynomialTransform (*)(SipReverseTransform const &)) &PolynomialTransform::convert,
                   "other"_a);

    cls.def("__call__", &PolynomialTransform::operator(), "in"_a, py::is_operator());

    cls.def("getOrder", &PolynomialTransform::getOrder);
    cls.def("getXCoeffs", &PolynomialTransform::getXCoeffs);
    cls.def("getYCoeffs", &PolynomialTransform::getYCoeffs);
    cls.def("linearize", &PolynomialTransform::linearize);
}

void declareScaledPolynomialTransform(py::module & mod) {
    py::class_<ScaledPolynomialTransform,
               std::shared_ptr<ScaledPolynomialTransform>> cls(mod, "ScaledPolynomialTransform");

    cls.def(py::init<PolynomialTransform const &,
                     afw::geom::AffineTransform const &,
                     afw::geom::AffineTransform const &>(),
            "poly"_a, "inputScaling"_a, "outputScalingInverse"_a);

    cls.def(py::init<ScaledPolynomialTransform const &>(), "other"_a);

    cls.def_static("convert",
                   (ScaledPolynomialTransform (*)(PolynomialTransform const &))
                        &ScaledPolynomialTransform::convert,
                   "other"_a);
    cls.def_static("convert",
                   (ScaledPolynomialTransform (*)(SipForwardTransform const &))
                        &ScaledPolynomialTransform::convert,
                   "other"_a);
    cls.def_static("convert",
                   (ScaledPolynomialTransform (*)(SipReverseTransform const &))
                        &ScaledPolynomialTransform::convert,
                   "other"_a);

    cls.def("__call__", &ScaledPolynomialTransform::operator(), "in"_a, py::is_operator());

    cls.def("getPoly", &ScaledPolynomialTransform::getPoly, py::return_value_policy::reference_internal);
    cls.def("getInputScaling", &ScaledPolynomialTransform::getInputScaling,
            py::return_value_policy::reference_internal);
    cls.def("getOutputScalingInverse", &ScaledPolynomialTransform::getOutputScalingInverse,
            py::return_value_policy::reference_internal);

    cls.def("linearize", &ScaledPolynomialTransform::linearize);
}

}  // namespace lsst::meas::astrom::<anonymous>

PYBIND11_PLUGIN(_polynomialTransform) {
    py::module mod("_polynomialTransform", "Python wrapper for afw _polynomialTransform library");

    if (_import_array() < 0) {
            PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
            return nullptr;
    };

    declarePolynomialTransform(mod);
    declareScaledPolynomialTransform(mod);

    mod.def("compose", (PolynomialTransform (*)(afw::geom::AffineTransform const &,
                                                PolynomialTransform const &)) &compose,
            "t1"_a, "t2"_a);
    mod.def("compose", (PolynomialTransform (*)(PolynomialTransform const &,
                                                afw::geom::AffineTransform const &)) &compose,
            "t1"_a, "t2"_a);

    return mod.ptr();
}

}}}  // namespace lsst::meas::astrom
