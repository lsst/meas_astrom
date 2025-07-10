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
#include "pybind11/pybind11.h"
#include "lsst/cpputils/python.h"
#include "pybind11/stl.h"

#include <memory>

#include "ndarray/pybind11.h"

#include "lsst/geom/AffineTransform.h"
#include "lsst/meas/astrom/SipTransform.h"
#include "lsst/meas/astrom/PolynomialTransform.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {

namespace {

void declarePolynomialTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPolynomialTransform = py::class_<PolynomialTransform>;

    wrappers.wrapType(PyPolynomialTransform(wrappers.module, "PolynomialTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<ndarray::Array<double const, 2, 0> const &,
                        ndarray::Array<double const, 2, 0> const &>(),
                "xCoeffs"_a, "yCoeffs"_a);
        cls.def(py::init<PolynomialTransform const &>(), "other"_a);

        cls.def_static("convert",
                       (PolynomialTransform(*)(ScaledPolynomialTransform const &)) &PolynomialTransform::convert,
                       "other"_a);
        cls.def_static("convert",
                       (PolynomialTransform(*)(SipForwardTransform const &)) &PolynomialTransform::convert,
                       "other"_a);
        cls.def_static("convert",
                       (PolynomialTransform(*)(SipReverseTransform const &)) &PolynomialTransform::convert,
                       "other"_a);

        cls.def("__call__", &PolynomialTransform::operator(), "in"_a);

        cls.def("getOrder", &PolynomialTransform::getOrder);
        cls.def("getXCoeffs", &PolynomialTransform::getXCoeffs);
        cls.def("getYCoeffs", &PolynomialTransform::getYCoeffs);
        cls.def("linearize", &PolynomialTransform::linearize);
    });
}

void declareScaledPolynomialTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyClass = py::class_<ScaledPolynomialTransform>;

    wrappers.wrapType(PyClass(wrappers.module, "ScaledPolynomialTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<PolynomialTransform const &, geom::AffineTransform const &,
                        geom::AffineTransform const &>(),
                "poly"_a, "inputScaling"_a, "outputScalingInverse"_a);
        cls.def(py::init<ScaledPolynomialTransform const &>(), "other"_a);

        cls.def_static(
                "convert",
                (ScaledPolynomialTransform(*)(PolynomialTransform const &)) &ScaledPolynomialTransform::convert,
                "other"_a);
        cls.def_static(
                "convert",
                (ScaledPolynomialTransform(*)(SipForwardTransform const &)) &ScaledPolynomialTransform::convert,
                "other"_a);
        cls.def_static(
                "convert",
                (ScaledPolynomialTransform(*)(SipReverseTransform const &)) &ScaledPolynomialTransform::convert,
                "other"_a);

        cls.def("__call__", &ScaledPolynomialTransform::operator(), "in"_a);

        cls.def("getPoly", &ScaledPolynomialTransform::getPoly, py::return_value_policy::reference_internal);
        cls.def("getInputScaling", &ScaledPolynomialTransform::getInputScaling,
                py::return_value_policy::reference_internal);
        cls.def("getOutputScalingInverse", &ScaledPolynomialTransform::getOutputScalingInverse,
                py::return_value_policy::reference_internal);
        cls.def("linearize", &ScaledPolynomialTransform::linearize);
    });
}

}  // namespace

void wrapPolynomialTransform(lsst::cpputils::python::WrapperCollection &wrappers){
    declarePolynomialTransform(wrappers);
    declareScaledPolynomialTransform(wrappers);

    wrappers.module.def("compose",
            (PolynomialTransform(*)(geom::AffineTransform const &, PolynomialTransform const &)) & compose,
            "t1"_a, "t2"_a);
    wrappers.module.def("compose",
            (PolynomialTransform(*)(PolynomialTransform const &, geom::AffineTransform const &)) & compose,
            "t1"_a, "t2"_a);
}

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
