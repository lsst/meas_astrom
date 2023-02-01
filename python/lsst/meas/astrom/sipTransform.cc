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

#include "lsst/geom/Point.h"
#include "lsst/geom/LinearTransform.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/meas/astrom/SipTransform.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {
namespace {

void declareSipTransformBase(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PySipTransformBase = py::class_<SipTransformBase, std::shared_ptr<SipTransformBase>>;

    wrappers.wrapType(PySipTransformBase(wrappers.module, "_SipTransformBase"), [](auto &mod, auto &cls) {
        cls.def("getPixelOrigin", &SipTransformBase::getPixelOrigin, py::return_value_policy::copy);
        cls.def("getCdMatrix", &SipTransformBase::getCdMatrix, py::return_value_policy::copy);
        cls.def("getPoly", &SipTransformBase::getPoly, py::return_value_policy::copy);
    });
}

void declareSipForwardTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PySipForwardTransform =
            py::class_<SipForwardTransform, std::shared_ptr<SipForwardTransform>, SipTransformBase>;

    wrappers.wrapType(PySipForwardTransform(wrappers.module, "SipForwardTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<geom::Point2D const &, geom::LinearTransform const &, PolynomialTransform const &>(),
                "pixelOrigin"_a, "cdMatrix"_a, "forwardSipPoly"_a);
        cls.def(py::init<SipForwardTransform const &>(), "other"_a);

        cls.def_static("convert",
                       (SipForwardTransform(*)(PolynomialTransform const &, geom::Point2D const &,
                                               geom::LinearTransform const &)) &
                               SipForwardTransform::convert,
                       "poly"_a, "pixelOrigin"_a, "cdMatrix"_a);
        cls.def_static("convert",
                       (SipForwardTransform(*)(ScaledPolynomialTransform const &, geom::Point2D const &,
                                               geom::LinearTransform const &)) &
                               SipForwardTransform::convert,
                       "scaled"_a, "pixelOrigin"_a, "cdMatrix"_a);
        cls.def_static("convert",
                       (SipForwardTransform(*)(ScaledPolynomialTransform const &)) &SipForwardTransform::convert,
                       "scaled"_a);

        cls.def("__call__", &SipForwardTransform::operator(), "in"_a);
        cls.def("transformPixels", &SipForwardTransform::transformPixels, "s"_a);

        cls.def("linearize", &SipForwardTransform::linearize);
    });
}

void declareSipReverseTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PySipReverseTransform =
            py::class_<SipReverseTransform, std::shared_ptr<SipReverseTransform>, SipTransformBase>;

    wrappers.wrapType(PySipReverseTransform(wrappers.module, "SipReverseTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<geom::Point2D const &, geom::LinearTransform const &, PolynomialTransform const &>(),
                "pixelOrigin"_a, "cdMatrix"_a, "reverseSipPoly"_a);
        cls.def(py::init<SipReverseTransform const &>(), "other"_a);

        cls.def_static("convert",
                       (SipReverseTransform(*)(PolynomialTransform const &, geom::Point2D const &,
                                               geom::LinearTransform const &)) &
                               SipReverseTransform::convert,
                       "poly"_a, "pixelOrigin"_a, "cdMatrix"_a);
        cls.def_static("convert",
                       (SipReverseTransform(*)(ScaledPolynomialTransform const &, geom::Point2D const &,
                                               geom::LinearTransform const &)) &
                               SipReverseTransform::convert,
                       "scaled"_a, "pixelOrigin"_a, "cdMatrix"_a);
        cls.def_static("convert",
                       (SipReverseTransform(*)(ScaledPolynomialTransform const &)) &SipReverseTransform::convert,
                       "scaled"_a);

        cls.def("__call__", &SipReverseTransform::operator(), "in"_a);
        cls.def("transformPixels", &SipReverseTransform::transformPixels, "s"_a);

        cls.def("linearize", &SipReverseTransform::linearize);
    });
}

}  // namespace

void wrapSipTransform(lsst::cpputils::python::WrapperCollection &wrappers){
    declareSipTransformBase(wrappers);
    declareSipForwardTransform(wrappers);
    declareSipReverseTransform(wrappers);

    wrappers.module.def("makeWcs", makeWcs, "sipForward"_a, "sipReverse"_a, "skyOrigin"_a);
    wrappers.module.def("transformWcsPixels", transformWcsPixels, "wcs"_a, "s"_a);
    wrappers.module.def("rotateWcsPixelsBy90", rotateWcsPixelsBy90, "wcs"_a, "nQuarter"_a, "dimensions"_a);
}

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
