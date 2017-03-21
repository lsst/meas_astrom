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
#include "pybind11/stl.h"

#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/astrom/SipTransform.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {
namespace {

static void declareSipTransformBase(py::module &mod) {
    py::class_<SipTransformBase, std::shared_ptr<SipTransformBase>> cls(mod, "_SipTransformBase");

    cls.def("getPixelOrigin", &SipTransformBase::getPixelOrigin, py::return_value_policy::copy);
    cls.def("getCDMatrix", &SipTransformBase::getCDMatrix, py::return_value_policy::copy);
    cls.def("getPoly", &SipTransformBase::getPoly, py::return_value_policy::copy);
}

static void declareSipForwardTransform(py::module &mod) {
    py::class_<SipForwardTransform, std::shared_ptr<SipForwardTransform>, SipTransformBase> cls(
            mod, "SipForwardTransform");

    cls.def(py::init<afw::geom::Point2D const &, afw::geom::LinearTransform const &,
                     PolynomialTransform const &>(),
            "pixelOrigin"_a, "cdMatrix"_a, "forwardSipPoly"_a);
    cls.def(py::init<SipForwardTransform const &>(), "other"_a);

    cls.def_static("convert", (SipForwardTransform(*)(PolynomialTransform const &, afw::geom::Point2D const &,
                                                      afw::geom::LinearTransform const &)) &
                                      SipForwardTransform::convert,
                   "poly"_a, "pixelOrigin"_a, "cdMatrix"_a);
    cls.def_static("convert",
                   (SipForwardTransform(*)(ScaledPolynomialTransform const &, afw::geom::Point2D const &,
                                           afw::geom::LinearTransform const &)) &
                           SipForwardTransform::convert,
                   "scaled"_a, "pixelOrigin"_a, "cdMatrix"_a);
    cls.def_static("convert",
                   (SipForwardTransform(*)(ScaledPolynomialTransform const &)) & SipForwardTransform::convert,
                   "scaled"_a);
    cls.def_static("extract", &SipForwardTransform::extract, "wcs"_a);

    cls.def("__call__", &SipForwardTransform::operator(), "in"_a);
    cls.def("transformPixels", &SipForwardTransform::transformPixels, "s"_a);

    cls.def("linearize", &SipForwardTransform::linearize);
}

static void declareSipReverseTransform(py::module &mod) {
    py::class_<SipReverseTransform, std::shared_ptr<SipReverseTransform>, SipTransformBase> cls(
            mod, "SipReverseTransform");

    cls.def(py::init<afw::geom::Point2D const &, afw::geom::LinearTransform const &,
                     PolynomialTransform const &>(),
            "pixelOrigin"_a, "cdMatrix"_a, "reverseSipPoly"_a);
    cls.def(py::init<SipReverseTransform const &>(), "other"_a);

    cls.def_static("convert", (SipReverseTransform(*)(PolynomialTransform const &, afw::geom::Point2D const &,
                                                      afw::geom::LinearTransform const &)) &
                                      SipReverseTransform::convert,
                   "poly"_a, "pixelOrigin"_a, "cdMatrix"_a);
    cls.def_static("convert",
                   (SipReverseTransform(*)(ScaledPolynomialTransform const &, afw::geom::Point2D const &,
                                           afw::geom::LinearTransform const &)) &
                           SipReverseTransform::convert,
                   "scaled"_a, "pixelOrigin"_a, "cdMatrix"_a);
    cls.def_static("convert",
                   (SipReverseTransform(*)(ScaledPolynomialTransform const &)) & SipReverseTransform::convert,
                   "scaled"_a);
    cls.def_static("extract", &SipReverseTransform::extract, "wcs"_a);

    cls.def("__call__", &SipReverseTransform::operator(), "in"_a);
    cls.def("transformPixels", &SipReverseTransform::transformPixels, "s"_a);

    cls.def("linearize", &SipReverseTransform::linearize);
}

}  // namespace lsst::meas::astrom::<anonymous>

PYBIND11_PLUGIN(sipTransform) {
    py::module mod("sipTransform");

    declareSipTransformBase(mod);
    declareSipForwardTransform(mod);
    declareSipReverseTransform(mod);

    mod.def("makeWcs", &makeWcs);
    mod.def("transformWcsPixels", &transformWcsPixels, "wcs"_a, "s"_a);

    return mod.ptr();
}
}
}
}  // namespace lsst::meas::astrom
