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

#include <memory>

#include "lsst/geom/Point.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/pex/config/python.h"  // defines LSST_DECLARE_CONTROL_FIELD
#include "lsst/meas/astrom/matchOptimisticB.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {
namespace {

static void declareRecordProxy(py::module &mod) {
    py::class_<RecordProxy, std::shared_ptr<RecordProxy>> cls(mod, "RecordProxy");

    cls.def_readwrite("record", &RecordProxy::record);
    cls.def_readwrite("position", &RecordProxy::position);
    cls.def_readwrite("used", &RecordProxy::used);

    cls.def(py::init<std::shared_ptr<afw::table::SimpleRecord>, geom::Point2D const &>(), "record"_a,
            "position"_a);

    // TO DO: decide if we need to wrap operator PTR(lsst::afw::table::SimpleRecord)()

    cls.def("__eq__", &RecordProxy::operator==, py::is_operator());
    cls.def("__ne__", &RecordProxy::operator!=, py::is_operator());

    cls.def("getX", &RecordProxy::getX);
    cls.def("getY", &RecordProxy::getY);
}

static void declareProxyPair(py::module &mod) {
    py::class_<ProxyPair, std::shared_ptr<ProxyPair>> cls(mod, "ProxyPair");

    cls.def_readwrite("first", &ProxyPair::first);
    cls.def_readwrite("second", &ProxyPair::second);
    cls.def_readwrite("distance", &ProxyPair::distance);
    cls.def_readwrite("pa", &ProxyPair::pa);

    cls.def(py::init<RecordProxy const &, RecordProxy const &>(), "s1"_a, "s2"_a);
}

static void declareMatchOptimisticBControl(py::module &mod) {
    py::class_<MatchOptimisticBControl> cls(mod, "MatchOptimisticBControl");

    cls.def(py::init<>());

    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, refFluxField);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, sourceFluxField);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, numBrightStars);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, minMatchedPairs);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, matchingAllowancePix);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, maxOffsetPix);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, maxRotationDeg);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, allowedNonperpDeg);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, numPointsForShape);
    LSST_DECLARE_CONTROL_FIELD(cls, MatchOptimisticBControl, maxDeterminant);

    cls.def("validate", &MatchOptimisticBControl::validate);
}

}  // namespace lsst::meas::astrom::<anonymous>

PYBIND11_PLUGIN(matchOptimisticB) {
    py::module mod("matchOptimisticB");

    declareRecordProxy(mod);
    declareProxyPair(mod);
    declareMatchOptimisticBControl(mod);

    mod.def("makeProxies", (ProxyVector(*)(afw::table::SourceCatalog const &, afw::geom::SkyWcs const &,
                                           afw::geom::SkyWcs const &)) &
                                   makeProxies,
            "sourceCat"_a, "distortedWcs"_a, "tanWcs"_a);
    mod.def("makeProxies",
            (ProxyVector(*)(afw::table::SimpleCatalog const &, afw::geom::SkyWcs const &)) & makeProxies,
            "posRefCat"_a, "tanWcs"_a);

    mod.def("matchOptimisticB", &matchOptimisticB, "posRefCat"_a, "sourceCat"_a, "control"_a, "wcs"_a,
            "posRefBegInd"_a = 0, "verbose"_a = false);

    return mod.ptr();
}
}
}
}  // namespace lsst::meas::astrom
