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
#include <vector>

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/table/Match.h"
#include "lsst/meas/astrom/sip/CreateWcsWithSip.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {
namespace sip {
namespace {

template <typename MatchT>
static void declareCreateWcsWithSip(py::module &mod, std::string const &name) {
    py::class_<CreateWcsWithSip<MatchT>, std::shared_ptr<CreateWcsWithSip<MatchT>>> cls(mod, name.c_str());

    cls.def(py::init<std::vector<MatchT> const &, afw::image::Wcs const &, int const,
                     afw::geom::Box2I const &, int const>(),
            "matches"_a, "linearWcs"_a, "order"_a, "bbox"_a = afw::geom::Box2I(), "ngrid"_a = 0);

    cls.def("getNewWcs", &CreateWcsWithSip<MatchT>::getNewWcs);
    cls.def("getScatterInPixels", &CreateWcsWithSip<MatchT>::getScatterInPixels);
    cls.def("getScatterOnSky", &CreateWcsWithSip<MatchT>::getScatterOnSky);
    cls.def("getLinearScatterInPixels", &CreateWcsWithSip<MatchT>::getLinearScatterInPixels);
    cls.def("getLinearScatterOnSky", &CreateWcsWithSip<MatchT>::getLinearScatterOnSky);
    cls.def("getOrder", &CreateWcsWithSip<MatchT>::getOrder);
    cls.def("getNPoints", &CreateWcsWithSip<MatchT>::getNPoints);
    cls.def("getNGrid", &CreateWcsWithSip<MatchT>::getNGrid);
    cls.def("getSipA", &CreateWcsWithSip<MatchT>::getSipA, py::return_value_policy::copy);
    cls.def("getSipB", &CreateWcsWithSip<MatchT>::getSipB, py::return_value_policy::copy);
    cls.def("getSipAp", &CreateWcsWithSip<MatchT>::getSipAp, py::return_value_policy::copy);
    cls.def("getSipBp", &CreateWcsWithSip<MatchT>::getSipBp, py::return_value_policy::copy);

    mod.def("makeCreateWcsWithSip", &makeCreateWcsWithSip<MatchT>, "matches"_a, "linearWcs"_a, "order"_a,
            "bbox"_a = afw::geom::Box2I(), "ngrid"_a = 0);
}

}  // namespace lsst::meas::astrom::sip::<anonymous>

PYBIND11_PLUGIN(createWcsWithSip) {
    py::module mod("createWcsWithSip");

    // Need to import numpy for ndarray and eigen conversions
    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    declareCreateWcsWithSip<afw::table::ReferenceMatch>(mod, "CreateWcsWithSipReferenceMatch");
    declareCreateWcsWithSip<afw::table::SourceMatch>(mod, "CreateWcsWithSipSourceMatch");

    return mod.ptr();
}
}
}
}
}  // namespace lsst::meas::astrom::sip
