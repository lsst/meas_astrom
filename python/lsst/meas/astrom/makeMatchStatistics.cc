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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

#include "lsst/meas/astrom/makeMatchStatistics.h"

namespace lsst {
namespace meas {
namespace astrom {

namespace {

template<typename MatchT>
void declareMakeMatchStatistics(py::module & mod) {
    mod.def("makeMatchStatistics", &makeMatchStatistics<MatchT>,
            "matchList"_a, "flags"_a, "sctrl"_a=afw::math::StatisticsControl());
    mod.def("makeMatchStatisticsInPixels", &makeMatchStatisticsInPixels<MatchT>,
            "wcs"_a, "matchList"_a, "flags"_a, "sctrl"_a=afw::math::StatisticsControl());
    mod.def("makeMatchStatisticsInRadians", &makeMatchStatisticsInRadians<MatchT>,
            "wcs"_a, "matchList"_a, "flags"_a, "sctrl"_a=afw::math::StatisticsControl());
}

}  // namespace lsst::meas::astrom::<anonymous>

PYBIND11_PLUGIN(_makeMatchStatistics) {
    py::module mod("_makeMatchStatistics", "Python wrapper for afw _makeMatchStatistics library");

    declareMakeMatchStatistics<afw::table::ReferenceMatch>(mod);
    declareMakeMatchStatistics<afw::table::SourceMatch>(mod);

    return mod.ptr();
}

}}}  // namespace lsst::meas::astrom
