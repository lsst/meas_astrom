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

#include <string>

#include "ndarray/pybind11.h"

#include "lsst/meas/astrom/sip/LeastSqFitter1d.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {
namespace sip {
namespace {

template <typename FittingFunc>
static void declareLeastSqFitter1d(py::module &mod, std::string const &name) {
    py::class_<LeastSqFitter1d<FittingFunc>, std::shared_ptr<LeastSqFitter1d<FittingFunc>>> cls(mod,
                                                                                                name.c_str());

    cls.def(py::init<std::vector<double> const &, std::vector<double> const &, std::vector<double> const &,
                     int>(),
            "x"_a, "y"_a, "s"_a, "order"_a);

    cls.def("getParams", &LeastSqFitter1d<FittingFunc>::getParams);
    cls.def("getErrors", &LeastSqFitter1d<FittingFunc>::getErrors);
    cls.def("getBestFitFunction", &LeastSqFitter1d<FittingFunc>::getBestFitFunction);
    cls.def("valueAt", &LeastSqFitter1d<FittingFunc>::valueAt, "x"_a);
    cls.def("residuals", &LeastSqFitter1d<FittingFunc>::residuals);
    cls.def("getChiSq", &LeastSqFitter1d<FittingFunc>::getChiSq);
    cls.def("getReducedChiSq", &LeastSqFitter1d<FittingFunc>::getReducedChiSq);
}

}  // namespace

PYBIND11_PLUGIN(leastSqFitter1d) {
    py::module mod("leastSqFitter1d");

    declareLeastSqFitter1d<afw::math::PolynomialFunction1<double>>(mod, "LeastSqFitter1dPoly");

    return mod.ptr();
}
}  // namespace sip
}  // namespace astrom
}  // namespace meas
}  // namespace lsst
