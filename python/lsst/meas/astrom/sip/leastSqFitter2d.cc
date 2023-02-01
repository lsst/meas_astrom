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
#include "pybind11/eigen.h"
#include "pybind11/stl.h"

#include <string>

#include "ndarray/pybind11.h"

#include "lsst/meas/astrom/sip/LeastSqFitter2d.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {
namespace sip {
namespace {

template <typename FittingFunc>
void declareLeastSqFitter2d(lsst::cpputils::python::WrapperCollection &wrappers, std::string const &name) {
    using PyLeastSqFitter2d = py::class_<LeastSqFitter2d<FittingFunc>, std::shared_ptr<LeastSqFitter2d<FittingFunc>>>;

    wrappers.wrapType(PyLeastSqFitter2d(wrappers.module, name.c_str()), [](auto &mod, auto &cls) {
        cls.def(py::init<std::vector<double> const &, std::vector<double> const &, std::vector<double> const &,
                        std::vector<double> const &, int>(),
                "x"_a, "y"_a, "z"_a, "s"_a, "order"_a);

        cls.def("getParams", &LeastSqFitter2d<FittingFunc>::getParams);
        cls.def("getErrors", &LeastSqFitter2d<FittingFunc>::getErrors);
        cls.def("valueAt", &LeastSqFitter2d<FittingFunc>::valueAt, "x"_a, "y"_a);
        cls.def("residuals", &LeastSqFitter2d<FittingFunc>::residuals);
        cls.def("getChiSq", &LeastSqFitter2d<FittingFunc>::getChiSq);
        cls.def("getReducedChiSq", &LeastSqFitter2d<FittingFunc>::getReducedChiSq);
    });
}

}  // namespace

void wrapLeastSqFitter2d(lsst::cpputils::python::WrapperCollection &wrappers){
    declareLeastSqFitter2d<afw::math::PolynomialFunction1<double>>(wrappers, "LeastSqFitter2dPoly");
}

}  // namespace sip
}  // namespace astrom
}  // namespace meas
}  // namespace lsst
