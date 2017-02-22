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

#include "lsst/afw/table/fwd.h"
#include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {
namespace sip {

PYBIND11_PLUGIN(matchSrcToCatalogue) {
    py::module mod("matchSrcToCatalogue");

    py::class_<MatchSrcToCatalogue> cls(mod, "MatchSrcToCatalogue");

    cls.def(py::init<afw::table::SimpleCatalog const &, afw::table::SourceCatalog const &,
                     std::shared_ptr<afw::image::Wcs const>, afw::geom::Angle>(),
            "catSet"_a, "imgSet"_a, "wcs"_a, "dist"_a);

    cls.def("setDist", &MatchSrcToCatalogue::setDist, "dist"_a);
    cls.def("setWcs", &MatchSrcToCatalogue::setWcs, "wcs"_a);
    cls.def("setCatSrcSet", &MatchSrcToCatalogue::setCatSrcSet, "catSet"_a);
    cls.def("setImgSrcSet", &MatchSrcToCatalogue::setImgSrcSet, "srcSet"_a);
    cls.def("findMatches", &MatchSrcToCatalogue::findMatches);
    cls.def("getMatches", &MatchSrcToCatalogue::getMatches);

    return mod.ptr();
}
}
}
}
}  // namespace lsst::meas::astrom::sip
