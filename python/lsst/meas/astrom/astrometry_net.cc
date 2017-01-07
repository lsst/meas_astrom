// -*- lsst-C++ -*-

// Astrometry.net include files...
extern "C" {
#include "astrometry/solver.h"
#include "astrometry/index.h"
#include "astrometry/multiindex.h"
#include "astrometry/starkd.h"
#include "astrometry/fitsioutils.h"
#include "astrometry/fitstable.h"
#include "astrometry/log.h"
#include "astrometry/tic.h"
#include "astrometry/healpix.h"

#undef ATTRIB_FORMAT
#undef FALSE
#undef TRUE

#undef logdebug
#undef debug
}

#include <string>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lsst/log/Log.h"
#include "lsst/utils/pybind11.h"
#include "lsst/meas/astrom/astrometry_net.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace astrom {

namespace {

/*
 * Wrap index_t
 *
 * This class has no constructor; the wrapper is only used for objects returned by MultiIndex::operator[]
 */
void declareIndex(py::module & mod) {
    py::class_<index_t> cls(mod, "index_t");

    cls.def("overlapsScaleRange", [](index_t & self, double qlo, double qhi) {
        return index_overlaps_scale_range(&self, qlo, qhi);
    }, "qlow"_a, "qhigh"_a);
    cls.def("reload", [](index_t & self) {
        if(index_reload(&self)) {
            std::ostringstream os;
            os << "Failed to reload multi-index file " << self.indexname;
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, os.str());
        }
    });

    cls.def_readonly("indexid", &index_t::indexid);
    cls.def_readonly("healpix", &index_t::healpix);
    cls.def_readonly("hpnside", &index_t::hpnside);
    cls.def_readonly("nstars", &index_t::nstars);
    cls.def_readonly("nquads", &index_t::nquads);
}

/**
 * Wrap MultiIndex, a thin shim around multiindex_t
 */
void declareMultiIndex(py::module & mod) {
    py::class_<MultiIndex> cls(mod, "MultiIndex");

    cls.def(py::init<std::string const &>(), "filepath"_a);

    cls.def("__getitem__", [](MultiIndex const & self, int i) {
        auto cind = lsst::utils::cppIndex(self.getLength(), i);
        return self[cind];
    }, py::return_value_policy::reference_internal, py::is_operator());

    cls.def("addIndex", &MultiIndex::addIndex, "filepath"_a, "metadataOnly"_a);
    cls.def("isWithinRange", &MultiIndex::isWithinRange, "ra"_a, "dec"_a, "radius"_a);
    cls.def("unload", &MultiIndex::unload);
    cls.def_property_readonly("name", &MultiIndex::getName);
    cls.def("__len__", &MultiIndex::getLength);
    cls.def("reload", &MultiIndex::reload);
}

/**
 * Wrap Solver, a thin shim around solver_t
 */
void declareSolver(py::module & mod) {
    py::class_<Solver> cls(mod, "Solver");

    cls.def(py::init<>());

    /**
    Load reference objects in a region of the sky described by a center coordinate and a radius

    @param[in] inds  list of star kd-trees from astrometry.net
    @param[in] ctrCoord  center of search region
    @param[in] radius  search radius
    @param[in] idCol  name of ID column in astrometry.net data
    @param[in] filterNameList  names of filters in astrometry.net data
    @param[in] magColList  names of magnitude columns in astrometry.net data
    @param[in] magErrColList  names of magnitude uncertainty (sigma) columns in astrometry.net data
    @param[in] starGalCol  name of "starGal" column (true if object is a star) in astrometry.net data
    @param[in] varCol  name of "var" column (true if brightness is variable) in astrometry.net data
    @param[in] uniqueIds  if true then only return unique IDs (the first of each seen)

    Returned schema:
    - id
    - coord: sky position (an lsst::afw::coord::IcrsCoord)
    - centroid: centroid on some exposure, if relevant (an lsst::afw::geom::Point2D); returned value is not set
    - hasCentroid: if true then centroid has been set; returned value is false
    - <filterName>_flux: flux in the specified filter (double)
    - <filterName>_fluxSigma: flux uncertainty in the specified filter (double)
    - resolved (if starGalCol specified): true if object is not resolved
    - variable (if varCol specified): true if brightness is variable
    - photometric: true if not resolved (or starGalCol blank) and not variable (or varCol blank);
        note that if starGalCol and varCol both blank then all objects are claimed to be photometric
    */
    cls.def("getCatalog", &Solver::getCatalog,
            "inds"_a, "ctrCoord"_a, "radius"_a, "idCol"_a, "filterNameList"_a, "magColList"_a,
            "magErrColList"_a, "starGalCol"_a, "varCol"_a, "uniqueIds"_a=true);
    cls.def("getSolveStats", &Solver::getSolveStats);
    cls.def("getWcs", &Solver::getWcs);
    cls.def("didSolve", &Solver::didSolve);
    cls.def("run", &Solver::run, "cpulimit"_a);
    cls.def("getQuadSizeRangeArcsec", &Solver::getQuadSizeRangeArcsec);
    cls.def("addIndices", &Solver::addIndices, "indices"_a);
    cls.def("setParity", &Solver::setParity, "setParityFlipped", "parity"_a);
    cls.def("setMatchThreshold", &Solver::setMatchThreshold, "threshold"_a);
    cls.def("setPixelScaleRange", &Solver::setPixelScaleRange, "low"_a, "high"_a);
    cls.def("setRaDecRadius", &Solver::setRaDecRadius, "ra"_a, "dec"_a, "rad"_a);
    cls.def("setImageSize", &Solver::setImageSize, "width"_a, "height"_a);
    cls.def("setMaxStars", &Solver::setMaxStars, "maxStars"_a);
    cls.def("setStars", &Solver::setStars, "sourceCat"_a, "x0"_a, "y0"_a);
}

// declare logging functions for use by the Python

LOG_LOGGER an_log = LOG_GET("meas.astrom.astrometry_net");

void an_log_callback(void* baton, enum log_level level,
                            const char* file,
                            int line, const char* func, const char* format,
                            va_list va) {
    // translate between logging levels
    int levelmap[5];
    levelmap[LOG_NONE ] = LOG_LVL_FATAL;
    levelmap[LOG_ERROR] = LOG_LVL_FATAL;
    levelmap[LOG_MSG  ] = LOG_LVL_INFO;
    levelmap[LOG_VERB ] = LOG_LVL_DEBUG;
    levelmap[LOG_ALL  ] = LOG_LVL_DEBUG;
    int lsstlevel = levelmap[level];

    va_list vb;
    // find out how long the formatted string will be
    va_copy(vb, va);
    const int len = vsnprintf(NULL, 0, format, va) + 1; // "+ 1" for the '\0'
    va_end(va);
    // allocate a string of the appropriate length
    char msg[len];
    (void)vsnprintf(msg, len, format, vb);
    va_end(vb);

    // trim trailing \n
    if (msg[len-2] == '\n') {
        msg[len-2] = '\0';
    }

    // not using the LOGS macro because the original location info is wanted
    if (an_log.isEnabledFor(lsstlevel)) {
        an_log.logMsg(log4cxx::Level::toLevel(lsstlevel),
                       log4cxx::spi::LocationInfo(file, func, line), msg);
    }
}

/// start astrometry_net logging
void start_an_logging() {
    // NOTE, this has to happen before the log_use_function!
    log_init(LOG_VERB);
    log_use_function(an_log_callback, NULL);
    log_to(NULL);
}

/// stop astrometry_net logging and perform any other necessary cleanup
void finalize() {
    log_use_function(NULL, NULL);
    log_to(stdout);
}

}  // namespace lsst::meas::astrom::<anonymous>

PYBIND11_PLUGIN(_astrometry_net) {
    py::module mod("_astrometry_net",
                   "Python wrapper for astrometry_net code needed by meas_astrom");

    // code that is run at import time
    fits_use_error_system();
    start_an_logging();

    mod.def("healpixDistance", &healpixDistance, "hp"_a, "nside"_a, "coord"_a);

    mod.def("an_log_init", [](int level) {
        log_init(static_cast<log_level>(level));
    }, "level"_a);

    mod.def("an_log_get_level", []() {
        return static_cast<int>(log_get_level());
    });
    mod.def("an_log_set_level", [](int level) {
        log_set_level(static_cast<log_level>(level));
    }, "level"_a);
    mod.def("finalize", &finalize);

    declareMultiIndex(mod);
    declareIndex(mod);
    declareSolver(mod);

    return mod.ptr();
}

}}}  // namespace lsst::meas::astrom