// -*- lsst-C++ -*-
%define astrometry_net_DOCSTRING
"
Python interface to Astrometry.net
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.astrometry_net",
        docstring=astrometry_net_DOCSTRING) astrometry_net

%pythonnondynamic;
%naturalvar;  // use const reference typemaps

%include "std_string.i"
%include "std_vector.i"

%{
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

#include <vector>
#include <set>
#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/format.hpp"

#include "lsst/meas/astrom/detail/utils.h"
#include "lsst/base.h"
#include "lsst/pex/logging.h"
#include "lsst/daf/persistence.h"
#include "lsst/daf/base/Persistable.h"
#include "lsst/daf/base/PropertyList.h"
#include "lsst/afw/table.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/coord.h"

namespace afwCoord = lsst::afw::coord;
namespace afwTable = lsst::afw::table;
namespace afwGeom  = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace dafBase  = lsst::daf::base;
namespace pexLog   = lsst::pex::logging;

struct timer_baton {
    solver_t* s;
    double timelimit;
};

static time_t timer_callback(void* baton) {
    struct timer_baton* tt = static_cast<struct timer_baton*>(baton);
    solver_t* solver = tt->s;
    //printf("Timer callback; time used %f, limit %f\n",
    //       solver->timeused, tt->timelimit);
    if (solver->timeused > tt->timelimit)
        solver->quit_now = 1;
    return 1;
}

// Global logger to which Astrometry.net will go.
static PTR(pexLog::Log) an_log;

PTR(pexLog::Log) get_an_log() {
    return an_log;
}
void set_an_log(PTR(pexLog::Log) newlog) {
    an_log = newlog;
}

static void an_log_callback(void* baton, enum log_level level,
                            const char* file,
                            int line, const char* func, const char* format,
                            va_list va) {
    // translate between logging levels
    int levelmap[5];
    levelmap[LOG_NONE ] = pexLog::Log::FATAL;
    levelmap[LOG_ERROR] = pexLog::Log::FATAL;
    levelmap[LOG_MSG  ] = pexLog::Log::INFO;
    levelmap[LOG_VERB ] = pexLog::Log::DEBUG;
    levelmap[LOG_ALL  ] = pexLog::Log::DEBUG;
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

    dafBase::PropertySet ps;
    ps.set("an_file", file);
    ps.set("an_line", line);
    ps.set("an_func", func);

    an_log->log(lsstlevel, std::string(msg), ps);
}

static void start_an_logging() {
    an_log = PTR(pexLog::Log)(new pexLog::Log(pexLog::Log::getDefaultLog(),
                                              "meas.astrom.astrometry_net"));
    an_log->markPersistent();
    // NOTE, this has to happen before the log_use_function!
    log_init(LOG_VERB);
    log_use_function(an_log_callback, NULL);
    log_to(NULL);
}

static void stop_an_logging() {
    log_use_function(NULL, NULL);
    log_to(stdout);
    an_log.reset();
}

void finalize() {
    stop_an_logging();
}
    %}

%init %{
    // Astrometry.net logging
    fits_use_error_system();
    start_an_logging();
    %}

%include "boost_shared_ptr.i"
%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"
%import "lsst/daf/base/baseLib.i"

void finalize();
PTR(pexLog::Log) get_an_log();
void set_an_log(PTR(pexLog::Log) newlog);

%lsst_exceptions();

%template(VectorOfString) std::vector<std::string>;

%import "lsst/afw/table/Source.i"

%shared_ptr(lsst::afw::image::Wcs);
%import "lsst/afw/image/Wcs.h"

%template(VectorOfIndexPtr) std::vector<index_t*>;

%newobject solver_new;
%newobject multiindex_new;

%include "astrometry/solver.h"
%include "astrometry/index.h"
%include "astrometry/multiindex.h"

%inline %{
    void an_log_init(int level) {
        log_init((log_level)level);
    }
    int an_log_get_level() {
        return (int)log_get_level();
    }
    void an_log_set_level(int lvl) {
        log_set_level((log_level)lvl);
    }

    /// Calculate the distance from coordinates to a healpix
    ///
    /// Note that this assumes that the astrometry.net catalog reference system is ICRS.
    lsst::afw::geom::Angle healpixDistance(int hp, int Nside, lsst::afw::coord::Coord const& coord) {
        lsst::afw::coord::IcrsCoord icrs = coord.toIcrs();
        return lsst::afw::geom::Angle(healpix_distance_to_radec(hp, Nside, icrs.getLongitude().asDegrees(),
                                                                icrs.getLatitude().asDegrees(), NULL),
                                      lsst::afw::geom::degrees);
    }
%}

%extend multiindex_t {
    // An index being within range is a property of the star kd-tree, hence of
    // the multi-index as a whole
    int isWithinRange(double ra, double dec, double radius) {
        return index_is_within_range(multiindex_get($self, 0), ra, dec, radius);
    }

    void unload() {
        multiindex_unload($self);
    }

    ~multiindex_t() {
        multiindex_free($self);
    }

    char* get_name() {
        return $self->fits->filename;
    }

    %pythoncode %{
    def reload(self):
        if multiindex_reload_starkd(self):
            raise RuntimeError('Failed to reload multi-index star file %s' % self.name)
    __swig_getmethods__['name'] = get_name
    if _newclass: x = property(get_name)

    def __len__(self):
        return multiindex_n(self)
    def __getitem__(self, i):
        if i < 0 or i > multiindex_n(self):
            raise IndexError('Index %i out of bound for multiindex_t' % i)
        return multiindex_get(self, i)
    %}
}

%extend index_t {
    int overlapsScaleRange(double qlo, double qhi) {
        return index_overlaps_scale_range($self, qlo, qhi);
    }

    %pythoncode %{
    def reload(self):
        if index_reload(self):
            raise RuntimeError('Failed to reload multi-index file %s' % self.indexname)
    %}
}

%extend solver_t {
    ~solver_t() {
        // Working around a bug in Astrometry.net: doesn't take ownership of the
        // field.
        // unseemly familiarity with the innards... but valgrind-clean.
        starxy_free($self->fieldxy);
        $self->fieldxy = NULL;
        solver_free($self);
    }

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
    lsst::afw::table::SimpleCatalog getCatalog(
        std::vector<index_t*> inds,
        lsst::afw::coord::Coord const &ctrCoord,
        lsst::afw::geom::Angle const &radius,
        const char* idCol,
        std::vector<std::string> const& filterNameList,
        std::vector<std::string> const& magColList,
        std::vector<std::string> const& magErrColList,
        const char* starGalCol,
        const char* varCol,
        bool uniqueIds=true)
    {
        if ((filterNameList.size() != magColList.size()) || (filterNameList.size() != magErrColList.size())) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterError,
                "Filter name, mag column, and mag error column vectors must be the same length.");
        }
        std::vector<lsst::meas::astrom::detail::MagColInfo> magColInfoList;
        for (size_t i=0; i<filterNameList.size(); ++i) {
            lsst::meas::astrom::detail::MagColInfo mc;
            mc.filterName = filterNameList[i];
            mc.magCol = magColList[i];
            mc.magErrCol = magErrColList[i];
            magColInfoList.push_back(mc);
        }
        return lsst::meas::astrom::detail::getCatalogImpl(inds, ctrCoord, radius,
            idCol, magColInfoList, starGalCol, varCol, uniqueIds);
    }

    PTR(lsst::daf::base::PropertyList) getSolveStats() {
        // Gather solve stats...
        PTR(dafBase::PropertyList) qa = boost::make_shared<dafBase::PropertyList>();
        // FIXME -- Ticket #1875 prevents dotted-names from working with toString().
        qa->set("meas_astrom*an*n_tried", $self->numtries);
        qa->set("meas_astrom*an*n_matched", $self->nummatches);
        qa->set("meas_astrom*an*n_scaleok", $self->numscaleok);
        qa->set("meas_astrom*an*n_cxdxcut", $self->num_cxdx_skipped);
        qa->set("meas_astrom*an*n_meanxcut", $self->num_meanx_skipped);
        qa->set("meas_astrom*an*n_radeccut", $self->num_radec_skipped);
        qa->set("meas_astrom*an*n_scalecut", $self->num_abscale_skipped);
        qa->set("meas_astrom*an*n_verified", $self->num_verified);
        qa->set("meas_astrom*an*time_used", $self->timeused);
        qa->set("meas_astrom*an*best_logodds", $self->best_logodds);
        if ($self->best_index) {
            index_t* ind = $self->best_index;
            qa->set("meas_astrom*an*best_index*id", ind->indexid);
            qa->set("meas_astrom*an*best_index*hp", ind->healpix);
            qa->set("meas_astrom*an*best_index*nside", ind->hpnside);
            qa->set("meas_astrom*an*best_index*name", std::string(ind->indexname));
        }
        if ($self->have_best_match) {
            MatchObj* mo = &($self->best_match);
            std::string s = boost::str(boost::format("%i") % mo->star[0]);
            for (int i=1; i<mo->dimquads; i++)
                s = s + boost::str(boost::format(", %i") % mo->star[i]);
            qa->set("meas_astrom*an*best_match*starinds", s);
            qa->set("meas_astrom*an*best_match*coderr", std::sqrt(mo->code_err));
            qa->set("meas_astrom*an*best_match*nmatch", mo->nmatch);
            qa->set("meas_astrom*an*best_match*ndistract", mo->ndistractor);
            qa->set("meas_astrom*an*best_match*nconflict", mo->nconflict);
            qa->set("meas_astrom*an*best_match*nfield", mo->nfield);
            qa->set("meas_astrom*an*best_match*nindex", mo->nindex);
            qa->set("meas_astrom*an*best_match*nbest", mo->nbest);
            qa->set("meas_astrom*an*best_match*logodds", mo->logodds);
            qa->set("meas_astrom*an*best_match*parity", mo->parity ? 0 : 1);
            qa->set("meas_astrom*an*best_match*nobjs", mo->objs_tried);
        }
        return qa;
    }

    PTR(lsst::afw::image::Wcs) getWcs() {
        MatchObj* match = solver_get_best_match($self);
        if (!match)
            return PTR(afwImage::Wcs)();
        tan_t* wcs = &(match->wcstan);

        afwGeom::Point2D crpix(wcs->crpix[0], wcs->crpix[1]);
        afwCoord::Coord::ConstPtr crval
            (new afwCoord::Coord(wcs->crval[0] * afwGeom::degrees,
                                 wcs->crval[1] * afwGeom::degrees));
        return afwImage::makeWcs(*crval, crpix,
                                 wcs->cd[0][0], wcs->cd[0][1],
                                 wcs->cd[1][0], wcs->cd[1][1]);
    }

    bool didSolve() {
        return solver_did_solve($self);
    }

    void run(double cpulimit) {
        solver_log_params($self);
        struct timer_baton tt;
        if (cpulimit > 0.) {
            tt.s = $self;
            tt.timelimit = cpulimit;
            $self->userdata = &tt;
            $self->timer_callback = timer_callback;
        }
        solver_run($self);
        if (cpulimit > 0.) {
            $self->timer_callback = NULL;
            $self->userdata = NULL;
        }
    }

    double getQuadSizeLow() {
        double qlo,qhi;
        solver_get_quad_size_range_arcsec($self, &qlo, &qhi);
        return qlo;
    }
    double getQuadSizeHigh() {
        double qlo,qhi;
        solver_get_quad_size_range_arcsec($self, &qlo, &qhi);
        return qhi;
    }

    void addIndices(std::vector<index_t*> inds) {
        for (std::vector<index_t*>::iterator pind = inds.begin();
             pind != inds.end(); ++pind) {
            lsst::meas::astrom::detail::IndexManager man(*pind);
//            printf("Checking index \"%s\"\n", man.index->indexname);
            if ($self->use_radec) {
                double ra,dec,radius;
                xyzarr2radecdeg($self->centerxyz, &ra, &dec);
                radius = distsq2deg($self->r2);
                if (!index_is_within_range(man.index, ra, dec, radius)) {
                                     //printf("Not within RA,Dec range\n");
                    continue;
                }
            }
            // qlo,qhi in arcsec
            double qlo, qhi;
            solver_get_quad_size_range_arcsec($self, &qlo, &qhi);
            if (!index_overlaps_scale_range(man.index, qlo, qhi)) {
//                printf("Not within quad scale range\n");
                continue;
            }
//            printf("Adding index.\n");
            if (index_reload(man.index)) {
                throw LSST_EXCEPT(lsst::pex::exceptions::IoError,
                                  "Failed to index_reload() an astrometry_net_data index file -- out of file descriptors?");
            }

            solver_add_index($self, man.index);
        }
    }

    void setParity(bool p) {
        if (p)
            $self->parity = PARITY_FLIP;
        else
            $self->parity = PARITY_NORMAL;
    }

    void setMatchThreshold(double t) {
        solver_set_keep_logodds($self, t);
    }

    void setPixelScaleRange(double lo, double hi) {
        $self->funits_lower = lo;
        $self->funits_upper = hi;
    }

    void setRaDecRadius(double ra, double dec, double rad) {
        solver_set_radec($self, ra, dec, rad);
    }

    void setImageSize(int W, int H) {
        solver_set_field_bounds($self, 0, W, 0, H);
        double hi = hypot(W,H);
        double lo = 0.1 * std::min(W,H);
        solver_set_quad_size_range($self, lo, hi);
    }

    void setMaxStars(int N) {
        $self->endobj = N;
    }

    void setStars(lsst::afw::table::SourceCatalog const & srcs, int x0, int y0) {
        // convert to Astrometry.net "starxy_t"
        starxy_free($self->fieldxy);
        const size_t N = srcs.size();
        starxy_t *starxy = starxy_new(N, true, false);
        for (size_t i=0; i<N; ++i) {
            double const x    = srcs[i].getX();
            double const y    = srcs[i].getY();
            double const flux = srcs[i].getPsfFlux();
            starxy_set(starxy, i, x - x0, y - y0);
            starxy_set_flux(starxy, i, flux);
        }
        // Sort the array
        starxy_sort_by_flux(starxy);

        starxy_free(solver_get_field($self));
        solver_free_field($self);
        solver_set_field($self, starxy);
        solver_reset_field_size($self);
        // Find field boundaries and precompute kdtree
        solver_preprocess_field($self);
    }

 }
