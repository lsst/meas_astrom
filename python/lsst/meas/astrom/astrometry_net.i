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
        solver_free($self);
    }

    lsst::afw::table::SimpleCatalog
       getCatalog(std::vector<index_t*> inds,
                  double ra, double dec, double radius,
                  const char* idcol,
                  std::vector<std::string> const& magnameVec,
                  std::vector<std::string> const& magcolVec,
                  std::vector<std::string> const& magerrcolVec,
                  const char* stargalcol,
                  const char* varcol,
                  bool unique_ids=true)
    {
        if ((magnameVec.size() != magcolVec.size()) || (magnameVec.size() != magerrcolVec.size())) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                              "Mag name, mag column, and mag error column vectors must be the same length.");
        }
        std::vector<lsst::meas::astrom::detail::mag_column_t> magcols;
        for (size_t i=0; i<magnameVec.size(); ++i) {
            lsst::meas::astrom::detail::mag_column_t mc;
            mc.name = magnameVec[i];
            mc.magcol = magcolVec[i];
            mc.magerrcol = magerrcolVec[i];
            magcols.push_back(mc);
        }
        return lsst::meas::astrom::detail::getCatalogImpl(inds, ra, dec, radius,
                                                          idcol, magcols, stargalcol, varcol,
                                                          unique_ids);
    }

    lsst::afw::table::SimpleCatalog
        getCatalog(std::vector<index_t*> inds,
                   double ra, double dec, double radius,
                   const char* idcol,
                   const char* magcol,
                   const char* magerrcol,
                   const char* stargalcol,
                   const char* varcol,
                   bool unique_ids=true)
    {
        std::vector<lsst::meas::astrom::detail::mag_column_t> cols;

        if (magcol || magerrcol) {
            lsst::meas::astrom::detail::mag_column_t mc;
            mc.name = (magcol ? magcol : magerrcol);
            mc.magcol = (magcol ? magcol : "");
            mc.magerrcol = (magerrcol ? magerrcol : "");
            cols.push_back(mc);

            // Also add plain old "flux"/"flux.err" schema entries.
            if (magcol) {
                mc.name = "flux";
                mc.magcol = magcol;
                mc.magerrcol = (magerrcol ? magerrcol : "");
                cols.push_back(mc);
            }
        }
        return lsst::meas::astrom::detail::getCatalogImpl(inds, ra, dec, radius,
                                                          idcol, cols, stargalcol, varcol,
                                                          unique_ids);
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

    void addIndex(index_t* ind) {
        solver_add_index($self, ind);
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
