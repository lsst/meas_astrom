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
#include "solver.h"
#include "index.h"
#include "starkd.h"
#include "fitsioutils.h"
#include "fitstable.h"
#include "log.h"
#include "tic.h"

#undef ATTRIB_FORMAT
#undef FALSE
#undef TRUE
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
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/geom.h"

namespace afwCoord = lsst::afw::coord;
namespace afwTable = lsst::afw::table;
namespace afwGeom  = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace dafBase  = lsst::daf::base;

struct timer_baton {
    solver_t* s;
    double timelimit;
};

static time_t timer_callback(void* baton) {
    struct timer_baton* tt = static_cast<struct timer_baton*>(baton);
    solver_t* solver = tt->s;

    // Unfortunately, a bug in Astrometry.net 0.30 means the timeused is not updated before
    // calling the timer callback.  Doh!

    // solver.c : update_timeused (which is static) does:
    double usertime, systime;
    get_resource_stats(&usertime, &systime, NULL);
    solver->timeused = std::max(0., (usertime + systime) - solver->starttime);
    //printf("Timer callback; time used %f, limit %f\n", solver->timeused, tt->timelimit);
    if (solver->timeused > tt->timelimit)
        solver->quit_now = 1;
    return 1;
}
    %}

%init %{
    // Astrometry.net logging
    fits_use_error_system();
    log_init(LOG_MSG);
    %}

%include "boost_shared_ptr.i"
%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"
%import "lsst/daf/base/baseLib.i"

%lsst_exceptions();

%template(VectorOfString) std::vector<std::string>;

%import "lsst/afw/table/tableLib.i"

%shared_ptr(lsst::afw::image::Wcs);
%import "lsst/afw/image/Wcs.h"

%template(VectorOfIndexPtr) std::vector<index_t*>;
%newobject solver_new;
%newobject index_load;
%include "solver.h"
%include "index.h"

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

// index_t is a typedef of index_s, but swig doesn't notice the typedef, grumble
// grumble.
%extend index_s {
    ~index_s() {
        //printf("Deleting index_s %s\n", $self->indexname);
        index_free($self);
    }
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
//        printf("Solver run...\n");
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
//        printf("solver_run returned.\n");
    }

    std::vector<index_t*> getActiveIndexFiles() {
        std::vector<index_t*> inds;
        int N = solver_n_indices($self);
        for (int i=0; i<N; i++) {
            inds.push_back(solver_get_index($self, i));
        }
        return inds;
    }

    void addIndices(std::vector<index_t*> inds) {
        for (std::vector<index_t*>::iterator pind = inds.begin();
             pind != inds.end(); ++pind) {
            index_t* ind = *pind;
//            printf("Checking index \"%s\"\n", ind->indexname);
            if ($self->use_radec) {
                double ra,dec,radius;
                xyzarr2radecdeg($self->centerxyz, &ra, &dec);
                radius = distsq2deg($self->r2);
                if (!index_is_within_range(ind, ra, dec, radius)) {
                                     //printf("Not within RA,Dec range\n");
                    continue;
                }
            }
            // qlo,qhi in arcsec
            double qlo, qhi;
            solver_get_quad_size_range_arcsec($self, &qlo, &qhi);
            if (!index_overlaps_scale_range(ind, qlo, qhi)) {
//                printf("Not within quad scale range\n");
                continue;
            }
//            printf("Adding index.\n");
            if (index_reload(ind)) {
                throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException,
                                  "Failed to index_reload() an astrometry_net_data index file -- out of file descriptors?");
            }

            // Change once astrometry.net-0.40+ is in...
            if (ind->quads->fb->fid) {
                if (fclose(ind->quads->fb->fid)) {
                    //SYSERROR("Failed to fclose() astrometry_net_data quadfile");
                    throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException,
                                      "Failed to fclose() an astrometry_net_data quadfile");
                }
                ind->quads->fb->fid = NULL;
            }

            kdtree_fits_t* io;
            io = (kdtree_fits_t*)ind->codekd->tree->io;
            if (io->fid) {
                if (fclose(io->fid)) {
                    throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException,
                                      "Failed to fclose() an astrometry_net_data code kdtree");
                }
                io->fid = NULL;
            }
                
            io = (kdtree_fits_t*)ind->starkd->tree->io;
            if (io->fid) {
                if (fclose(io->fid)) {
                    throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException,
                                      "Failed to fclose() an astrometry_net_data star kdtree");
                }
                io->fid = NULL;
            }

            solver_add_index($self, ind);
        }
    }

    void setParity(bool p) {
        if (p)
            $self->parity = PARITY_FLIP;
        else
            $self->parity = PARITY_NORMAL;
    }

    void setMatchThreshold(double t) {
        // AN 0.30:
        $self->logratio_record_threshold = t;
        // AN 0.40:
        //solver_set_keep_logodds($self, t)
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
