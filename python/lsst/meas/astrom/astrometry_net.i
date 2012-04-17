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
#undef FALSE
#undef TRUE
    }

#include <vector>
#include <set>
#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/format.hpp"

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

/*
 * Implementation for index_s::getCatalog method
 */
static lsst::afw::table::SimpleCatalog
getCatalogImpl(std::vector<index_t*> inds,
	       double ra, double dec, double radius,
	       const char* idcol,
	       std::vector<std::string> const & magcolVec,
	       std::vector<std::string> const & magerrcolVec,
	       const char* stargalcol,
	       const char* varcol,
	       bool unique_ids=true)
{
   assert(!magcolVec.empty());
   char const* magcol = magcolVec[0] == "" ? NULL : magcolVec[0].c_str();
   assert(!magerrcolVec.empty());
   char const* magerrcol = magerrcolVec[0] == "" ? NULL : magerrcolVec[0].c_str();

   /*
     If unique_ids == true: return only reference sources with unique IDs;
     arbitrarily keep the first star found with each ID.
   */

   // FIXME -- cut on indexid?
   // FIXME -- cut on healpix?

   // additional margin on healpixes, in deg.
   double margin = 1.0;

   double xyz[3];
   radecdeg2xyzarr(ra, dec, xyz);
   double r2 = deg2distsq(radius);

   afwTable::Schema schema = afwTable::SimpleTable::makeMinimalSchema(); // contains ID, ra, dec.
   afwTable::Key<double> fluxKey;    // these are double for consistency with measured fluxes;
   afwTable::Key<double> fluxErrKey; // may be unnecessary, but less surprising.
   if (magcol) {
      fluxKey = schema.addField<double>("flux", "flux");
      if (magerrcol) {
	 fluxErrKey = schema.addField<double>("flux.err", "flux uncertainty");
      }
   }
   afwTable::Key<afwTable::Flag> stargalKey;
   if (stargalcol) {
      stargalKey = schema.addField<afwTable::Flag>(
	 "stargal", "set if the reference object is a star"
	 );
   }
   afwTable::Key<afwTable::Flag> varKey;
   if (varcol) {
      varKey = schema.addField<afwTable::Flag>("var", "set if the reference object is variable");
   }
   afwTable::Key<afwTable::Flag> photometricKey = schema.addField<afwTable::Flag>(
      "photometric", "set if the reference object can be used in photometric calibration"
      );
            
   // make catalog with no IdFactory, since IDs are external
   afwTable::SimpleCatalog cat(afwTable::SimpleTable::make(schema, PTR(afwTable::IdFactory)()));

   // for unique_ids: keep track of the IDs we have already added to the result set.
   std::set<boost::int64_t> uids;

   for (std::vector<index_t*>::iterator pind = inds.begin();
	pind != inds.end(); ++pind) {
      index_t* ind = (*pind);
      //printf("checking index \"%s\"\n", ind->indexname);
      if (!index_is_within_range(ind, ra, dec, radius + margin)) {
	 //printf(" skipping: not within range\n");
	 continue;
      }
      // Ensure the index is loaded...
      index_reload(ind);

      // Find nearby stars
      double *radecs = NULL;
      int *starinds = NULL;
      int nstars = 0;
      startree_search_for(ind->starkd, xyz, r2, NULL,
			  &radecs, &starinds, &nstars);
      //printf("found %i\n", nstars);
      if (nstars == 0)
	 continue;

      float* mag = NULL;
      float* magerr = NULL;
      boost::int64_t* id = NULL;
      bool* stargal = NULL;
      bool* var = NULL;
      if (idcol || magcol || magerrcol || stargalcol || varcol) {
	 fitstable_t* tag = startree_get_tagalong(ind->starkd);
	 tfits_type flt = fitscolumn_float_type();
	 tfits_type boo = fitscolumn_boolean_type();
	 tfits_type i64 = fitscolumn_i64_type();

	 if (idcol) {
	    id = static_cast<boost::int64_t*>(fitstable_read_column_inds(tag, idcol, i64, starinds, nstars));
	    assert(id);
	 }

	 if (id && unique_ids) {
	    // remove duplicate IDs.

	    // FIXME -- this shouldn't be necessary once we get astrometry_net 0.40
	    // multi-index functionality in place.

	    if (uids.empty()) {
	       uids = std::set<boost::int64_t>(id, id+nstars);
	    } else {
	       int nkeep = 0;
	       for (int i=0; i<nstars; i++) {
		  //std::pair<std::set<boost::int64_t>::iterator, bool> 
		  if (uids.insert(id[i]).second) {
		     // inserted; keep this one.
		     if (nkeep != i) {
			// compact the arrays.
			starinds[nkeep] = starinds[i];
			radecs[nkeep*2+0] = radecs[i*2+0];
			radecs[nkeep*2+1] = radecs[i*2+1];
			id[nkeep] = id[i];
		     }
		     nkeep++;
		  } else {
		     // did not insert (this id has already been found);
		     // drop this star.
		  }
	       }
	       nstars = nkeep;
	       // if they were all duplicate IDs...
	       if (nstars == 0) {
		  free(starinds);
		  free(radecs);
		  free(id);
		  continue;
	       }
	    }
	 }

	 if (magcol) {
	    mag = static_cast<float*>(fitstable_read_column_inds(tag, magcol, flt, starinds, nstars));
	    assert(mag);
	 }
	 if (magerrcol) {
	    magerr = static_cast<float*>(fitstable_read_column_inds(tag, magerrcol, flt, starinds, nstars));
	    assert(magerr);
	 }
	 if (stargalcol) {
	    /*  There is something weird going on with handling of bools; maybe "T" vs "F"?
		stargal = static_cast<bool*>(fitstable_read_column_inds(tag, stargalcol, boo, starinds, nstars));
		for (int j=0; j<nstars; j++) {
		printf("  sg %i = %i, vs %i\n", j, (int)sg[j], stargal[j] ? 1:0);
		}
	    */
	    uint8_t* sg = static_cast<uint8_t*>(fitstable_read_column_inds(tag, stargalcol, fitscolumn_u8_type(), starinds, nstars));
	    stargal = static_cast<bool*>(malloc(nstars));
	    assert(stargal);
	    for (int j=0; j<nstars; j++) {
	       stargal[j] = (sg[j] > 0);
	    }
	    free(sg);
	 }
	 if (varcol) {
	    var = static_cast<bool*>(fitstable_read_column_inds(tag, varcol, boo, starinds, nstars));
	    assert(var);
	 }
      }

      for (int i=0; i<nstars; i++) {
	 PTR(afwTable::SimpleRecord) src = cat.addNew();

	 // Note that all coords in afwTable catalogs are ICRS; hopefully that's what the 
	 // reference catalogs are (and that's what the code assumed before JFB modified it).
	 src->setCoord(
	    afwCoord::IcrsCoord(
	       radecs[i * 2 + 0] * afwGeom::degrees,
	       radecs[i * 2 + 1] * afwGeom::degrees
	       )
	    );

	 if (id)    src->setId(id[i]);

	 if (mag) {
	    // Dustin thinks converting to flux is 'LAME!';
	    // Jim thinks it's nice for consistency (and photocal wants fluxes
	    // as inputs, so we'll continue to go with that for now) even though
	    // we don't need to anymore.
	    // Dustin rebuts that photocal immediately converts those fluxes into mags,
	    // so :-P
	    double flux = pow(10.0, -mag[i] / 2.5);
	    src->set(fluxKey, flux);
	    if (magerr)    src->set(fluxErrKey, magerr[i] * flux * -std::log(10.0) / 2.5);
	 }

	 bool ok = true;
	 if (stargal) {
	    src->set(stargalKey, stargal[i]);
	    ok &= stargal[i];
	 }
	 if (var) {
	    src->set(varKey, var[i]);
	    ok &= (!var[i]);
	 }
	 src->set(photometricKey, ok);
      }

      free(id);
      free(mag);
      free(magerr);
      free(stargal);
      free(var);
      free(radecs);
      free(starinds);
   }
   return cat;
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
/*
 %typemap(newfree) solver_t* {
 printf("I'm the newfree typemap for solver_t\n");
 };
 */
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

/* swig doesn't notice the typedef? grumble grumble.
 %extend index_t {
 ~index_t() {
 printf("Deleting index_t %s\n", $self->indexname);
 index_free($self);
 }
 }
 */

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
		  std::vector<std::string> const& magcolVec,
		  std::vector<std::string> const& magerrcolVec,
		  const char* stargalcol,
		  const char* varcol,
		  bool unique_ids=true)
    {
       return getCatalogImpl(inds, ra, dec, radius,
			     idcol, magcolVec, magerrcolVec, stargalcol, varcol, unique_ids);
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
       std::vector<std::string> magcolVec;
       std::vector<std::string> magerrcolVec;

       magcolVec.push_back(std::string(magcol ? magcol : ""));
       magerrcolVec.push_back(std::string(magerrcol ? magerrcol : ""));

       return getCatalogImpl(inds, ra, dec, radius,
			     idcol, magcolVec, magerrcolVec, stargalcol, varcol, unique_ids);
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
                assert(0);
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

    void setStars(lsst::afw::table::SourceCatalog const & srcs) {
        // convert to Astrometry.net "starxy_t"
        starxy_free($self->fieldxy);
        const size_t N = srcs.size();
        starxy_t *starxy = starxy_new(N, true, false);
        for (size_t i=0; i<N; ++i) {
            double const x    = srcs[i].getX();
            double const y    = srcs[i].getY();
            double const flux = srcs[i].getPsfFlux();
            starxy_set(starxy, i, x, y);
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
