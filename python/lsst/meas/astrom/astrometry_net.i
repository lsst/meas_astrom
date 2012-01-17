%define astrometry_net_DOCSTRING
"
Python interface to Astrometry.net
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.astrometry_net",
		docstring=astrometry_net_DOCSTRING) astrometry_net


%{
	extern "C" {
#include "solver.h"
#include "index.h"
#include "fitsioutils.h"
#include "log.h"
#undef FALSE
#undef TRUE
	}
#include <vector>
#include "boost/shared_ptr.hpp"

#include "lsst/base.h"
#include "lsst/pex/logging.h"
#include "lsst/daf/persistence.h"
#include "lsst/daf/base/Persistable.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/DiaSource.h"
#include "lsst/afw/detection/AperturePhotometry.h"

#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom.h"

	namespace afwCoord = lsst::afw::coord;
	namespace afwGeom  = lsst::afw::geom;
	namespace afwImage = lsst::afw::image;

	%}

%init %{
	fits_use_error_system();
	log_init(LOG_MSG);
	%}

%include "std_string.i"
%include "std_vector.i"
%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"
%shared_ptr(lsst::daf::base::Persistable);
%import "lsst/daf/base/Persistable.h"

%lsst_exceptions();

//%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/detection/source.i"

%include "solver.h"
 //%include "log.h"

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

%extend solver_t {

	/*
	 void getTanWcs(double* crpix1, double* crpix2,
	 double* crval1, double* crval2,
	 double* cd1, double* cd2, double* cd3, double* cd4) {
	 }
	 */
	PTR(afwImage::Wcs) getWcs() {
		MatchObj* match = solver_get_best_match($self);
		if (!match)
			return PTR(afwImage::Wcs)();
		tan_t* wcs = &(match->wcstan);

		afwGeom::Point2D crpix(wcs->crpix[0], wcs->crpix[1]);
		afwCoord::Coord::ConstPtr crval
			(new afwCoord::Coord(wcs->crval[0] * afwGeom::degrees,
								 wcs->crval[1] * afwGeom::degrees));
		return afwImage::makeWcs(crval, crpix,
								 wcs->cd[0][0], wcs->cd[0][1],
								 wcs->cd[1][0], wcs->cd[1][1]);
	}

	bool didSolve() {
		return solver_did_solve($self);
	}

	void run() {
		printf("Solver run...\n");
		solver_run($self);
		printf("solver_run returned.\n");
	}

	void addIndex(const char* fn) {
		printf("Loading index file...\n");
		index_t* ind = index_load(fn, INDEX_ONLY_LOAD_METADATA, NULL);
		if (!ind) {
			assert(0);
		}
		printf("  index %i, hp %i (nside %i), nstars %i, nquads %i\n",
			   ind->indexid, ind->healpix, ind->hpnside,
			   ind->nstars, ind->nquads);
		if ($self->use_radec) {
			double ra,dec,radius;
			xyzarr2radecdeg($self->centerxyz, &ra, &dec);
			radius = distsq2deg($self->r2);
			if (!index_is_within_range(ind, ra, dec, radius)) {
				printf("Not within RA,Dec range\n");
				index_free(ind);
				return;
			}
		}
		// qlo,qhi in arcsec
		double qlo, qhi;
		solver_get_quad_size_range_arcsec($self, &qlo, &qhi);
		if (!index_overlaps_scale_range(ind, qlo, qhi)) {
			printf("Not within quad scale range\n");
			index_free(ind);
			return;
		}
		printf("Added index.\n");
		solver_add_index($self, ind);
	}

	void setRaDecRadius(double ra, double dec, double rad) {
		solver_set_radec($self, ra, dec, rad);
	}

	void setImageSize(int W, int H) {
		solver_set_field_bounds($self, 0, W, 0, H);
	}

	void setMaxStars(int N) {
		$self->endobj = N;
	}

	void setStars(lsst::afw::detection::SourceSet srcs) {
		// convert to Astrometry.net "starxy_t"
		const size_t N = srcs.size();
		starxy_t *starxy = starxy_new(N, true, false);
		for (size_t i=0; i<N; ++i) {
			double const x    = srcs[i]->getXAstrom();
			double const y    = srcs[i]->getYAstrom();
			double const flux = srcs[i]->getPsfFlux();
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
