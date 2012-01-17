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
#include "starkd.h"
#include "fitsioutils.h"
#include "fitstable.h"
#include "log.h"
#undef FALSE
#undef TRUE
	}
#include <vector>
#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"

#include "lsst/base.h"
#include "lsst/pex/logging.h"
#include "lsst/daf/persistence.h"
#include "lsst/daf/base/Persistable.h"
#include "lsst/daf/base/PropertyList.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/DiaSource.h"
#include "lsst/afw/detection/AperturePhotometry.h"

#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom.h"

	namespace afwCoord = lsst::afw::coord;
	namespace afwDet   = lsst::afw::detection;
	namespace afwGeom  = lsst::afw::geom;
	namespace afwImage = lsst::afw::image;
	namespace dafBase  = lsst::daf::base;

	static afwCoord::Coord::Ptr radectocoord(afwCoord::CoordSystem coordsys,
											 const double* radec) {
		afwCoord::Coord::Ptr rd = afwCoord::makeCoord(coordsys,
													  radec[0] * afwGeom::degrees,
													  radec[1] * afwGeom::degrees);
		return rd;
	}


	%}

%init %{
	fits_use_error_system();
	log_init(LOG_MSG);
	%}

%include "std_string.i"
%include "std_vector.i"
%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"
%import "lsst/daf/base/baseLib.i"

%lsst_exceptions();

//%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/detection/source.i"

%shared_ptr(lsst::afw::image::Wcs);
%import "lsst/afw/image/Wcs.h"

%template(VectorOfIndexPtr) std::vector<index_t*>;
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

%extend solver_t {

	lsst::afw::detection::SourceSet
		getCatalog(std::vector<index_t*> inds,
				   double ra, double dec, double radius,
				   const char* magcol,
				   const char* magerrcol,
				   const char* stargalcol,
				   const char* varcol,
				   boost::int64_t starflag) {
		lsst::afw::detection::SourceSet cat;

        // FIXME -- cut on indexid?
		// FIXME -- cut on healpix?

		// additional margin on healpixes, in deg.
		double margin = 1.0;

		double xyz[3];
		radecdeg2xyzarr(ra, dec, xyz);
		double r2 = deg2distsq(radius);

		for (std::vector<index_t*>::iterator pind = inds.begin();
			 pind != inds.end(); ++pind) {
			index_t* ind = (*pind);
			printf("checking index \"%s\"\n", ind->indexname);
			if (!index_is_within_range(ind, ra, dec, radius + margin)) {
				printf(" skipping: not within range\n");
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
			printf("found %i\n", nstars);
			if (nstars == 0)
				continue;
			// FIXME -- handle duplicates here, or in python?

			float* mag = NULL;
			float* magerr = NULL;
			bool* stargal = NULL;
			bool* var = NULL;
			if (magcol || magerrcol || stargalcol || varcol) {
				fitstable_t* tag = startree_get_tagalong(ind->starkd);
				tfits_type flt = fitscolumn_float_type();
				tfits_type boo = fitscolumn_boolean_type();

				if (magcol) {
					//mag = startree_get_data_column(ind->starkd, magcol, starinds, nstars);
					mag = static_cast<float*>(fitstable_read_column_inds(tag, magcol, flt, starinds, nstars));
					assert(mag);
				}
				if (magerrcol) {
					//mag = startree_get_data_column(ind->starkd, magerrcol, starinds, nstars);
					magerr = static_cast<float*>(fitstable_read_column_inds(tag, magerrcol, flt, starinds, nstars));
					assert(magerr);
				}
				if (stargalcol) {
					stargal = static_cast<bool*>(fitstable_read_column_inds(tag, stargalcol, boo, starinds, nstars));
					assert(stargal);
				}
				if (varcol) {
					var = static_cast<bool*>(fitstable_read_column_inds(tag, varcol, boo, starinds, nstars));
					assert(var);
				}
			}

			// FIXME -- allow this to be set?
			afwCoord::CoordSystem coordsys = afwCoord::ICRS;

			for (int i=0; i<nstars; i++) {
				afwDet::Source::Ptr src(new afwDet::Source());
				src->setAllRaDecFields(radectocoord(coordsys, radecs + i*2));

				if (mag) {
					// LAME!
					double flux = pow(10.0, -mag[i]/2.5);
					src->setPsfFlux(flux);
					if (magerr)
						src->setPsfFluxErr(magerr[i] * flux * -std::log(10.)/2.5);
				}
				if (starflag) {
					bool ok = true;
					if (stargal)
						ok &= stargal[i];
					if (var)
						ok &= (!var[i]);
					if (ok)
						src->setFlagForDetection(src->getFlagForDetection() | starflag);
				}

				cat.push_back(src);
			}

			free(mag);
			free(magerr);
			free(stargal);
			free(var);

			free(radecs);
			free(starinds);
		}

		return cat;
	}

	PTR(lsst::daf::base::PropertyList) getSolveStats() {
		// Gather solve stats...
		PTR(dafBase::PropertyList) qa(new dafBase::PropertyList);
		std::string s = qa->toString();
		printf("empty QA: %s\n", s.c_str());
		qa->set("meas_astrom.an.n_tried", $self->numtries);
		qa->set("meas_astrom.an.n_matched", $self->nummatches);
		qa->set("meas_astrom.an.n_scaleok", $self->numscaleok);
		qa->set("meas_astrom.an.n_cxdxcut", $self->num_cxdx_skipped);
		qa->set("meas_astrom.an.n_meanxcut", $self->num_meanx_skipped);
		qa->set("meas_astrom.an.n_radeccut", $self->num_radec_skipped);
		qa->set("meas_astrom.an.n_scalecut", $self->num_abscale_skipped);
		qa->set("meas_astrom.an.n_verified", $self->num_verified);
		qa->set("meas_astrom.an.time_used", $self->timeused);
		qa->set("meas_astrom.an.best_logodds", $self->best_logodds);
		if ($self->best_index) {
			index_t* ind = $self->best_index;
			qa->set("meas_astrom.an.best_index.id", ind->indexid);
			qa->set("meas_astrom.an.best_index.hp", ind->healpix);
			qa->set("meas_astrom.an.best_index.nside", ind->hpnside);
			qa->set("meas_astrom.an.best_index.name", std::string(ind->indexname));
		}
		s = qa->toString();
		printf("returning QA: %s\n", s.c_str());
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

	void addIndices(const std::vector<index_t*> & inds) {
		for (std::vector<index_t*>::iterator pind = inds.begin();
			 pind != inds.end(); ++pind) {
			index_t* ind = *pind;
			printf("Checking index \"%s\"\n", ind->indexname);
			if ($self->use_radec) {
				double ra,dec,radius;
				xyzarr2radecdeg($self->centerxyz, &ra, &dec);
				radius = distsq2deg($self->r2);
				if (!index_is_within_range(ind, ra, dec, radius)) {
					printf("Not within RA,Dec range\n");
					continue;
				}
			}
			// qlo,qhi in arcsec
			double qlo, qhi;
			solver_get_quad_size_range_arcsec($self, &qlo, &qhi);
			if (!index_overlaps_scale_range(ind, qlo, qhi)) {
				printf("Not within quad scale range\n");
				continue;
			}
			printf("Adding index.\n");
			if (index_reload(ind)) {
				assert(0);
			}
			solver_add_index($self, ind);
		}
	}

	/*
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
		if (index_reload(ind)) {
			assert(0);
		}
		solver_add_index($self, ind);
	}
	 */

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
