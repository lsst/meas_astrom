// -*- lsst-C++ -*-
// Astrometry.net include files...
extern "C" {
#include "solver.h"
#include "index.h"
#include "starkd.h"
#include "fitsioutils.h"
#include "fitstable.h"

#undef ATTRIB_FORMAT
#undef FALSE
#undef TRUE
}

#include <set>
#include "boost/cstdint.hpp"
#include "boost/format.hpp"

#include "lsst/meas/astrom/detail/utils.h"
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/Source.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "lsst/afw/geom/Angle.h"
#pragma clang diagnostic pop

namespace afwCoord = lsst::afw::coord;
namespace afwTable = lsst::afw::table;
namespace afwGeom  = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace astrom {
namespace detail {

static float *
read_column(fitstable_t* tag,
            char const *colName,
            tfits_type type,
            int *starinds,
            int nstars,
            char const *indexName)
{
    float *col =
        static_cast<float*>(fitstable_read_column_inds(tag, colName, type, starinds, nstars));
    if (!col) {
        throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                          str(boost::format("Unable to read data for %s from %s") % colName % indexName));
    }

    return col;
}

    
/*
 * Implementation for index_s::getCatalog method
 */
afwTable::SimpleCatalog
getCatalogImpl(std::vector<index_t*> inds,
	       double ra, double dec, double radius,
	       char const* idcol,
	       std::vector<std::string> const & magcolVec,
	       std::vector<std::string> const & magerrcolVec,
	       char const* stargalcol,
	       char const* varcol,
	       bool unique_ids)
{
   unsigned int const nMag = magcolVec.size();       /* number of magnitude columns */
   unsigned int const nMagErr = magerrcolVec.size(); /* number of magnitude error columns */
   if (nMagErr > nMag) {
       throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
			"You may not specify more mag errors than mags");
   }
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
   std::vector<afwTable::Key<double> > fluxKey;	// these are double for consistency with measured fluxes;
   fluxKey.reserve(1 + nMag);
   std::vector<afwTable::Key<double> > fluxErrKey; // may be unnecessary, but less surprising.
   fluxErrKey.reserve(1 + nMagErr);

   if (nMag) {
      fluxKey.push_back(schema.addField<double>("flux", "flux"));
      if (nMagErr) {
          fluxErrKey.push_back(schema.addField<double>("flux.err", "flux uncertainty"));
      }
   }
   for (unsigned int j = 0; j != nMag; ++j) {
      std::string const& magCol = magcolVec[j];
      fluxKey.push_back(schema.addField<double>(magCol, magCol + std::string(" flux")));
      if (j < nMagErr) {
          fluxErrKey.push_back(schema.addField<double>(magCol + ".err", magCol + " flux uncertainty"));
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

   afwTable::SimpleCatalog cat;
   if (idcol) {
       // make catalog with no IdFactory, since IDs are external
       cat = afwTable::SimpleCatalog(afwTable::SimpleTable::make(schema, PTR(afwTable::IdFactory)()));
   } else {
       // let the catalog assign IDs
       cat = afwTable::SimpleCatalog(afwTable::SimpleTable::make(schema));
   }

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

      std::vector<float*> mag;
      mag.reserve(nMag + 1);
      std::vector<float*> magerr;
      magerr.reserve(nMagErr + 1);
      boost::int64_t* id = NULL;
      bool* stargal = NULL;
      bool* var = NULL;
      if (idcol || nMag || nMagErr || stargalcol || varcol) {
	 fitstable_t* tag = startree_get_tagalong(ind->starkd);
	 tfits_type flt = fitscolumn_float_type();
	 tfits_type boo = fitscolumn_boolean_type();
	 tfits_type i64 = fitscolumn_i64_type();

         if (!tag) {
             std::string msg = boost::str(boost::format("astrometry_net_data index file %s does not contain a tag-along table, "
                                                        "so can't retrieve extra columns.  idcol=%s, stargalcol=%s, varcol=%s") %
                                          ind->indexname % idcol % stargalcol % varcol);
             msg += ", mags=[";
             for (unsigned int i=0; i<nMag; i++) {
                 if (i) {
                     msg += ",";
                 }
                 msg += " '" + magcolVec[i] + "'";
             }
             msg += " ], magErrs=[";
             for (unsigned int i=0; i<nMagErr; i++) {
                 if (i) {
                     msg += ",";
                 }
                 msg += " '" + magerrcolVec[i] + "'";
             }
             msg += " ]";
             throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException, msg);
         }

	 if (idcol) {
             id = static_cast<int64_t*>(fitstable_read_column_inds(tag, idcol, i64, starinds, nstars));
             if (!id) {
                 throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                                   str(boost::format("Unable to read data for %s from %s") %
                                       idcol % ind->indexname));
             }
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

	 for (unsigned int j = 0; j != nMag; ++j) {
             if (magcolVec[j] == "") {
                 continue;
             }
             
             char const* magcol = magcolVec[j].c_str();

             if (j == 0) {		/* we need an extra copy for "flux" */
                 mag.push_back(read_column(tag, magcol, flt, starinds, nstars, ind->indexname));
             }
             
             mag.push_back(read_column(tag, magcol, flt, starinds, nstars, ind->indexname));
	 }
	 for (unsigned int j = 0; j != nMagErr; ++j) {
             if (magerrcolVec[j] == "") {
                 continue;
             }

             char const* magerrcol = magerrcolVec[j].c_str();
             if (j == 0) {		/* we need an extra copy for "flux" */
                 magerr.push_back(read_column(tag, magerrcol, flt, starinds, nstars, ind->indexname));
             }
             
             magerr.push_back(read_column(tag, magerrcol, flt, starinds, nstars, ind->indexname));
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
            if (!stargal) {
                throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                                  str(boost::format("Unable to read data for %s from %s") %
                                      stargalcol % ind->indexname));
            }
	    for (int j=0; j<nstars; j++) {
	       stargal[j] = (sg[j] > 0);
	    }
	    free(sg);
	 }
	 if (varcol) {
	    var = static_cast<bool*>(fitstable_read_column_inds(tag, varcol, boo, starinds, nstars));
            if (!var) {
                throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                                  str(boost::format("Unable to read data for %s from %s") %
                                      varcol % ind->indexname));
            }
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

	 if (id) {
             src->setId(id[i]);
         }

         if (nMag) {
             for (unsigned int j = 0; j != nMag + 1; ++j) { // +1 for "flux" column
                 // Dustin thinks converting to flux is 'LAME!';
                 // Jim thinks it's nice for consistency (and photocal wants fluxes
                 // as inputs, so we'll continue to go with that for now) even though
                 // we don't need to anymore.
                 // Dustin rebuts that photocal immediately converts those fluxes into mags,
                 // so :-P
                 double flux = pow(10.0, -0.4*mag[j][i]);
                 src->set(fluxKey[j], flux);
                 if (j < nMagErr) {
                     src->set(fluxErrKey[j], -0.4*magerr[j][i]*flux*std::log(10.0));
                 }
             }
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
      for (unsigned int j = 0; j < mag.size(); ++j) {
         free(mag[j]);
      }
      for (unsigned int j = 0; j < magerr.size(); ++j) {
          free(magerr[j]);
      }
      free(stargal);
      free(var);
      free(radecs);
      free(starinds);
   }
   return cat;
}
}}}}

