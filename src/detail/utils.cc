// -*- lsst-C++ -*-
// Astrometry.net include files:
extern "C" {
#include "astrometry/solver.h"
#include "astrometry/index.h"
#include "astrometry/starkd.h"
#include "astrometry/fitsioutils.h"
#include "astrometry/fitstable.h"

#undef ATTRIB_FORMAT
#undef FALSE
#undef TRUE
}

#include <set>
#include <cstdint>
#include "boost/format.hpp"

#include "lsst/meas/astrom/detail/utils.h"
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Calib.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/geom/Point.h"
#pragma clang diagnostic pop

namespace afwTable = lsst::afw::table;
namespace afwGeom  = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace astrom {
namespace detail {

static float *
read_column(
    fitstable_t* tag,
    char const *colName,
    tfits_type type,
    int *starinds,
    int nstars,
    char const *indexName)
{
    float *col = static_cast<float*>(fitstable_read_column_inds(tag, colName, type, starinds, nstars));
    if (!col) {
        throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundError,
                          str(boost::format("Unable to read data for %s from %s") % colName % indexName));
    }

    return col;
}

    
afwTable::SimpleCatalog
getCatalogImpl(std::vector<index_t*> inds,
    lsst::afw::coord::Coord const &ctrCoord,
    lsst::afw::geom::Angle const &radius,
    char const* idCol,
    std::vector<MagColInfo> const& magColInfoList,
    char const* isStarCol,
    char const* isVarCol,
    bool uniqueIds)
{
    /*
     If uniqueIds == true: return only reference sources with unique IDs;
     arbitrarily keep the first star found with each ID.
     */

    size_t const nMag = magColInfoList.size();   /* number of magnitude[error] columns */

    for (auto mc = magColInfoList.cbegin(); mc != magColInfoList.cend(); ++mc) {
        if (mc->filterName.empty()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterError,
                "Magnitude names cannot be empty strings.");
        }
        // We enforce this condition because we convert the mags to fluxes, and
        // we need a flux to compute a flux error!
        if (mc->magCol.empty()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterError,
                "Magnitude column names cannot be empty strings.");
        }
        //printf("mag col \"%s\", \"%s\", \"%s\"\n", mc->filterName.c_str(), mc->magCol.c_str(), mc->magErrCol.c_str());
    }
    
    auto icrsCoord = ctrCoord.toIcrs();
    double raDeg = icrsCoord.getLongitude().asDegrees();
    double decDeg = icrsCoord.getLatitude().asDegrees();
    double xyz[3];
    radecdeg2xyzarr(raDeg, decDeg, xyz);
    double r2 = deg2distsq(radius.asDegrees());

    afwTable::Schema schema = afwTable::SimpleTable::makeMinimalSchema(); // contains id and coord

    afw::table::PointKey<double>::addFields(schema, "centroid",
        "centroid on some exposure; invalid unless \"hasCentroid\" is true)", "pixels");
    auto hasCentroidKey = schema.addField<afwTable::Flag>("hasCentroid",
        "true if centroid field has been set");

    std::vector<afwTable::Key<double> > fluxKey;    // these are double for consistency with measured fluxes;
    std::vector<afwTable::Key<double> > fluxErrKey; // double may be unnecessary, but less surprising.
    fluxKey.reserve(nMag);
    fluxErrKey.reserve(nMag);

    for (auto mc = magColInfoList.cbegin(); mc != magColInfoList.cend(); ++mc) {
        // Add schema elements for each requested flux (and optionally flux error)
        // avoid the comment "flux flux"
        fluxKey.push_back(
            schema.addField<double>(
                mc->filterName + "_flux",
                mc->filterName + " flux"));
        if (mc->hasErr()) {
            fluxErrKey.push_back(
                schema.addField<double>(
                    mc->filterName + "_fluxSigma",
                    mc->filterName + " flux uncertainty (sigma)"));
        }
    }

    afwTable::Key<afwTable::Flag> resolvedKey;
    if (isStarCol) {
        resolvedKey = schema.addField<afwTable::Flag>(
            "resolved",
            "set if the reference object is resolved");
    }
    afwTable::Key<afwTable::Flag> variableKey;
    if (isVarCol) {
        variableKey = schema.addField<afwTable::Flag>(
            "variable",
            "set if the reference object is variable");
    }
    afwTable::Key<afwTable::Flag> photometricKey = schema.addField<afwTable::Flag>(
        "photometric",
        "set if the reference object can be used in photometric calibration");

    afwTable::SimpleCatalog cat;
    if (idCol) {
        // make catalog with no IdFactory, since IDs are external
        cat = afwTable::SimpleCatalog(afwTable::SimpleTable::make(schema, PTR(afwTable::IdFactory)()));
    } else {
        // let the catalog assign IDs
        cat = afwTable::SimpleCatalog(afwTable::SimpleTable::make(schema));
    }

    // for uniqueIds: keep track of the IDs we have already added to the result set.
    std::set<std::int64_t> uids;

    for (std::vector<index_t*>::iterator pind = inds.begin(); pind != inds.end(); ++pind) {
        index_t* ind = (*pind);
        // Find nearby stars
        double *radecs = NULL;
        int *starinds = NULL;
        int nstars = 0;
        startree_search_for(ind->starkd, xyz, r2, NULL, &radecs, &starinds, &nstars);
        //printf("found %i in \"%s\"\n", nstars, ind->indexname);
        if (nstars == 0) {
            continue;
        }

        std::vector<float*> mag;
        std::vector<float*> magErr;
        mag.reserve(nMag);
        magErr.reserve(nMag);
        std::int64_t* id = NULL;
        bool* stargal = NULL;
        bool* var = NULL;
        if (idCol || nMag || isStarCol || isVarCol) {
            fitstable_t* tag = startree_get_tagalong(ind->starkd);
            tfits_type flt = fitscolumn_float_type();
            tfits_type boo = fitscolumn_boolean_type();
            tfits_type i64 = fitscolumn_i64_type();

            if (!tag) {
                std::string msg = boost::str(boost::format(
                    "astrometry_net_data index file %s does not contain a tag-along table, "
                    "so can't retrieve extra columns.  idCol=%s, isStarCol=%s, isVarCol=%s") %
                    ind->indexname % idCol % isStarCol % isVarCol);
                 msg += ", mag columns=[";
                 for (unsigned int i=0; i<nMag; i++) {
                     if (i) {
                         msg += ",";
                     }
                     msg += " name='" + magColInfoList[i].filterName +
                         "', mag='" + magColInfoList[i].magCol +
                         "', magErr='" + magColInfoList[i].magErrCol + "'";
                 }
                 msg += " ].  You may need to edit the $ASTROMETRY_NET_DATA_DIR/andConfig.py file to set idColumn=None, etc.";
                 throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundError, msg);
            }

            if (idCol) {
                id = static_cast<int64_t*>(fitstable_read_column_inds(tag, idCol, i64, starinds, nstars));
                if (!id) {
                    throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundError,
                        str(boost::format("Unable to read data for %s from %s") % idCol % ind->indexname));
                }
            }
            if (id && uniqueIds) {
                // remove duplicate IDs.

                // FIXME -- this shouldn't be necessary once we get astrometry_net 0.40
                // multi-index functionality in place.

                if (uids.empty()) {
                    uids = std::set<std::int64_t>(id, id+nstars);
                } else {
                    int nkeep = 0;
                    for (int i=0; i<nstars; i++) {
                        //std::pair<std::set<std::int64_t>::iterator, bool> 
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

            for (auto mc = magColInfoList.cbegin(); mc != magColInfoList.cend(); ++mc) {
                char const* col = mc->magCol.c_str();
                mag.push_back(read_column(tag, col, flt, starinds, nstars, ind->indexname));
                if (mc->hasErr()) {
                    char const* col = mc->magErrCol.c_str();
                    magErr.push_back(read_column(tag, col, flt, starinds, nstars, ind->indexname));
                }
            }
            if (isStarCol) {
                /*  There is something weird going on with handling of bools; maybe "T" vs "F"?
                stargal = static_cast<bool*>(fitstable_read_column_inds(tag, isStarCol, boo, starinds, nstars));
                for (int j=0; j<nstars; j++) {
                printf("  sg %i = %i, vs %i\n", j, (int)sg[j], stargal[j] ? 1:0);
                }
                */
                uint8_t* sg = static_cast<uint8_t*>(fitstable_read_column_inds(
                    tag, isStarCol, fitscolumn_u8_type(), starinds, nstars));
                stargal = static_cast<bool*>(malloc(nstars));
                if (!stargal) {
                    throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundError,
                        str(boost::format("Unable to read data for %s from %s") % isStarCol % ind->indexname));
                }
                for (int j=0; j<nstars; j++) {
                    stargal[j] = (sg[j] > 0);
                }
                free(sg);
            }
            if (isVarCol) {
                var = static_cast<bool*>(fitstable_read_column_inds(tag, isVarCol, boo, starinds, nstars));
                if (!var) {
                    throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundError,
                        str(boost::format("Unable to read data for %s from %s") % isVarCol % ind->indexname));
                }
            }
        }

        for (int i=0; i<nstars; i++) {
            PTR(afwTable::SimpleRecord) src = cat.addNew();

            // Note that all coords in afwTable catalogs are ICRS; hopefully that's what the 
            // reference catalogs are (and that's what the code assumed before JFB modified it).
            src->setCoord(
                lsst::afw::coord::IcrsCoord(
                    radecs[i * 2 + 0] * afwGeom::degrees,
                    radecs[i * 2 + 1] * afwGeom::degrees
                )
            );

            if (id) {
                src->setId(id[i]);
            }

            src->set(hasCentroidKey, false);

            assert(fluxKey.size() == nMag);
            // only non-empty error columns are populated in these vectors.
            assert(fluxErrKey.size() == magErr.size());
            // index into non-empty error columns.
            size_t ej = 0;
            size_t j = 0;
            for (auto mc = magColInfoList.cbegin(); mc != magColInfoList.cend(); ++mc, ++j) {
                // There is some dispute about returning flux or magnitudes
                // flux is not conventional, but is convenient for several reasons
                // - it simplifies matching to sources for astrometric calibration
                // - flux errors are easier to combine than magnitude errors
                //printf("mag %s = %g\n", mc->filterName.c_str(), mag[j][i]);
                double flux = lsst::afw::image::fluxFromABMag(mag[j][i]);
                //printf("flux %g\n", flux);
                src->set(fluxKey[j], flux);
                if (mc->hasErr()) {
                    //printf("mag err = %g\n", magErr[ej][i]);
                    double fluxErr = lsst::afw::image::fluxErrFromABMagErr(magErr[ej][i], mag[j][i]);
                    //printf("flux err = %g\n", fluxErr);
                    src->set(fluxErrKey[ej], fluxErr);
                    ej++;
                }
            }
            assert(ej == fluxErrKey.size());

            bool photometric = true;
            if (stargal) {
                src->set(resolvedKey, !stargal[i]);
                photometric &= stargal[i];
            }
            if (var) {
                src->set(variableKey, var[i]);
                photometric &= (!var[i]);
            }
            src->set(photometricKey, photometric);
        }

        free(id);
        for (size_t j=0; j<mag.size(); ++j) {
            free(mag[j]);
        }
        for (size_t j=0; j<magErr.size(); ++j) {
            free(magErr[j]);
        }
        free(stargal);
        free(var);
        free(radecs);
        free(starinds);
    }
    return cat;
}

}}}}

