// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
#include "lsst/afw/geom.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
#include "lsst/daf/base/PropertyList.h"

#include <cstdio>
#include <iostream>

extern "C" {
#include "fitsioutils.h"
}

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace net {


using namespace std;
namespace afwGeom = lsst::afw::geom;
namespace afwImg = lsst::afw::image;
namespace afwCoord = lsst::afw::coord;
namespace afwGeom = lsst::afw::geom;
namespace Except = lsst::pex::exceptions;
namespace afwDet = lsst::afw::detection;
namespace pexLog = lsst::pex::logging;
namespace dafBase = lsst::daf::base;

int const USE_ALL_STARS_FOR_SOLUTION = -1;

static vector<double> getTagAlongFromIndex(index_t* index, string fieldName, int *ids, int numIds);

static afwCoord::Coord::Ptr xyztocoord(afwCoord::CoordSystem coordsys, const double* xyz) {
    afwCoord::Coord::Ptr radec = afwCoord::makeCoord(coordsys, afwGeom::Point3D(xyz[0], xyz[1], xyz[2]));
    return radec;
}

//
//Constructors, Destructors
//
    GlobalAstrometrySolution::GlobalAstrometrySolution(const std::string policyPath, pexLog::Log mylog):
    _mylog(mylog),
    _indexList(), 
    _solver(NULL), 
    _starxy(NULL), 
    _numBrightObjects(USE_ALL_STARS_FOR_SOLUTION),
    _isSolved(false) {

        // make qfits errors bubble up to astrometry.net to get logged.
        fits_use_error_system();
    
    _solver = solver_new();

    setDefaultValues();
    
    lsst::pex::policy::Policy pol(policyPath);
    _equinox = pol.getDouble("equinox");
    _raDecSys = pol.getString("raDecSys");
    string pkgDir = lsst::utils::eups::productDir("ASTROMETRY_NET_DATA");

    //Add meta information about every index listed in the policy file
    std::vector<std::string> indexArray = pol.getStringArray("indexFile");

    _mylog.log(pexLog::Log::DEBUG, "Loading astrometry_net_data metadata...");    
    for (unsigned int i = 0; i < indexArray.size(); ++i) {
        index_t *meta = _loadIndexMeta(pkgDir + "/" + indexArray[i]);
        bool duplicate = false;
        if (!meta)
            continue;
        // Check for duplicates.
        for (unsigned int j=0; j<_indexList.size(); j++) {
            index_t* other = _indexList[j];
            //These three values uniquely identify an index
            if (meta->indexid == other->indexid &&
                meta->healpix == other->healpix &&
                meta->hpnside == other->hpnside) {
                string msg = boost::str(boost::format("Index file \"%s\" is a duplicate (has same index id, healpix and healpix nside) as index file \"%s\"")
                                        % meta->indexname % other->indexname);
                _mylog.log(pexLog::Log::WARN, msg);
                duplicate = true;
                break;
            }
        }
        if (duplicate)
            continue;
        
        _indexList.push_back(meta);
    }
    _mylog.log(pexLog::Log::DEBUG, "astrometry_net_data metadata loaded.");
}
    
dafBase::PropertySet::Ptr GlobalAstrometrySolution::getMatchedIndexMetadata() {
    MatchObj* mo = solver_get_best_match(_solver);
    assert(mo);
    index_t* index = mo->index;
    assert(index);
    startree_t* starkd = index->starkd;
    assert(starkd);
    qfits_header* hdr = startree_header(starkd);
    char* val;

    dafBase::PropertyList::Ptr props(new dafBase::PropertyList());

    // reference catalog name
    val = fits_get_dupstring(hdr, "REFCAT");
    props->add("REFCAT", std::string(val ? val : "none"),
               "Reference catalog name");
    free(val);

    // reference catalog md5sum
    val = fits_get_dupstring(hdr, "REFCAMD5");
    props->add("REFCAMD5", std::string(val ? val : "none"),
               "Reference catalog MD5 checksum");
    free(val);

    // Astrometry.net index
    props->add("ANINDID", index->indexid, "Astrometry.net index id");
    props->add("ANINDHP", index->healpix, "Astrometry.net index HEALPix");
    props->add("ANINDNM", std::string(index->indexname),
               "Astrometry.net index name");

    return props;
}

void GlobalAstrometrySolution::loadIndices() {
    for (unsigned int i=0; i<_indexList.size(); i++)
        index_reload(_indexList[i]);
}

std::vector<const index_t*> GlobalAstrometrySolution::getIndexList() {
    std::vector<const index_t*> rtn;
    for (unsigned int i=0; i<_indexList.size(); i++)
        rtn.push_back(_indexList[i]);
    return rtn;
}

std::vector<int> GlobalAstrometrySolution::getIndexIdList() {
    std::vector<int> rtn;
    for (unsigned int i=0; i<_indexList.size(); i++)
        rtn.push_back(_indexList[i]->indexid);
    return rtn;
}

index_t *GlobalAstrometrySolution::_loadIndexMeta(std::string filename){
    return index_load(filename.c_str(), INDEX_ONLY_LOAD_METADATA, NULL);
}

GlobalAstrometrySolution::~GlobalAstrometrySolution() {
    for (unsigned int i=0; i<_indexList.size(); i++) {
        index_free(_indexList[i]);
    }
    _indexList.clear();
    if ( _starxy != NULL) {
        starxy_free(_starxy);
        _starxy = NULL;
    }
    if (_solver != NULL){
        solver_free(_solver);
        _solver = NULL;
    }
}

void GlobalAstrometrySolution::setDefaultValues() {
    // Among other things, this sets the parity, positional uncertainty
    // (_solver->verify_pix) and matching accuracy (_solver->codetol)
    solver_set_default_values(_solver);
    // Set image scale limits (in arcseconds per pixel) to non-zero and
    // non-infinity.
    setMinimumImageScale(1e-6 * afwGeom::arcseconds);
    setMaximumImageScale(360 * afwGeom::degrees);
    // How good must a match be to be considered good enough?  Log-odds.
    setMatchThreshold(log(1e12));
    // Handedness of the CCD coordinate system
    setParity(UNKNOWN_PARITY);
    // Reset counters and record of best match found so far.
    solver_cleanup_field(_solver);
}

    
//
// Setup functions
//
void GlobalAstrometrySolution::setImageSize(int W, int H) {
	 solver_set_field_bounds(_solver, 0, W, 0, H);
}

///Set the image to be solved. The image is abstracted as a list of positions in pixel space
void GlobalAstrometrySolution::setStarlist(afwDet::SourceSet vec ///<List of Sources
                                          ) {
    if (vec.empty()) {
        throw(LSST_EXCEPT(pexExcept::LengthErrorException, "Source list contains no objects"));
    }

    // Step 1. Copy every valid element of the input vector into a tempory starlist structure
    // Valid means all of x, y and psfFlux are positive and finite.
    int const size = vec.size();
    starxy_t *tmpStarxy = starxy_new(size, true, false);   

    int i = 0;
    for (afwDet::SourceSet::iterator ptr = vec.begin(); ptr != vec.end(); ++ptr) {
        double const x    = (*ptr)->getXAstrom();
        double const y    = (*ptr)->getYAstrom();
        double const flux = (*ptr)->getPsfFlux();

        // Only include objects where positions and fluxes are positive finite values
        if( isfinite(x)    && (x>=0) &&
            isfinite(y)    && (y>=0) &&
            isfinite(flux) && (flux>0)
          ) {
            starxy_set(tmpStarxy, i, x, y);
            starxy_set_flux(tmpStarxy, i, flux);
            ++i;
        }
    }

    int nwarn = 15;
    if (i < nwarn) {
        string msg = boost::str(boost::format("Source list only has %i valid objects; probably need %i or more\n") % i % nwarn);
        msg += "Valid objects have positive, finite values for x, y and psfFlux";
        _mylog.log(pexLog::Log::WARN, msg);
    }

    // Step 2. Copy these elements into a new starxy structure of the correct
    // size. If we don't do this, starxy will report its size incorrectly, and
    // we'll get confused later on.
    if (_starxy)
        starxy_free(_starxy);
    
    if (i == size) {
        // Every part of tmpStarxy is valid
        _starxy = tmpStarxy;
    } else {
        _starxy = starxy_new(i, true, false);
        for(int j=0; j<i; ++j) {
            starxy_set_x( _starxy, j, starxy_get_x( tmpStarxy, j));
            starxy_set_y( _starxy, j, starxy_get_y( tmpStarxy, j));
            starxy_set_flux( _starxy, j, starxy_get_flux( tmpStarxy, j));
        }
        starxy_free(tmpStarxy);
    }

    // Sort the array
    starxy_sort_by_flux(_starxy);
    _solverSetField();
}


///\brief Only find a solution using the brightest N objects
///Reducing the number of sources in the solution list can reduce the time taken to solve the image
///The truncated list is used by solve() to get the linear wcs, but all input sources are used when
///calculating the distortion terms of the SIP matrix.    
void GlobalAstrometrySolution::setNumBrightObjects(int N) {
    if (N <= 0)
        throw(LSST_EXCEPT(pexExcept::RangeErrorException, "setNumBrightObjects: must be positive"));

    _numBrightObjects = N;
    if (_starxy)
        _solverSetField();
}


void GlobalAstrometrySolution::_solverSetField() {
    if ( ! _starxy) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    int starxySize = starxy_n(_starxy);
    if (starxySize == 0) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist has zero elements"));
    }

    int N = _numBrightObjects;
    //The default value, -1, indicates that all objects should be used
    if (N == USE_ALL_STARS_FOR_SOLUTION) {
        N = starxySize;
    }

    starxy_t *shortlist = starxy_subset(_starxy, N);
    assert(shortlist);

    // Set the pointer in the solver to the new, smaller field
    starxy_free(solver_get_field(_solver));
    solver_free_field(_solver);
    solver_set_field(_solver, shortlist);
    solver_reset_field_size(_solver);

    // Find field boundaries and precompute kdtree
    solver_preprocess_field(_solver);
}   

// Set the plate scale of the image in arcsec per pixel
void GlobalAstrometrySolution::setImageScaleArcsecPerPixel(lsst::afw::geom::Angle imgScale ///< Plate scale of image
                                                          ) {
    //Note that the solver will fail if min==max, so we make them different by a small amount.    
    setMinimumImageScale(0.99 * imgScale);
    setMaximumImageScale(1.01 * imgScale);
}


///Set the verbosity level for astrometry.net. The higher the level the more information is returned.
///1 and 2 are typically good values to use. 4 will print so much to the screen that it slows execution
void GlobalAstrometrySolution::setLogLevel(int level) {
    if (level < 0 || level > 4) {
        throw(LSST_EXCEPT(pexExcept::DomainErrorException, "Logging level must be between 0 and 4"));
    }
    log_init((enum log_level) level);
    fits_use_error_system();
}

///How good does a match need to be to be accepted. Typical value is log(1e12) approximately 27
void GlobalAstrometrySolution::setMatchThreshold(double threshold) {
    solver_set_record_logodds(_solver, threshold);
}

///You can double the speed of a match if you know the parity, i.e whether the image is flipped or not.
///North up and East right (or some rotation thereof) is parity==NORMAL_PARITY, the opposite is
///parity==FLIPPED_PARITY. The default is UNKNOWN_PARITY
void GlobalAstrometrySolution::setParity(int parity){
    //Insist on legal values.
    if (solver_set_parity(_solver, parity)) {
        throw LSST_EXCEPT(pexExcept::DomainErrorException, "Illegal parity setting");
    }
}

//
// Solve functions
//

// Solve using a wcs object as an initial guess
bool GlobalAstrometrySolution::solve(const lsst::afw::image::Wcs::Ptr wcsPtr, 
                                     double imageScaleUncertaintyPercent) {
    double unc = imageScaleUncertaintyPercent/100.0;
    double xc, yc;
    solver_get_field_center(_solver, &xc, &yc);

    // Get the central ra/dec and pixel scale
    afwCoord::Coord::ConstPtr raDec = wcsPtr->pixelToSky(xc, yc);
    afwGeom::Angle ra  = raDec->getLongitude();
    afwGeom::Angle dec = raDec->getLatitude();
    _mylog.log(pexLog::Log::DEBUG,
               boost::format("Solving using initial guess at position of\n %.7f %.7f deg\n") % ra.asDegrees() % dec.asDegrees());

    afwGeom::Angle pixelScale = wcsPtr->pixelScale();
    afwGeom::Angle lwr = pixelScale*(1 - unc);
    afwGeom::Angle upr = pixelScale*(1 + unc);
    setMinimumImageScale(lwr);
    setMaximumImageScale(upr);
    _mylog.log(pexLog::Log::DEBUG, boost::format("Exposure's WCS scale: %g arcsec/pix; setting scale range %.3f - %.3f arcsec/pixel\n") %
               pixelScale.asArcseconds() % lwr.asArcseconds() % upr.asArcseconds());

    if ( wcsPtr->isFlipped()) {
        setParity(FLIPPED_PARITY);
        _mylog.log(pexLog::Log::DEBUG, "Setting Flipped parity");        
    } else {
        setParity(NORMAL_PARITY);
        _mylog.log(pexLog::Log::DEBUG, "Setting Normal parity");        
    }
        
    return solve(ra, dec);
}

///Find a solution with an initial guess at the position    
bool GlobalAstrometrySolution::solve(afwCoord::Coord::ConstPtr raDec) {
    return solve(raDec->getLongitude(), raDec->getLatitude());
}

///Find a solution with an initial guess at the position.
bool GlobalAstrometrySolution::solve(afwGeom::Angle ra,   ///<Right ascension
                                     afwGeom::Angle dec   ///< Declination
                                          )  {    
    string msg;

    // Tell the solver to only consider matches within the image size of the supposed RA,Dec.
    // The magic number 2.0 out front says to accept matches within 2 radii of the given *center* position.
    afwGeom::Angle maxRadius = (2.0 * solver_get_max_radius_arcsec(_solver)) * afwGeom::arcseconds;
    msg = boost::str(boost::format("Setting RA,Dec = (%g, %g), radius = %g deg") % ra.asDegrees() % dec.asDegrees() % maxRadius.asDegrees());
    _mylog.log(pexLog::Log::DEBUG, msg);
    solver_set_radec(_solver, ra.asDegrees(), dec.asDegrees(), maxRadius.asDegrees());

    _mylog.log(pexLog::Log::DEBUG, "Doing solve step");
    _callSolver(ra, dec);

    if (_isSolved){
        const char* indexname = solver_get_best_match_index_name(_solver);
        msg = boost::str(boost::format("Solved index is %s") % indexname);
    } else {
        msg = boost::str(boost::format("Failed to verify position (%.7f %.7f)\n") % ra.asDegrees() % dec.asDegrees());
    }
    _mylog.log(pexLog::Log::DEBUG, msg);            

    return _isSolved;
}

// Find a solution blindly (with no initial guess)
bool GlobalAstrometrySolution::solve()  {    

    // Don't use any hints about the RA,Dec position.
    solver_clear_radec(_solver);

    _callSolver(afwGeom::NullAngle, afwGeom::NullAngle);

    if (_isSolved){
        const char* indexname = solver_get_best_match_index_name(_solver);
        string msg = boost::str(boost::format("Astrometric solution found with index %s") % indexname);
        _mylog.log(pexLog::Log::DEBUG, msg);
    } else {
        _mylog.log(pexLog::Log::DEBUG, "Failed to find astrometric solution");
    }

    return _isSolved;
}

// Check that all the setup was done correctly, then set the solver to work By
// default, _callSolver does a blind solve, unless the optional ra and dec are
// given values. Sets _isSolved to true if a solution is found
bool GlobalAstrometrySolution::_callSolver(afwGeom::Angle ra, afwGeom::Angle dec) {
    //Throw exceptions if setup is incorrect
    if ( !_starxy) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    if (_indexList.size() == 0) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No index files loaded yet"));
    }
    
    if (_isSolved) {
        string msg = "Solver indicated that a match has already been found. Do you need to reset?";
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
    }

    afwGeom::Angle lower = solver_get_pixscale_low (_solver) * afwGeom::arcseconds;
    afwGeom::Angle upper = solver_get_pixscale_high(_solver) * afwGeom::arcseconds;

    if (lower >= upper) {
        string msg = boost::str(boost::format("Minimum image scale (%g) must be strictly less than max scale (%g)") % lower.asArcseconds() % upper.asArcseconds());
        throw(LSST_EXCEPT(pexExcept::DomainErrorException, msg));
    }

    //Calculate the best guess at image size
    double xSizePixels = solver_field_width(_solver);
    double ySizePixels = solver_field_height(_solver);
    double minSizePixels = min(xSizePixels, ySizePixels);
    double maxSizePixels = max(xSizePixels, ySizePixels);
    assert(maxSizePixels >= minSizePixels && minSizePixels > 0);

    //Set the range of sizes of quads to examine, in pixels.
    solver_set_quad_size_fraction(_solver, 0.1, 1.0);

    //Output some useful debugging info
    _mylog.format(pexLog::Log::DEBUG, "Image size %.0f x %.0f pixels", xSizePixels, ySizePixels);
    _mylog.format(pexLog::Log::DEBUG, "Searching plate scale range %.3f -- %.3f arcsec/pixel",
                  lower.asArcseconds(), upper.asArcseconds());
    _mylog.format(pexLog::Log::DEBUG, "--> Image size %.3f x %.3f to %.3f x %.3f arcmin",
                  lower.asArcminutes() * xSizePixels, upper.asArcminutes() * ySizePixels,
                  lower.asArcminutes() * xSizePixels, upper.asArcminutes() * ySizePixels);

    double qlo, qhi;
    solver_get_quad_size_range_arcsec(_solver, &qlo, &qhi);
    _mylog.format(pexLog::Log::DEBUG, "Using indices with quads in the range %.2f to %.2f arcmin\n",
                  arcsec2arcmin(qlo), arcsec2arcmin(qhi));

    _mylog.log(pexLog::Log::DEBUG, "Setting indices");
    _addSuitableIndicesToSolver(qlo * afwGeom::arcseconds, qhi * afwGeom::arcseconds, ra, dec);
    _mylog.log(pexLog::Log::DEBUG, "Doing solve step");
    solver_run(_solver);
    
    if (solver_did_solve(_solver)) {
        _isSolved = true;
        MatchObj* match = solver_get_best_match(_solver);
        _mylog.format(pexLog::Log::DEBUG, "Solved: %i matches, %i conflicts, %i unmatched",
                      (int)match->nmatch, (int)match->nconflict, (int)match->ndistractor);
        _mylog.log(pexLog::Log::DEBUG, "Calling tweak2() to tune up match...");
        _mylog.format(pexLog::Log::DEBUG, "Starting log-odds: %g", match->logodds);
        // Use "tweak2" to tune up this match, resulting in a better WCS and more catalog matches.
        // magic 1: only go to linear order (no SIP distortions).

	// HACK -- astrometry_net 0.30
        solver_tweak2(_solver, match, 1);
	// -- astrometry_net > 0.30:
        //solver_tweak2(_solver, match, 1, NULL);
        _mylog.format(pexLog::Log::DEBUG, "After tweak2: %i matches, %i conflicts, %i unmatched",
                      (int)match->nmatch, (int)match->nconflict, (int)match->ndistractor);
        _mylog.format(pexLog::Log::DEBUG, "After tweak2: log-odds: %g", match->logodds);
    } else {
        _isSolved = false;
    }
    
    return _isSolved;
}



/// \brief Find indices that may contain a the correct solution, and add them to the solver.
/// 
/// Find indices that cover a suitable range of plate scales and optionally a suitable position
/// on the sky.
/// If these indices have been previously loaded from disk, add them to the solver, otherwise
/// load them from disk, then add them to the solver.
/// Because these files are typically large, caching them in memory can save a lot of 
/// initialisation time.
/// Because each index can take a long time to search, only using suitable ones speeds
/// the matching process.
/// 
/// If ra and dec are not supplied, the default values (see the header file) indicate to the 
/// function that position should not be used to determine the suitablity of an index
/// 
/// \param ra Optional right ascension of inital guess at solution position
/// \param dec Optional declination of inital guess at solution position
/// 
/// \return Number of indices loaded
int GlobalAstrometrySolution::_addSuitableIndicesToSolver(afwGeom::Angle quadSizeLow,
                                                          afwGeom::Angle quadSizeHigh,
                                                          afwGeom::Angle ra, afwGeom::Angle dec) {
    bool hasAtLeastOneIndexOfSuitableScale = false;
    int nMeta = _indexList.size();
    int nSuitable = 0;
    bool blind = (ra == afwGeom::NullAngle) || (dec == afwGeom::NullAngle);

    for (int i = 0; i<nMeta; ++i){
        index_t* index = _indexList[i];

        if (!index_overlaps_scale_range(index, quadSizeLow.asArcseconds(), quadSizeHigh.asArcseconds()))
            continue;

        hasAtLeastOneIndexOfSuitableScale = true;

        //Is this either a blind solve, or does the index cover a suitable
        //patch of sky
        if (!(blind || index_is_within_range(index, ra.asDegrees(), dec.asDegrees(),
                                             quadSizeHigh.asDegrees())))
            continue;

        // Found a good one!

        // Load the data (not just metadata), if it hasn't been already...
        index_reload(index);

        string msg = boost::str(boost::format("Adding index %s") % index->indexname);
    	_mylog.log(pexLog::Log::DEBUG, msg);

        solver_add_index(_solver, index);
        nSuitable++;
    }
    
    if(nSuitable == 0) {
        string msg = "No suitable indices found for given input parameters:";
        
        if(hasAtLeastOneIndexOfSuitableScale) {
            msg += "Probably the ra/dec range isn't covered";
        } else {
            msg += "No indices of a suitable scale were found";
        }
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
    }
    
    return nSuitable;
}


//
//Return the solution
//              

///After solving, return a linear Wcs (i.e without distortion terms)
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getWcs()  {

    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    MatchObj* match = solver_get_best_match(_solver);

    afwGeom::Point2D crpix = afwGeom::Point2D(match->wcstan.crpix[0], match->wcstan.crpix[1]);
    afwGeom::Point2D crval = afwGeom::Point2D(match->wcstan.crval[0], match->wcstan.crval[1]);
    
    int naxis = 2;
    Eigen::Matrix2d CD;
    for (int i = 0; i<naxis; ++i) {
        for (int j = 0; j<naxis; ++j) {
            CD(i, j) = match->wcstan.cd[i][j];
        }
    }
    std::string const ctype1 = "RA---TAN";
    std::string const ctype2 = "DEC--TAN";
    return lsst::afw::image::Wcs::Ptr(new lsst::afw::image::Wcs(crval, crpix, CD,
                                                                ctype1, ctype2, _equinox, _raDecSys)); 
}

///After solving, return a full Wcs including SIP distortion matrics
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getDistortedWcs(int order)  {
    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }
    if (!_starxy) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist isn't set"));
    }

    //Generate an array of radec of positions in the field
    MatchObj* mo = solver_get_best_match(_solver);

    //Call the tweaking algorthim to generate the distortion coeffecients
    
    //jitter is a measure of how much we can expect the xy of stars to scatter from the expected
    //radec due to noise in our measurments.
    // don't use Angle here
    double jitterArcsec = tan_pixel_scale(&mo->wcstan) * solver_get_field_jitter(_solver);
    jitterArcsec = hypot(jitterArcsec, mo->index_jitter);
    int inverseOrder = order;
    int iterations = 5;        //blind.c:628 uses 5
    bool isWeighted = true;
    int skipShift = true;

    sip_t *sip = tweak_just_do_it(&mo->wcstan, _starxy, mo->refxyz,
                                  NULL, NULL, NULL, mo->nindex,
                                  jitterArcsec, order, inverseOrder,
                                  iterations, isWeighted, skipShift);

    //Check that tweaking worked.
    if (sip == NULL) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Tweaking failed"));
    }

    afwGeom::Point2D crpix = afwGeom::Point2D(sip->wcstan.crpix[0], sip->wcstan.crpix[1]);
    afwGeom::Point2D crval = afwGeom::Point2D(sip->wcstan.crval[0], sip->wcstan.crval[1]);

    //Linear conversion matrix
    int naxis = 2;
    Eigen::Matrix2d CD;
    for (int i = 0; i<naxis; ++i) {
        for (int j = 0; j<naxis; ++j) {
            CD(i, j) = sip->wcstan.cd[i][j];
        }
    }

    //Forward distortion terms. In the SIP notation, these matrices are referred to
    //as A and B. I can find no documentation that insists that these matrices be
    //the same size, so I assume they aren't.
    int aSize = sip->a_order + 1;
    Eigen::MatrixXd sipA(aSize, aSize);
    for (int i = 0; i<aSize; ++i){
        for (int j = 0; j<aSize; ++j){
            sipA(i, j) = sip->a[i][j];
        }
    }

    //Repeat for B
    int bSize = sip->b_order + 1;
    Eigen::MatrixXd sipB(bSize, bSize);
    for (int i = 0; i<bSize; ++i){
        for (int j = 0; j<bSize; ++j){
            sipB(i, j) = sip->b[i][j];
        }
    }

    //Repeat for Ap, for the reverse transform
    int apSize = sip->ap_order + 1;
    Eigen::MatrixXd sipAp(apSize, apSize);
    for (int i = 0; i<apSize; ++i){
        for (int j = 0; j<apSize; ++j){
            sipAp(i, j) = sip->ap[i][j];
        }
    }

    //And finally, Bp, also part of the reverse transform
    int bpSize = sip->bp_order + 1;
    Eigen::MatrixXd sipBp(bpSize, bpSize);
    for (int i = 0; i<bpSize; ++i){
        for (int j = 0; j<bpSize; ++j){
            sipBp(i, j) = sip->bp[i][j];
        }
    }    
        
    lsst::afw::image::Wcs::Ptr
        wcsPtr(new lsst::afw::image::TanWcs(crval, crpix, CD, sipA, sipB, sipAp, sipBp, _equinox, _raDecSys));
    
    sip_free(sip);
    return wcsPtr;
}    

static vector<boost::int64_t> getInt64TagAlongFromIndex(const index_t* index, string fieldName, const int *ids, int numIds) {
    int64_t *tagAlong=NULL;
    string msg;
    vector<boost::int64_t> out(0);
    
    if (fieldName != "") {
        // Grab tag-along data here. If it's not there, throw an exception
        if (!startree_has_tagalong(index->starkd) ) {
            msg = boost::str(boost::format("Index file \"%s\" has no tag-along data") % index->indexname);
            throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
        }
        
	// in later Astrometry.net versions:
#if 0
        tagAlong = startree_get_data_column_int64(index->starkd, fieldName.c_str(), ids, numIds);
#else
	fitstable_t* table = startree_get_tagalong(index->starkd);
	if (!table) {
	  throw LSST_EXCEPT(pexExcept::RuntimeErrorException, 
			    boost::str(boost::format("Index file \"%s\": failed to startree_get_tagalong") % index->indexname));
	}
	tagAlong = (int64_t*)fitstable_read_column_inds(table, fieldName.c_str(), fitscolumn_i64_type(), ids, numIds);
#endif

        if (!tagAlong) {
	  /*
            msg = boost::str(boost::format("No meta data (tag-along) column called \"%s\" found in index %s") %
	    fieldName % index->indexname);
            throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
	  */
	  return out;
        }

        out = vector<boost::int64_t>(&tagAlong[0], &tagAlong[numIds]);
        // the vector constructor makes a copy of the data.
        free(tagAlong);
        return out;
    }
    return out;
}

static vector<boost::int64_t> getIds(const string idName, const index_t* index, const int* starinds, int nstars) {
    vector<boost::int64_t> ids;
    if (idName != "") {
      ids = getInt64TagAlongFromIndex(index, idName, starinds, nstars);
    }
    /*
      if (!ids.size()) {
      // Define the ID to be the index -- better than nothing!
      for (int j=0; j<nstars; j++) {
      ids.push_back(starinds[j]);
      }
      }
    */
    return ids;
}

lsst::afw::coord::CoordSystem GlobalAstrometrySolution::_getCoordSys() {
    return afwCoord::makeCoordEnum(_raDecSys);
}


///\brief Return a list of the sources that match and the catalogue objects they matched to.
///
///For each object in the catalogue that matches an input source return a SourceMatch object
///containing a) the catalogue object, b) the source object, c) the distance between them in pixels
///For the two objects, X,YAstrom and Ra/Dec are set.
vector<afwDet::SourceMatch> GlobalAstrometrySolution::getMatchedSources(string filterName, string idName){
    if (!_isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }
    
    string msg = "";
    MatchObj* match = solver_get_best_match(_solver);
    assert(match->nfield > 0);
    
    vector<afwDet::SourceMatch> sourceMatchSet;
    afwImg::Wcs::Ptr wcsPtr = this->getWcs();

    // FIXME -- wait, doesn't match->theta possibly have non-matched (negative) entries?!
    // Load magnitude information from catalogue
    vector<double> refMag = getTagAlongFromIndex(match->index, filterName, match->theta, match->nfield);
    vector<boost::int64_t> refId = getIds(idName, match->index, match->theta, match->nfield);

    const starxy_t* fieldxy = solver_get_field(_solver);

    // from metadata.paf via _raDecSys
    afwCoord::CoordSystem coordsys = _getCoordSys();
        
    for (int i=0; i<match->nfield; i++) {
        // "theta" is the mapping from image (aka field) stars to index (aka reference) stars.
        // negative means no match.
        if (match->theta[i] < 0)
            continue;
        
        //Matching input sources    
        afwDet::Source::Ptr src(new afwDet::Source());
        double x = starxy_getx(fieldxy, i);
        double y = starxy_gety(fieldxy, i);
        double flux = starxy_get_flux(fieldxy, i);
        
        src->setXAstrom(x);
        src->setYAstrom(y);
        src->setPsfFlux(flux);
        
        //Weight positions by confidence in the fact that they match
        // FIXME -- this is insane.  What units are these in?
        double confidence = verify_logodds_to_weight(match->matchodds[i]); //Can be == 0
        src->setXAstromErr( 1/(confidence + DBL_EPSILON));
        src->setYAstromErr( 1/(confidence + DBL_EPSILON));
        
        src->setRaDecAstromFromXy(wcsPtr);
        
        // Matching catalogue objects
        afwDet::Source::Ptr ref(new afwDet::Source());
        double* xyz = match->refxyz + match->theta[i] * 3;
        afwCoord::Coord::Ptr radec = xyztocoord(coordsys, xyz);
        ref->setAllRaDecFields(radec);
        ref->setAllXyFromRaDec(wcsPtr);

        // FIXME -- should set the other error fields -- but not to these values!
        ref->setXAstromErr( 1/(confidence + DBL_EPSILON));
        ref->setYAstromErr( 1/(confidence + DBL_EPSILON));

        if (refMag.size()) {
            ref->setPsfFlux( pow(10.0, -refMag[i]/2.5) );
        }
        if (refId.size()) {
	  ref->setSourceId(refId[i]);
        }

        double dx = ref->getXAstrom() - src->getXAstrom();
        double dy = ref->getYAstrom() - src->getYAstrom();
        double dist = hypot(dx, dy);
        sourceMatchSet.push_back(afwDet::SourceMatch(ref, src, dist));
    }

    return sourceMatchSet;
}


///Astrometry.net catalogues store additional data about objects in addition to their position (usually 
///magnitudes. This function returns the names of the strings used to describe those fields
///
///\note This function makes a potentially untrue assumption. It assumes that the meta data fields
///are uniform in all indices read from disk, so it only checks for metadata in one index file
vector<string> GlobalAstrometrySolution::getCatalogueMetadataFields() {

    vector<string> output;
    
    //Reload the index if necssary
    index_reload(_indexList[0]);
    if (! startree_has_tagalong(_indexList[0]->starkd) ) {
        _mylog.log(pexLog::Log::DEBUG, "No tag-along data found for index");
        return output;
    }
    
    fitstable_t *table = startree_get_tagalong(_indexList[0]->starkd);
    assert(table != NULL);
    sl *nameList = fitstable_get_fits_column_names(table, NULL);
    int numNames = sl_size(nameList);
    for(int i=0; i< numNames; ++i) {
        string name(sl_pop(nameList));
        output.push_back(name);
    }
    sl_free2(nameList);
        
    return output;
}

///Returns a sourceSet of objects that are nearby in an raDec sense to the best match solution
ReferenceSources
GlobalAstrometrySolution::getCatalogueForSolvedField(string filterName, string idName, double margin) {
    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }
    MatchObj* match = solver_get_best_match(_solver);
    double* xyz;
    int* starinds;
    int nstars;
    double scale;
    double r2;
    int i, W, H;
    int outi;
    
    ReferenceSources refs;
    afwDet::SourceSet ss;

    // arcsec/pix
    scale = tan_pixel_scale(&(match->wcstan));
    // add margin
    r2 = deg2distsq(match->radius_deg + arcsec2deg(scale * margin));

    refs.indexid = match->index->indexid;
    
    startree_search_for(match->index->starkd, match->center, r2, &xyz, NULL, &starinds, &nstars);
    if (nstars == 0)
        return refs;

    afwCoord::CoordSystem coordsys = _getCoordSys();

    W = (int)(match->wcstan.imagew);
    H = (int)(match->wcstan.imageh);
    outi = 0;
    for (i=0; i<nstars; i++) {
        double* xyzi;
        double px, py;
        xyzi = xyz + 3*i;
        if (!tan_xyzarr2pixelxy(&(match->wcstan), xyzi, &px, &py))
            continue;
        // in bounds (+ margin) ?
        if (px < -margin || px > W+margin || py < -margin || py > H+margin)
            continue;

        afwDet::Source::Ptr src(new afwDet::Source());
        src->setRaDec(xyztocoord(coordsys, xyzi));
        ss.push_back(src);
        starinds[outi] = starinds[i];
        outi++;
    }

    vector<double> mag = getTagAlongFromIndex(match->index, filterName, starinds, outi);
    if (mag.size()) {
        for (unsigned int i=0; i<ss.size(); i++) {
            // It seems crazy to convert back to flux...
            ss[i]->setPsfFlux(pow(10.0, -mag[i]/2.5));
        }
    }

    vector<boost::int64_t> ids = getIds(idName, match->index, starinds, outi);
    if (ids.size()) {
        for (unsigned int i=0; i<ss.size(); i++) {
            ss[i]->setSourceId(ids[i]);
        }
    }

    vector<int> inds(starinds, starinds + ss.size());

    free(starinds);
    free(xyz);

    refs.refsources = ss;
    refs.inds = inds;

    return refs;
}
    

///Returns a sourceSet of objects that are nearby in an RA,Dec sense to the best match solution
afwDet::SourceSet 
GlobalAstrometrySolution::getCatalogue(afwGeom::Angle radius, string filterName, string idName) {
    if (!_isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    // Don't free this!
    MatchObj* match = solver_get_best_match(_solver);
    double *center = match->center;
    double ra_deg, dec_deg;
    xyzarr2radecdeg(center, &ra_deg, &dec_deg);
    afwGeom::Angle ra  = ra_deg  * afwGeom::degrees;
    afwGeom::Angle dec = dec_deg * afwGeom::degrees;
    ReferenceSources refs = getCatalogue(ra, dec, radius, filterName, idName);
    return refs.refsources;
}

index_t* GlobalAstrometrySolution::_getIndex(int indexId) {
    for (unsigned int i=0; i<_indexList.size(); i++) {
        index_t* ind = _indexList[i];
        if (ind->indexid != indexId)
            continue;
        return ind;
    }
    return NULL;
}

template <typename T>
std::vector<T> GlobalAstrometrySolution::_getTagAlongData(int indexId, std::string columnName,
                                                          tfits_type ctype, std::vector<int> inds) {

    index_t* index = _getIndex(indexId);
    if (!index)
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          boost::str(boost::format("Astrometry.net index with ID %i was not found") % indexId)));

    fitstable_t* tag = startree_get_tagalong(index->starkd);
    if (!tag)
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          boost::str(boost::format("Astrometry.net index with ID %i: no tag-along table was found") % indexId)));

    int* cinds = new int[inds.size()];
    for (size_t i=0; i<inds.size(); i++)
        cinds[i] = inds[i];
    
    int arraysize = -1;
    void* vdata = fitstable_read_column_array_inds(tag, columnName.c_str(), ctype, cinds, inds.size(), &arraysize);
    delete[] cinds;
    std::vector<T> vals;

    if (ctype != fitscolumn_boolean_type()) {
        T* x = static_cast<T*>(vdata);
        vals = std::vector<T>(x + 0, x + inds.size() * arraysize);
    } else {
        // Workaround a bug in Astrometry.net 0.30: for boolean type, the "vdata"
        // values are 'T' and 'F' -- can't just cast them to 'bool'!!
        unsigned char* cdata = static_cast<unsigned char*>(vdata);
        T* x = new T[inds.size()];
        for (int i=0, N=inds.size(); i<N*arraysize; i++) {
            if (cdata[i] == 'T')
                x[i] = true;
            else if (cdata[i] == 'F')
                x[i] = false;
            else
                throw(LSST_EXCEPT(pexExcept::RuntimeErrorException,
                                  boost::str(boost::format("Error retrieving boolean column \"%s\" from FITS file: entry %i has value 0x%x")
                                             % columnName.c_str() % i % (int)cdata[i])));
        }
        vals = std::vector<T>(x + 0, x + inds.size() * arraysize);
        delete[] x;
    }
    free(vdata);
    return vals;
}

std::vector<double> GlobalAstrometrySolution::getTagAlongDouble(int indexId, std::string columnName,
                                                                std::vector<int> inds) {
    std::vector<double> vals = _getTagAlongData<double>(indexId, columnName, fitscolumn_double_type(), inds);
    return vals;
}

std::vector<int> GlobalAstrometrySolution::getTagAlongInt(int indexId, std::string columnName,
                                                          std::vector<int> inds) {
    std::vector<int> vals = _getTagAlongData<int>(indexId, columnName, fitscolumn_int_type(), inds);
    return vals;
}

std::vector<boost::int64_t> GlobalAstrometrySolution::getTagAlongInt64(int indexId, std::string columnName,
                                                                       std::vector<int> inds) {
    std::vector<boost::int64_t> vals = _getTagAlongData<boost::int64_t>(indexId, columnName, fitscolumn_i64_type(), inds);
    return vals;
}

std::vector<bool> GlobalAstrometrySolution::getTagAlongBool(int indexId, std::string columnName,
                                  std::vector<int> inds) {
    std::vector<bool> vals = _getTagAlongData<bool>(indexId, columnName, fitscolumn_boolean_type(), inds);
    return vals;
}

std::vector<TagAlongColumn> GlobalAstrometrySolution::getTagAlongColumns(int indexId) {
    index_t* index;
    if (indexId == -1)
        index = _indexList[0];
    else {
        index = _getIndex(indexId);
        if (!index)
            throw(LSST_EXCEPT(pexExcept::RuntimeErrorException,
                              boost::str(boost::format("Astrometry.net index with ID %i was not found") % indexId)));
    }

    fitstable_t* tag = startree_get_tagalong(index->starkd);
    if (!tag)
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          boost::str(boost::format("Astrometry.net index \"%s\": no tag-along table was found") % index->indexname)));
    
    std::vector<TagAlongColumn> columns;

    int N = fitstable_get_N_fits_columns(tag);
    for (int i=0; i<N; i++) {
        const char* colname = fitstable_get_fits_column_name(tag, i);
        if (!colname) {
            throw(LSST_EXCEPT(pexExcept::RuntimeErrorException,
                              boost::str(boost::format("Astrometry.net index \"%s\": couldn't get name for tag-along table column %i") % index->indexname % i)));
        }
        char* units;
        tfits_type fitstype;
        int arraysize;
        if (fitstable_find_fits_column(tag, colname, &units, &fitstype, &arraysize)) {
            throw(LSST_EXCEPT(pexExcept::RuntimeErrorException,
                              boost::str(boost::format("Astrometry.net index \"%s\": couldn't get tag-along table information for column \"%s\"") % index->indexname % colname)));
        }
        std::string ctype = "";
        switch (fitstype) {
        case TFITS_BIN_TYPE_D:
        case TFITS_BIN_TYPE_E:
            ctype = "Double";
            break;
        case TFITS_BIN_TYPE_I:
        case TFITS_BIN_TYPE_J:
        case TFITS_BIN_TYPE_A:
        case TFITS_BIN_TYPE_B:
            ctype = "Int";
            break;
        case TFITS_BIN_TYPE_K:
            ctype = "Int64";
            break;
        case TFITS_BIN_TYPE_L:
            ctype = "Bool";
            break;
        default:
            break;
        }

        TagAlongColumn col;
        col.name = std::string(colname);
        col.fitstype = fitstype;
        col.ctype = ctype;
        col.units = std::string(units);
        col.arraysize = arraysize;
        columns.push_back(col);
    }
    return columns;
}


std::vector<std::vector<double> > GlobalAstrometrySolution::getCatalogueExtra(afwGeom::Angle ra, afwGeom::Angle dec, afwGeom::Angle radius,
                                                                             std::vector<std::string> columns, int indexId) {
    index_t* index = NULL;
    for (unsigned int i=0; i<_indexList.size(); i++) {
        index_t* ind = _indexList[i];
        if (ind->indexid != indexId)
            continue;
        index = ind;
        break;
    }
    std::vector<std::vector<double> > x;
    if (!index) {
        return x;
    }

    int *starinds = NULL;
    int nstars = 0;
    double* radecs = NULL;
    double center[3];
    // this takes degrees
    radecdeg2xyzarr(ra.asDegrees(), dec.asDegrees(), center);
    double radius2 = radius.toUnitSphereDistanceSquared();
    startree_search_for(index->starkd, center, radius2, NULL, &radecs, &starinds, &nstars);

    std::vector<double> ras;
    std::vector<double> decs;
    // DEGREES
    for (int i=0; i<nstars; i++) {
        ras.push_back(radecs[2*i+0]);
        decs.push_back(radecs[2*i+1]);
    }
    x.push_back(ras);
    x.push_back(decs);

    for (std::vector<std::string>::iterator it = columns.begin(); it != columns.end(); it++) {
        x.push_back(getTagAlongFromIndex(index, *it, starinds, nstars));
    }
    free(starinds);
    return x;
}


///Returns a sourceSet of objects that are nearby in an raDec sense to the requested position. If
///filterName is not blank, we also extract out the magnitude information for that filter and 
///store (as a flux) in the returned SourceSet object. The value of filterName must match one of the
///strings returned by getCatalogueMetadataFields(). If you're not interested in fluxes, set
///filterName to ""
ReferenceSources
GlobalAstrometrySolution::getCatalogue(afwGeom::Angle ra,
                                       afwGeom::Angle dec,
                                       afwGeom::Angle radius,
                                       string filterName,
                                       string idName,
                                       int indexId) {

    double center[3];
    // degrees
    radecdeg2xyzarr(ra.asDegrees(), dec.asDegrees(), center);
    double radius2 = radius.toUnitSphereDistanceSquared();
    string msg;

    ReferenceSources refs;
    refs.indexid = -1;

    afwCoord::CoordSystem coordsys = _getCoordSys();

    for (unsigned int i=0; i<_indexList.size(); i++) {
        index_t* index = _indexList[i];
        if (indexId != -1 && index->indexid != indexId)
            continue;
        if (!index_is_within_range(index, ra.asDegrees(), dec.asDegrees(), radius.asDegrees()))
            continue;

        // Ensure the index is loaded...
        index_reload(index);

        //Find nearby stars
        double *xyz = NULL;
        int *starinds = NULL;
        int nstars = 0;
        startree_search_for(index->starkd, center, radius2, &xyz, NULL, &starinds, &nstars);

        if (nstars == 0)
            continue;

        vector<double> mag = getTagAlongFromIndex(index, filterName, starinds, nstars);
        vector<boost::int64_t> ids = getIds(idName, index, starinds, nstars);

        afwDet::SourceSet sources;
        std::vector<int> inds;

        // Create a source for every position stored
        for (int j = 0; j<nstars; ++j) {
            afwDet::Source::Ptr src(new afwDet::Source());
            src->setAllRaDecFields(xyztocoord(coordsys, xyz + j*3));
            inds.push_back(starinds[j]);
            if (mag.size()) {
                // convert mag to flux
                src->setPsfFlux( pow(10.0, -mag[j]/2.5) );
            }
            if (ids.size()) {
                src->setSourceId(ids[j]);
            }
            sources.push_back(src);
        }
        free(xyz);
        free(starinds);

        refs.indexid = index->indexid;
        refs.refsources = sources;
        refs.inds = inds;

        // NOTE change in behavior -- return only sources from first catalog...
        break;
    }
    return refs;
}


///A convenient interface to astrometry.net's startree_get_data_column. 
///Checks that the fieldName exists, and that something is returned.
///Returns an vector of length numIds
///
///\param index Astrometry.net index field to extract tagalong data from
///\param fieldName Name of tagalong column to extract. If this is empty (i.e ""), nothing is
///       extracted
///\param ids   Indices you want extracted. Get this from a match object, or startree_search_for
///\param numIds size of ids
static vector<double> getTagAlongFromIndex(index_t* index, string fieldName, int *ids, int numIds) {

    //Load magnitude information from catalogue
    double *tagAlong=NULL;
    string msg;
    
    if (fieldName != "") {
        // Grab tag-along data here. If it's not there, throw an exception
        if (! startree_has_tagalong(index->starkd) ) {
            msg = boost::str(boost::format("Index file \"%s\" has no metadata") % index->indexname);
            throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
        }

        for (int i=0; i<numIds; i++) {
            assert(ids[i] >= 0);
        }
        
        tagAlong = startree_get_data_column(index->starkd, fieldName.c_str(), ids, numIds);

        if( tagAlong == NULL) {
            msg = boost::str(boost::format("No tag-along column called %s found in index %s") %
                fieldName % index->indexname);
            throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
        }

        vector<double> out(&tagAlong[0], &tagAlong[numIds]);
        // the vector constructor makes a copy of the data.
        // http://www.sgi.com/tech/stl/Vector.html#1
        free(tagAlong);
        return out;
    }
    
    vector<double> out(0);
    return out;
    
}

///Plate scale of solution in arcsec/pixel. Note this is different than getMin(Max)ImageScale()
///which return the intial guesses of platescale.
lsst::afw::geom::Angle GlobalAstrometrySolution::getSolvedImageScale(){
    if (!_isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }
    MatchObj* match = getMatchObject();
    return match->scale * afwGeom::arcseconds;
} 

MatchObj* GlobalAstrometrySolution::getMatchObject() {
    return solver_get_best_match(_solver);
}

///Reset the object so it's ready to match another field.
void GlobalAstrometrySolution::reset() {
    if (_solver != NULL) {
        solver_free(_solver);
        _solver = solver_new();
    }
    if (_starxy != NULL) {
        starxy_free(_starxy);
        _starxy = NULL;
    }

    _numBrightObjects = USE_ALL_STARS_FOR_SOLUTION;
    _isSolved = false;
        
    //I should probably be smarter than this and remember the actual values of
    //the settings instead of just resetting the defaults
    setDefaultValues();
}
                                                         
}}}}
