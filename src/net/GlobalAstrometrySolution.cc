// -*- LSST-C++ -*-

#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
#include <cstdio>
#include <iostream>

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace net {


using namespace std;
namespace Except = lsst::pex::exceptions;
namespace Det = lsst::afw::detection;
namespace pexLog = lsst::pex::logging;

int const USE_ALL_STARS_FOR_SOLUTION = -1;
//
//Constructors, Destructors
//
GlobalAstrometrySolution::GlobalAstrometrySolution(const std::string policyPath):
    _mylog(pexLog::Log::getDefaultLog(), "meas.astrom.net", pexLog::Log::DEBUG),
    _indexList(NULL), 
    _solver(NULL), 
    _starxy(NULL), 
    _numBrightObjects(USE_ALL_STARS_FOR_SOLUTION) {
    
    _solver   = solver_new();

    setDefaultValues();
    
    lsst::pex::policy::Policy pol(policyPath);
    _equinox = pol.getDouble("equinox");
    _raDecSys = pol.getString("raDecSys");
    string pkgDir = lsst::utils::eups::productDir("ASTROMETRY_NET_DATA");

    //Add meta information about every index listed in the policy file
    std::vector<std::string> indexArray = pol.getStringArray("indexFile");


    _mylog.log(pexLog::Log::DEBUG, "Loading meta information on indices...");    
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
    _mylog.log(pexLog::Log::DEBUG, "Meta information loaded...");    
}


index_t *GlobalAstrometrySolution::_loadIndexMeta(std::string filename){
  //return index_load(filename.c_str(), INDEX_ONLY_LOAD_METADATA, NULL);
  double t0 = timenow();
  index_t* index = index_load(filename.c_str(), INDEX_ONLY_LOAD_METADATA, NULL);
  double dt = timenow() - t0;
  string msg = boost::str(boost::format("loading index file %s took %s sec.") % filename % dt);
  _mylog.log(pexLog::Log::DEBUG, msg);
  return index;
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


///astrometry.net intialises the solver with some default values that guarantee failure in any
///attempted match. These values are more reasonable. 
void GlobalAstrometrySolution::setDefaultValues() {

    //Among other things, this sets the parity, positional uncertainty (_solver->verify_pix)                 <
    //matching accuracy (_solver->codetol                    
    solver_set_default_values(_solver);
    
    //Set image scale boundaries (in arcseconds per pixel) to non-zero and non-infinity.
    //These values still exceed anything you'll find in a real image
    setMinimumImageScale(1e-6);
    setMaximumImageScale(3600*360);  //2pi radians per pixel

    //Do we allow the solver to assume the image may have some distortion in it?
    allowDistortion(true);

    //How good must a match be to be considered good enough?  Log-odds.
    setMatchThreshold(log(1e12));

    // Reset counters and record of best match found so far.
    solver_cleanup_field(_solver);

    setParity(UNKNOWN_PARITY);
}

    
//
// Setup functions
//

///Set the image to be solved. The image is abstracted as a list of positions in pixel space
void GlobalAstrometrySolution::setStarlist(lsst::afw::detection::SourceSet vec ///<List of Sources
                                          ) {
    if (vec.empty()) {
        throw(LSST_EXCEPT(Except::LengthErrorException, "Src list contains no objects"));
    }

    if (_starxy != NULL) {
        starxy_free(_starxy);
    }
    int const size = vec.size();
    _starxy = starxy_new(size, true, false);   

    int i = 0;
    for (lsst::afw::detection::SourceSet::iterator ptr = vec.begin(); ptr != vec.end(); ++ptr) {
        double const x    = (*ptr)->getXAstrom();
        double const y    = (*ptr)->getYAstrom();
        double const flux = (*ptr)->getPsfFlux();

        starxy_set(_starxy, i, x, y);
        starxy_set_flux(_starxy, i, flux);
        ++i;
    }

    //Sort the array
    starxy_sort_by_flux(_starxy);
    _solverSetField();
}


///\brief Only find a solution using the brightest N objects
///Reducing the number of sources in the solution list can reduce the time taken to solve the image
///The truncated list is used by solve() to get the linear wcs, but all input sources are used when
///calculating the distortion terms of the SIP matrix.    
void GlobalAstrometrySolution::setNumBrightObjects(int N) {

    if (N <= 0) {
        throw(LSST_EXCEPT(Except::RangeErrorException, "Illegal request. N must be greater than zero"));
    }

    _numBrightObjects = N;
    
    if (_starxy != NULL){
       _solverSetField();
    }
}


void GlobalAstrometrySolution::_solverSetField() {

    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    int starxySize = starxy_n(_starxy);
    if ( starxySize == 0){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist has zero elements"));
    }

    int N = _numBrightObjects;  //Because I'm a lazy typist
    
    //The default value, -1, indicates that all objects should be used
    if (N == USE_ALL_STARS_FOR_SOLUTION) {
        N = starxySize;
    }

    if (N > starxySize) {
        string msg = "numBrightObjects set to a larger value than number of stars";
        throw(LSST_EXCEPT(Except::RuntimeErrorException, msg));
    }

    starxy_t *shortlist = starxy_subset(_starxy, N);
    assert(shortlist);

    //Set the pointer in the solver to the new, smaller field
    starxy_free(_solver->fieldxy);
    solver_set_field(_solver, shortlist);
    solver_reset_field_size(_solver);

    //Find field boundaries and precompute kdtree
    solver_preprocess_field(_solver);
}   

///Set the plate scale of the image in arcsec per pixel
void GlobalAstrometrySolution::setImageScaleArcsecPerPixel(double imgScale ///< Plate scale of image
                                                          ) {
    //Note that the solver will fail if min==max, so we make them different by a small amount.    
    setMinimumImageScale(0.99*imgScale);
    setMaximumImageScale(1.01*imgScale);
}


///Inform the solver that the image may suffer from some distortion. Turn this on if a purely linear
///wcs solution will get the postions of some stars wrong by more than 1 pixel. Turned on by default
void GlobalAstrometrySolution::allowDistortion(bool hasDistortion) {
    _solver->distance_from_quad_bonus = (hasDistortion) ? true : false;
}

///Set the verbosity level for astrometry.net. The higher the level the more information is returned.
///1 and 2 are typically good values to use. 4 will print so much to the screen that it slows execution
void GlobalAstrometrySolution::setLogLevel(int level) {
    if (level < 0 || level > 4) {
        throw(LSST_EXCEPT(Except::DomainErrorException, "Logging level must be between 0 and 4"));
    }
    
    log_init((enum log_level) level);
}


///    
///How good does a match need to be to be accepted. Typical value is log(1e12) approximately 27
void GlobalAstrometrySolution::setMatchThreshold(double threshold) {
    _solver->logratio_record_threshold = threshold;
}


///You can double the speed of a match if you know the parity, i.e whether the image is flipped or not.
///North up and East right (or some rotation thereof) is parity==NORMAL_PARITY, the opposite is
///parity==FLIPPED_PARITY. The default is UNKNOWN_PARITY
void GlobalAstrometrySolution::setParity(int parity){

    //Insist on legal values.
    switch (parity){
        case UNKNOWN_PARITY:
        case FLIPPED_PARITY:
        case NORMAL_PARITY:
            _solver->parity = parity;
            break;
            return;
        default:
            throw LSST_EXCEPT(Except::DomainErrorException, "Illegal parity setting");
    }
}


//
// Solve functions
//

///Solve using a wcs object as an initial guess
bool GlobalAstrometrySolution::solve(const lsst::afw::image::Wcs::Ptr wcsPtr, 
                                     double imageScaleUncertaintyPercent) {

    //Rename the variable to something shorter to make the code easier to read
    double unc = imageScaleUncertaintyPercent/100.0;

    //This test is strictly unecessary as solverSetField throws the same exception
    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    double xc = (_solver->field_maxx + _solver->field_minx)/2.0;
    double yc = (_solver->field_maxy + _solver->field_miny)/2.0;

    //Get the central ra/dec and plate scale
    lsst::afw::image::PointD raDec = wcsPtr->xyToRaDec(xc, yc);
    double plateScaleArcsecPerPixel = sqrt(wcsPtr->pixArea(raDec))*3600;
    setMinimumImageScale(plateScaleArcsecPerPixel*(1 - unc));
    setMaximumImageScale(plateScaleArcsecPerPixel*(1 + unc));
    
    
    string msg = boost::str( boost::format("Solving using initial guess at position of\n %.7f %.7f\n") % raDec.getX() % raDec.getY());
    _mylog.log(pexLog::Log::DEBUG, msg);

    double lwr = plateScaleArcsecPerPixel*(1 - unc);
    double upr = plateScaleArcsecPerPixel*(1 + unc);
    msg = boost::str( boost::format("Scale range %.3f - %.3f arcsec/pixel\n") % lwr % upr);
    _mylog.log(pexLog::Log::DEBUG, msg);
    

    if ( wcsPtr->isFlipped()) {
        setParity(FLIPPED_PARITY);
        _mylog.log(pexLog::Log::DEBUG, "Setting Flipped parity");        
    } else {
        setParity(NORMAL_PARITY);
        _mylog.log(pexLog::Log::DEBUG, "Setting Normal parity");        
    }
        
    return(solve(raDec.getX(), raDec.getY()));
}


///Find a solution with an initial guess at the position    
bool GlobalAstrometrySolution::solve(const afw::image::PointD raDec   ///<Right ascension/declination
                                               ///in decimal degrees
                                          ) {
    return solve(raDec[0], raDec[1]);
}


    
///Find a solution with an initial guess at the position.
bool GlobalAstrometrySolution::solve(double ra,   ///<Right ascension in decimal degrees
                                     double dec   ///< Declination in decimal degrees
                                          )  {    
    string msg;

    // Tell the solver to only consider matches within the image size of the supposed RA,Dec.
    // The magic number 2.0 out front says to accept matches within 2 radii of the given *center* position.
    double maxRadius = 2.0 * arcsec2deg(_solver->funits_upper * _solver->field_diag / 2.0);
    msg = boost::str(boost::format("Setting RA,Dec = (%g, %g), radius = %g deg") % ra % dec % maxRadius);
    _mylog.log(pexLog::Log::DEBUG, msg);
    solver_set_radec(_solver, ra, dec, maxRadius);

    int success = _callSolver(ra, dec);

    if (success){
        char *indexname = _solver->index->indexname;
        msg = boost::str(boost::format("Position verified. Solved index is %s") % indexname);
        _mylog.log(pexLog::Log::DEBUG, msg);        
        
    }
    else {
        msg = boost::str(boost::format("Failed to verify position (%.7f %.7f)\n") % ra % dec);
    }

    _mylog.log(pexLog::Log::DEBUG, msg);            
    return(success);
}


///Find a solution blindly, with no initial guess. Go get a cup of tea, this function
///will take a while
bool GlobalAstrometrySolution::solve()  {    

    // Don't use any hints about the RA,Dec position.
    solver_clear_radec(_solver);

    int success = _callSolver(NO_POSITION_SET, NO_POSITION_SET);

    string msg;    
    if (success){
        _mylog.log(pexLog::Log::DEBUG, "Position Found");

        // Grab everything we need from the index file while it is still open!
        index_t* index = _solver->best_match.index;

        // FIXME -- refradec, fieldxy, tweak, tagalong.

    }
    else {
        _mylog.log(pexLog::Log::DEBUG, "Failed");
    }

    // Unload all index files?
    
    return(success);
}


///Check that all the setup was done correctly, then set the solver to work
///By default, _callSolver does a blind solve, unless the optional ra and dec
///are given values
bool GlobalAstrometrySolution::_callSolver(double ra, double dec) {
    //Throw exceptions if setup is incorrect
    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    if (_indexList.size() == 0) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No index files loaded yet"));
    }
    
    if (_solver->best_match_solves){
        string msg = "Solver indicated that a match has already been found. Do you need to reset?";
        throw(LSST_EXCEPT(Except::RuntimeErrorException, msg));
    }

    if ( _solver->funits_lower >= _solver->funits_upper) {
        string msg = "Minimum image scale must be strictly less than max scale";
        throw(LSST_EXCEPT(Except::DomainErrorException, msg));
    }

    //Calculate the best guess at image size
    double xSizePixels = solver_field_width(_solver);
    double ySizePixels = solver_field_height(_solver);
    double minSizePixels = min(xSizePixels, ySizePixels);
    double maxSizePixels = max(xSizePixels, ySizePixels);
    assert(maxSizePixels > minSizePixels && minSizePixels > 0);

    //Set the range of sizes of quads to examine
    //@FIXME the 10% and 90% should be parameters
    double imgSizeArcSecLwr = .10 * _solver->funits_lower*minSizePixels;
    double imgSizeArcSecUpr = .90 * _solver->funits_upper*maxSizePixels;

    //Output some useful debugging info
    string msg;
    msg = boost::str(boost::format("Image size %.0f x %.0f pixels") % xSizePixels % ySizePixels);
    _mylog.log(pexLog::Log::DEBUG, msg);
    
    msg = boost::str(boost::format("Platescale is %.3f -- %.3f arcsec/pixel") \
        % _solver->funits_lower % _solver->funits_upper);
    _mylog.log(pexLog::Log::DEBUG, msg);
    

    msg = boost::str(boost::format("Image size between %g and %g arcsec") 
        % (_solver->funits_lower*minSizePixels) % (_solver->funits_upper*maxSizePixels)  );
    _mylog.log(pexLog::Log::DEBUG, msg);

    _mylog.log(pexLog::Log::DEBUG, "Setting indices");
    _addSuitableIndicesToSolver(imgSizeArcSecLwr, imgSizeArcSecUpr, ra, dec);

    

    _mylog.log(pexLog::Log::DEBUG, "Doing solve step");

    //solver_print_to(_solver, stdout);

    solver_run(_solver);

    return(_solver->best_match_solves);
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
int GlobalAstrometrySolution::_addSuitableIndicesToSolver(double imgSizeArcSecLwr, double imgSizeArcSecUpr, double ra, double dec) {

    bool hasAtLeastOneIndexOfSuitableScale = false;
    int nMeta = _indexList.size();
    int nSuitable = 0;
    bool blind = (ra == NO_POSITION_SET) || (dec == NO_POSITION_SET);

    string msg;

    for (int i = 0; i<nMeta; ++i){
        index_t* index = _indexList[i];

        if (!index_overlaps_scale_range(index, imgSizeArcSecLwr, imgSizeArcSecUpr))
            continue;

        hasAtLeastOneIndexOfSuitableScale = true;

        //Is this either a blind solve, or does the index cover a suitable
        //patch of sky
        if (!(blind || index_is_within_range(index, ra, dec, arcsec2deg(imgSizeArcSecUpr))))
            continue;

        // Found a good one!

        // Load the data (not just metadata), if it hasn't been already...
        index_reload(index);

	msg = boost::str(boost::format("Adding index %s") % index->indexname);
	_mylog.log(pexLog::Log::DEBUG, msg);

        solver_add_index(_solver, index);
        nSuitable++;
    }
    
    if(nSuitable == 0) {
        string msg = "No suitable indices found for given input parameters:";
        
        if(hasAtLeastOneIndexOfSuitableScale) {
            msg += "Probably the ra/dec range isn't covered";
        }
        else {
            msg += "No indices of a suitable scale were found";
        }
        throw(LSST_EXCEPT(Except::RuntimeErrorException, msg));
    }
    
    return nSuitable;
}


//
//Return the solution
//              

///After solving, return a linear Wcs (i.e without distortion terms)
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getWcs()  {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    // CHECK THIS -- Astrometry.net probably doesn't add or subtract 1 from your coordinates;
    // this if you pass in zero-indexed source positions, that's what you'll get back.

    ///Astro.net conforms with wcslib in assuming that images are 1-indexed (i.e the bottom left-most pixel
    ///is (1,1). LSST is zero indexed, so we add 1 to the crpix values returned by _solver to convert
    lsst::afw::image::PointD crpix(_solver->best_match.wcstan.crpix[0] + 1,
                                   _solver->best_match.wcstan.crpix[1] + 1);   
    lsst::afw::image::PointD crval(_solver->best_match.wcstan.crval[0],
                                   _solver->best_match.wcstan.crval[1]);
    
    int naxis = 2;   //This is hardcoded into the sip_t structure
    Eigen::Matrix2d CD;
    for (int i = 0; i<naxis; ++i) {
        for (int j = 0; j<naxis; ++j) {
            CD(i, j) = _solver->best_match.wcstan.cd[i][j];
        }
    }

    lsst::afw::image::Wcs::Ptr wcsPtr(new lsst::afw::image::Wcs(crval, crpix, CD, _equinox, _raDecSys)); 
    return wcsPtr;
}


///
///After solving, return a full Wcs including SIP distortion matrics
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getDistortedWcs(int order)  {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    if (! _starxy){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist isn't set"));
    }

    //Generate an array of radec of positions in the field
    MatchObj* mo = &_solver->best_match;

    //Call the tweaking algorthim to generate the distortion coeffecients
    
    //jitter is a measure of how much we can expect the xy of stars to scatter from the expected
    //radec due to noise in our measurments.
    double jitterArcsec = tan_pixel_scale(&mo->wcstan)*_solver->verify_pix;
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
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Tweaking failed"));
    }

    ///Astro.net conforms with wcslib in assuming that images are 1-indexed (i.e the bottom left-most pixel
    ///is (1,1). LSST is zero indexed, so we add 1 to the crpix values returned by _solver to convert
    lsst::afw::image::PointD crpix(sip->wcstan.crpix[0] + 1,
                                   sip->wcstan.crpix[1] + 1);
    lsst::afw::image::PointD crval(sip->wcstan.crval[0],
                                   sip->wcstan.crval[1]);

    //Linear conversion matrix
    int naxis = 2;   //This is hardcoded into the sip_t structure
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
        
    
    lsst::afw::image::Wcs::Ptr wcsPtr(new lsst::afw::image::Wcs(crval, crpix, CD, 
                                      sipA, sipB, sipAp, sipBp, _equinox, _raDecSys));

    
    sip_free(sip);
    return wcsPtr;
}    


///\brief Return a list of the stars used to solve the image.
///After solving an image, use this function to return the set of objects that was used to determine a
///solution. Typically this list will be about 4 or 5 objects long. The ra dec of each object is 
///accessed using src.getRa() and src.getDec(), while the chip coords are accessed with  getXAstrom(),
///getYAstrom()
lsst::afw::detection::SourceSet GlobalAstrometrySolution::getMatchedSources(){
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }
    
    lsst::afw::detection::SourceSet set;

    MatchObj* match = &_solver->best_match;

    // Grab tag-along data.
    startree_t* skdt = match->index->starkd;
    double* umag = startree_get_data_column(skdt, "u", match->refstarid, match->nindex);
    double* gmag = startree_get_data_column(skdt, "g", match->refstarid, match->nindex);
    double* rmag = startree_get_data_column(skdt, "r", match->refstarid, match->nindex);
    double* imag = startree_get_data_column(skdt, "i", match->refstarid, match->nindex);
    double* zmag = startree_get_data_column(skdt, "z", match->refstarid, match->nindex);
    double* uerr = startree_get_data_column(skdt, "u_err", match->refstarid, match->nindex);
    double* gerr = startree_get_data_column(skdt, "g_err", match->refstarid, match->nindex);
    double* rerr = startree_get_data_column(skdt, "r_err", match->refstarid, match->nindex);
    double* ierr = startree_get_data_column(skdt, "i_err", match->refstarid, match->nindex);
    double* zerr = startree_get_data_column(skdt, "z_err", match->refstarid, match->nindex);
    double* reddening = startree_get_data_column(skdt, "reddening", match->refstarid, match->nindex);

    for (int i=0; i<match->nfield; i++) {
        // "theta" is the mapping from image (aka field) stars to index (aka reference) stars.
        // negative means no match.
        if (match->theta[i] < 0)
            continue;
        lsst::afw::detection::Source::Ptr ptr(new lsst::afw::detection::Source());
        ptr->setXAstrom(starxy_getx(_solver->fieldxy, i));
        ptr->setYAstrom(starxy_gety(_solver->fieldxy, i));

        double ra, dec;
        //This function is defined in astrometry.net. It converts the position of the star
        //as a three dimensional unit vector to ra,dec.
        xyzarr2radecdeg(match->refxyz + match->theta[i]*3, &ra, &dec);
        ptr->setRa(ra);
        ptr->setDec(dec);

        // int I = match->theta[i]
        // ptr->setMags(umag[I], gmag[I], rmag[I], imag[I], zmag[I])
        // ptr->setMagErrs(uerr[I], gerr[I], rerr[I], ierr[I], zerr[I])
        // ptr->setReddening(reddening[I]);  // e(b-v)

        set.push_back(ptr);
    }

    free(umag); free(gmag); free(rmag); free(imag); free(zmag);
    free(uerr); free(gerr); free(rerr); free(ierr); free(zerr);
    free(reddening);

    return set;
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
        _mylog.log(pexLog::Log::DEBUG, "No metadata found for index");        
        return output;
    }
    
    //Don't free this pointer, it points to memory that may be used later
    fitstable_t *table = startree_get_tagalong(_indexList[0]->starkd);
    assert(table != NULL);
    
    sl *nameList = fitstable_get_fits_column_names(table, NULL);

    for(int i=0; i< sl_size(nameList); ++i) {
        string name(sl_pop(nameList));
        output.push_back(name);
    }
    
    sl_free2(nameList);
        
    return output;
}    


///Returns a sourceSet of objects that are nearby in an raDec sense to the best match solution
lsst::afw::detection::SourceSet 
GlobalAstrometrySolution::getCatalogue(double radiusInArcsec, string filterName) {

    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    //Don't free this pointer. It points to memory that will still exist
    //when this function goes out of scope.
    double *center = _solver->best_match.center;
    double ra, dec;
    xyzarr2radecdeg(center, &ra, &dec);
    
    return getCatalogue(ra, dec, radiusInArcsec, field);
}


///Returns a sourceSet of objects that are nearby in an raDec sense to the requested position. If
///filterName is not blank, we also extract out the magnitude information for that filter and 
///store (as a flux) in the returned SourceSet object. The value of filterName must match one of the
///strings returned by getCatalogueMetadataFields(). If you're not interested in fluxes, set
///filterName to ""
lsst::afw::detection::SourceSet GlobalAstrometrySolution::getCatalogue(double ra,
    double dec,
    double radiusInArcsec,
    string filterName) {

    double center[3];
    radecdeg2xyzarr(ra, dec, center);
    double radius2 = arcsec2distsq(radiusInArcsec);
    string msg;

    Det::SourceSet out;


    for (unsigned int i=0; i<_indexList.size(); i++) {
        index_t* index = _indexList[i];
        if (!index_is_within_range(index, ra, dec, arcsec2deg(radiusInArcsec)))
            continue;

        // Ensure the index is loaded...
        index_reload(index);

        //Find nearby stars
        double *radec = NULL;
        int *starinds = NULL;
        int nstars = 0;
        startree_search_for(index->starkd, center, radius2, NULL, &radec, &starinds, &nstars);

        double *mag = NULL;
        if (metaName != "") {
            // Grab tag-along data here. If it's not there, throw an exception
            if (! startree_has_tagalong(index->starkd) ) {
                msg = boost::str(boost::format("Index file \"%s\" has no metadata") % index->indexname);
                throw(LSST_EXCEPT(Except::RuntimeErrorException, msg));
            }
            mag = startree_get_data_column(index->starkd, metaName.c_str(), starinds, nstars);

            if (mag == NULL) {
                msg = boost::str(boost::format("No meta data called %s found in index %s") %
                    metaName % index->indexname);
                throw(LSST_EXCEPT(Except::RuntimeErrorException, msg));            }
        }

        //Create a source for every position stored
        for (int j = 0; j<nstars; ++j) {
            Det::Source::Ptr ptr(new lsst::afw::detection::Source());

            ptr->setRa(radec[2*j]);
            ptr->setDec(radec[2*j + 1]);
            
            if(mag != NULL) {   //convert mag to flux
                ptr->setPsfFlux( exp10(-mag[j]/2.5) );
            }

            out.push_back(ptr);
        }
        free(radec);    
        free(mag);
    }

    return out;

}



///Plate scale of solution in arcsec/pixel. Note this is different than getMin(Max)ImageScale()
///which return the intial guesses of platescale.
double GlobalAstrometrySolution::getSolvedImageScale(){
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    return(_solver->best_match.scale);
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
        
    //I should probably be smarter than this and remember the actual values of
    //the settings instead of just resetting the defaults
    setDefaultValues();

}
                                                         
}}}}
