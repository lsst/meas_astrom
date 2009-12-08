
#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
#include <cstdio>
#include <iostream>

namespace lsst { namespace meas { namespace astrom { namespace net {
using namespace std;
namespace Except = lsst::pex::exceptions;
namespace Det = lsst::afw::detection;
namespace pexLog = lsst::pex::logging;

//
//Constructors, Destructors
//
GlobalAstrometrySolution::GlobalAstrometrySolution(const std::string policyPath):
    _mylog(pexLog::Log::getDefaultLog(), "meas.astrom.net", pexLog::Log::DEBUG),
    _indexList(NULL), _metaList(NULL), _solver(NULL), _starxy(NULL), _numBrightObjects(-1)
{
 
    

    _indexList = pl_new(sizeof(index_t));
    _metaList = pl_new(sizeof(index_meta_t));
    _solver   = solver_new();

    setDefaultValues();
    
    
    lsst::pex::policy::Policy pol(policyPath);
    _equinox = pol.getDouble("equinox");
    _raDecSys = pol.getString("raDecSys");
    string pkgDir = lsst::utils::eups::productDir("ASTROMETRY_NET_DATA");

    //Add meta information about every index listed in the policy file
    std::vector<std::string> indexArray = pol.getStringArray("indexFile");


    _mylog.log(pexLog::Log::DEBUG, "Loading meta information on indices...");    
    for(unsigned int i=0; i<indexArray.size(); ++i){
        string path = pkgDir+"/"+indexArray[i];
        pl_push(_metaList, _loadIndexMeta(path));
    }
    _mylog.log(pexLog::Log::DEBUG, "Meta information loaded...");    
}


index_meta_t *GlobalAstrometrySolution::_loadIndexMeta(string filename)
{
    int errValue=-2;
    index_meta_t *val = (index_meta_t*) malloc(sizeof(index_meta_t));
    
    
    //Copy a constant string into a non-const C style string
    //so it can be passed into a C function without complaint.
    int len = (int) filename.length(); 
    val->indexname = (char *) malloc((len+1)*sizeof(char));
    strncpy(val->indexname, filename.c_str(), len+1);

    
    qfits_header *hdr0 = qfits_header_readext(val->indexname, 0);
    if(hdr0 == NULL)
    {
        fprintf(stderr, "Error reading %s\n", val->indexname);
        return NULL;
    }
        
    val->index_scale_upper = qfits_header_getdouble(hdr0, "scale_u", (double) errValue);
    val->index_scale_lower = qfits_header_getdouble(hdr0, "scale_l", (double) errValue);
    val->indexid = qfits_header_getint(hdr0, "indexid", errValue);
    val->healpix = qfits_header_getint(hdr0, "healpix", errValue);
    val->nquads = qfits_header_getint(hdr0, "nquads", errValue);
    val->nstars = qfits_header_getint(hdr0, "nstars", errValue);
    //These two have default values. At least according to quadfile.c:
    //callback_read_header
    val->hpnside = qfits_header_getint(hdr0, "hpnside", 1);
    val->dimquads = qfits_header_getint(hdr0, "dimquads", 4);
    
    qfits_header_destroy(hdr0);

        
    qfits_header *hdr2 = qfits_header_readext(val->indexname, 2);
    if(hdr0 == NULL)
    {
        fprintf(stderr, "Error reading %s\n", val->indexname);
        return NULL;
    }
    
    val->circle = qfits_header_getboolean(hdr2, "CIRCLE", 0);
    val->cx_less_than_dx = qfits_header_getboolean(hdr2, "CXDX", 0);
    qfits_header_destroy(hdr2);

    //Sanity checking. These should never occur
    assert(val->index_scale_upper != (double) errValue);
    assert(val->index_scale_lower != (double) errValue);
    assert(val->indexid != errValue);
    assert(val->healpix != errValue); 
    assert(val->nquads != errValue);
    assert(val->nstars != errValue); 
    
    //This may occur in older index files
    if(!val->circle)
    {
        fprintf(stderr, "Code kdtree does not contain CIRCLE header card\n");
        return NULL;
    }
    
    return val;
}


GlobalAstrometrySolution::~GlobalAstrometrySolution() {

    if( _indexList != NULL) {
        pl_free_elements(_indexList);
        pl_free(_indexList);
        _indexList = NULL;
    }

    if( _metaList != NULL) {
        pl_free_elements(_metaList);
        pl_free(_metaList);
        _metaList = NULL;
    }

    if( _starxy != NULL) {
        starxy_free(_starxy);
        _starxy = NULL;
    }
    
    if(_solver != NULL){
        solver_free(_solver);
        _solver = NULL;
    }
}


///astrometry.net intialises the solver with some default values that garauntee failure in any
///attempted match. These values are more reasonable. 
void GlobalAstrometrySolution::setDefaultValues() {

    solver_set_default_values(_solver);
    
    //Set image scale boundaries (in arcseconds per pixel) to non-zero and non-infinity.
    //These values still exceed anything you'll find in a real image
    setMinimumImageScale(1e-6);
    setMaximumImageScale(3600*360);  //2pi radians per pixel

    //Do we allow the solver to assume the image may have some distortion in it?
    allowDistortion(true);

    //How good must a match be to be considered good enough? Chosen by referring to
    //control-program.c
    setMatchThreshold(30);


    //From blind_run() in blind.c
    _solver->numtries=0;
    _solver->nummatches=0;
    _solver->numscaleok=0;
    _solver->num_cxdx_skipped=0;
    _solver->num_verified=0;
    _solver->quit_now = FALSE;

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

    //This number is conservative. A bare minimum of 4 objects is needed, although the
    //search probably won't be unique with that few objects
    if (vec.size() < 20) {
        throw(LSST_EXCEPT(Except::LengthErrorException, "Src list should contain at least 20 objects"));
    }
    
    int const size = vec.size();
    if(_starxy != NULL) {
        starxy_free(_starxy);
    }
    _starxy = starxy_new(size, true, false);   

    int i=0;
    for (lsst::afw::detection::SourceSet::iterator ptr = vec.begin(); ptr != vec.end(); ++ptr) {
        
        double const x = (*ptr)->getXAstrom();
        double const y = (*ptr)->getYAstrom();
        double const flux= (*ptr)->getPsfFlux();

        starxy_set(_starxy, i, x, y);
        //There's no function to set the flux, so do it explicitly.
        //This would be a good improvement for the astrometry.net code
        _starxy->flux[i] = flux;
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
void GlobalAstrometrySolution::setNumBrightObjects(const int N) {

    if (N <= 0) {
        throw(LSST_EXCEPT(Except::RangeErrorException, "Illegal request. N must be greater than zero"));
    }

    _numBrightObjects = N;
    
    if(_starxy != NULL){
       _solverSetField();
    }
}


void GlobalAstrometrySolution::_solverSetField() {

    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    int starxySize = starxy_n(_starxy);
    if( starxySize == 0){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist has zero elements"));
    }
        

    int N = _numBrightObjects;  //Because I'm a lazy typist
    
    //The default value, -1, indicates that all objects should be used
    if (N == -1) {
        N = starxySize;
    }

    if(N > starxySize) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "numBrightObjects set to a larger value than number of stars"));
    }
    

    starxy_t *shortlist = starxy_new(N, true, true);
    for(int i=0; i<N; ++i) {
        double x = starxy_getx(_starxy, i);
        double y = starxy_gety(_starxy, i);
        double f = _starxy->flux[i];   //flux has no accessor function

        starxy_setx(shortlist, i, x);
        starxy_sety(shortlist, i, y);
        shortlist->flux[i] = f;
    }

    //Set the pointer in the solver to the new, smaller field
    starxy_free(_solver->fieldxy);
    solver_set_field(_solver, shortlist);

    //Find field boundaries and precompute kdtree
    solver_preprocess_field(_solver);

}   

///Set the plate scale of the image in arcsec per pixel
void GlobalAstrometrySolution::setImageScaleArcsecPerPixel(double imgScale ///< Plate scale of image
                                                          ) {
    //Note that the solver will fail if min==max, so we make them different by a small amount.    
    setMinimumImageScale(.99*imgScale);
    setMaximumImageScale(1.01*imgScale);
}


///Inform the solver that the image may suffer from some distortion. Turn this on if a purely linear
///wcs solution will get the postions of some stars wrong by more than 1 pixel. Turned on by default
void GlobalAstrometrySolution::allowDistortion(bool distort) {
    _solver->distance_from_quad_bonus = (distort) ? true : false;
}

///Set the verbosity level for astrometry.net. The higher the level the more information is returned.
///1 and 2 are typically good values to use. 4 will print so much to the screen that it slows execution
void GlobalAstrometrySolution::setLogLevel(const int level) {
    if (level < 0 || level > 4) {
        throw(LSST_EXCEPT(Except::DomainErrorException, "Logging level must be between 0 and 4"));
    }
    
    log_init((enum log_level) level);
}


///    
///How good does a match need to be to be accepted. Typical value is log(1e12) approximately 27
void GlobalAstrometrySolution::setMatchThreshold(const double threshold) {
    _solver->logratio_record_threshold = threshold;
}


///You can double the speed of a match if you know the parity, i.e whether the image is flipped or not.
///North up and East right (or some rotation thereof) is parity==NORMAL_PARITY, the opposite is parity==FLIPPED_PARITY.
///The default is UNKNOWN_PARITY
void GlobalAstrometrySolution::setParity(const int parity){

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
///Find a solution with an initial guess at the position    
bool GlobalAstrometrySolution::solve(const afw::image::PointD raDec   ///<Right ascension/declination
                                               ///in decimal degrees
                                          ) {
    return solve(raDec[0], raDec[1]);
}


    
///Find a solution with an initial guess at the position.
bool GlobalAstrometrySolution::solve(const double ra,   ///<Right ascension in decimal degrees
                                     const double dec   ///< Declination in decimal degrees
                                          )  {    
    //Throw exceptions if setup is incorrect
    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }
    
    if(pl_size(_metaList) == 0) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No index metas loaded yet"));
    }
    
    if(_solver->best_match_solves){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Solver indicated that a match has already been found. Do you need to reset?"));
    }

    if( _solver->funits_lower >= _solver->funits_upper) {
        throw(LSST_EXCEPT(Except::DomainErrorException, "Minimum image scale must be strictly less than max scale"));
    }

    int nMeta = pl_size(_metaList);
    for(int i=0; i<nMeta; ++i){
        //Which index does this meta point to?
        index_meta_t *meta = (index_meta_t*) pl_get(_metaList, i);
        int metaId = meta->indexid;
        int metaHealpix = meta->healpix;
        int metaHpnside = meta->hpnside;

        //Only look at indices that cover the appropriate region of sky, and the right
        //range of image scales.
        if( _isIndexMetaPossibleMatch(meta, ra, dec)) {       
                        
            //Have we already loaded this index from disk?
            int nIndex = pl_size(_indexList);
            index_t *trialIndex = NULL;
            
            for(int j=0; j<nIndex && trialIndex == NULL; ++j) {
                trialIndex = (index_t*) pl_get(_indexList, j);
                assert(trialIndex != NULL);

                //Three values uniquely identify an index
                bool isEqual = (metaId == trialIndex->meta.indexid);
                isEqual = isEqual && (metaHealpix == trialIndex->meta.healpix);
                isEqual = isEqual &&(metaHpnside == trialIndex->meta.hpnside);
                if(!isEqual) {
                    trialIndex = NULL;
                }
            }
            
            //If not loaded already, read from disk. 
            //This is a potentially slow operation
            if(trialIndex == NULL) {
                trialIndex = index_load(meta->indexname, 0);
                pl_push(_indexList, trialIndex);
            }
            solver_add_index(_solver, trialIndex);
        }
    }

    solver_run(_solver);
            
    if(_solver->best_match_solves){
        logmsg("Position (%.7f %.7f) verified\n", ra, dec);
    }
    else {
        logmsg("Failed to verify position (%.7f %.7f)\n", ra, dec);
    }
    
    return(_solver->best_match_solves);
}



///Find a solution blindly, with not initial guess. Go get a cup of tea, this function
///will take a while
bool GlobalAstrometrySolution::solve()  {    
    //Throw exceptions if setup is incorrect
    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }
    
    if(pl_size(_metaList) == 0) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No index metas loaded yet"));
    }
    
    if(_solver->best_match_solves){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Solver indicated that a match has already been found. Do you need to reset?"));
    }

    if( _solver->funits_lower >= _solver->funits_upper) {
        throw(LSST_EXCEPT(Except::DomainErrorException, "Minimum image scale must be strictly less than max scale"));
    }

    int nMeta = pl_size(_metaList);
    for(int i=0; i<nMeta; ++i){
        //Which index does this meta point to?
        index_meta_t *meta = (index_meta_t*) pl_get(_metaList, i);
        int metaId = meta->indexid;
        int metaHealpix = meta->healpix;
        int metaHpnside = meta->hpnside;

        //Only look at indices that cover the appropriate region of sky, and the right
        //range of image scales.
        if( _isMetaSuitableScale(meta)){  
            //Have we already loaded this index from disk?
            int nIndex = pl_size(_indexList);
            index_t *trialIndex = NULL;
            for(int j=0; j<nIndex && trialIndex == NULL; ++j) {
                trialIndex = (index_t*) pl_get(_indexList, j);
                assert(trialIndex != NULL);

                //Three values uniquely identify an index
                bool isEqual = (metaId == trialIndex->meta.indexid);
                isEqual = isEqual && (metaHealpix == trialIndex->meta.healpix);
                isEqual = isEqual &&(metaHpnside == trialIndex->meta.hpnside);
                if(!isEqual) {
                    trialIndex = NULL;
                }
            }

            //If not loaded already, read from disk. 
            //This is a potentially slow operation
            if(trialIndex == NULL) {
                trialIndex = index_load(meta->indexname, 0);
                pl_push(_indexList, trialIndex);
            }
            solver_add_index(_solver, trialIndex);
        }
    }

    solver_run(_solver);
            
    if(_solver->best_match_solves){
        logmsg("Position found\n");
    }
    else {
        logmsg("Failed\n");
    }
    
    return(_solver->best_match_solves);
}


bool GlobalAstrometrySolution::solve(const lsst::afw::image::Wcs::Ptr wcsPtr, const double imageScaleUncertaintyPercent) {

    //Rename the variable to something shorter to make the code easier to read
    double unc = imageScaleUncertaintyPercent/100.;

    //This test is strictly unecessary as solverSetField throws the same exception
    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    double xc = (_solver->field_maxx + _solver->field_minx)/2.0;
    double yc = (_solver->field_maxy + _solver->field_miny)/2.0;

    //Get the central ra/dec and plate scale
    lsst::afw::image::PointD raDec = wcsPtr->xyToRaDec(xc, yc);
    lsst::afw::image::PointD raDec2 = wcsPtr->xyToRaDec(xc+1, yc);
    double plateScale = hypot(raDec2.getX()-raDec.getX(), raDec2.getY()-raDec.getY());  //In degrees
    plateScale *= 3600;    //In arcseconds per pixel

    setMinimumImageScale(plateScale*(1-unc));
    setMaximumImageScale(plateScale*(1+unc));

    if( wcsPtr->isFlipped()) {
        setParity(FLIPPED_PARITY);
    } else {
        setParity(NORMAL_PARITY);
    }
        
    return(solve(raDec.getX(), raDec.getY()));
}


///Returns true if a meta points to an index at an appropriate scale for the image being solved
///and points to the correct region of the sky.
bool GlobalAstrometrySolution::_isIndexMetaPossibleMatch(index_meta_t *meta, double ra, double dec) {
    return  _isMetaNearby(meta, ra, dec) && _isMetaSuitableScale(meta);
}


///Does this meta point to an index that is close to the initial guess at position in radec space
bool GlobalAstrometrySolution::_isMetaNearby(index_meta_t *meta, double ra, double dec) {

    //-1 => an all sky index
    if(meta->healpix == -1)
    {   return 1;
    }

    //Create a unit vector from the postion 
    double xyz[3];
    radecdeg2xyzarr(ra, dec, xyz);
    vector<double> unitVector(3);
    for(int i=0; i<3; ++i){
        unitVector[i] = xyz[i];
    }

    //Size of image in units of the unit sphere
    double xSize = _solver->field_maxx - _solver->field_minx;
    double ySize = _solver->field_maxy - _solver->field_miny;
    //This function is defined in starutil.h
    double imgSize = arcsec2dist(_solver->funits_upper*hypot(xSize,ySize)/2.);
    
    
    //ids of nearby healpixes
    int hpArray[9];
    int nhp = healpix_get_neighbours_within_range(xyz, imgSize, hpArray, meta->hpnside);

    //Is meta->healpix one of the returned ids?
    for(int i=0; i<nhp; ++i) {   
        if(meta->healpix == hpArray[i]){
           return true;
        }
    }
    return false;

}


///Does this index contain asterisms that are similar in size (in arcsec) to the expected
///size of the image. We can speed up the match by ignoring files that contain only very
///large, or very small asterisms.
bool GlobalAstrometrySolution::_isMetaSuitableScale(index_meta_t *meta){

    //Size of image in pixels
    double xSize = _solver->field_maxx - _solver->field_minx;
    double ySize = _solver->field_maxy - _solver->field_miny;
    double imgSize = hypot(xSize,ySize)/2.;
    
    //Maximum possible size of the image. We don't want to look at asterisms bigger than this.
    double maxSizeArcsec = _solver->funits_upper*imgSize;
    
    //Neither do we want to look at asterisms that are significantly smaller than the
    //size of the image. 
    double minSizeArcsec = .1*_solver->funits_lower*imgSize;

    //Convert index scales from radians to arcsec
    double iLwr = meta->index_scale_lower*180./M_PI*3600.;
    double iUpr = meta->index_scale_upper*180./M_PI*3600.;
    
    if(iLwr < maxSizeArcsec && iUpr > minSizeArcsec){
        return true;
    }
    
    return true;    //Debugging: should be false.
}


//
//Return the solution
//              

///After solving, return a linear Wcs (i.e without distortion terms)
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getWcs()  {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    ///Astro.net conforms with wcslib in assuming that images are 1-indexed (i.e the bottom left-most pixel
    ///is (1,1). LSST is zero indexed, so we add 1 to the crpix values returned by _solver to convert
    lsst::afw::image::PointD crpix(_solver->best_match.wcstan.crpix[0]+1,
                                   _solver->best_match.wcstan.crpix[1]+1);   
    lsst::afw::image::PointD crval(_solver->best_match.wcstan.crval[0],
                                   _solver->best_match.wcstan.crval[1]);
    
    int naxis = 2;   //This is hardcoded into the sip_t structure
    Eigen::Matrix2d CD;
    for (int i=0; i<naxis; ++i) {
        for (int j=0; j<naxis; ++j) {
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
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    if(! _starxy){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist isn't set"));
    }

    //Generate an array of radec of positions in the field

    //radius of bounding circle of healpix of best match
    double radius = _solver->best_match.radius;
    radius = 1.05*radius;  //Add in a margin of safety
    double radius2 = radius*radius;
    
    //Don't free this pointer. It points to memory that will still exist
    //when this function goes out of scope.
    double *center = _solver->best_match.center;

    //Generate array
    double *radec=NULL;
    int nstars=0;
    startree_search(_solver->index->starkd, center, radius2, NULL, &radec, &nstars);

    //Call the tweaking algorthim to generate the distortion coeffecients

    //jitter is a measure of how much we can expect the xy of stars to scatter from the expected
    //radec due to noise in our measurments.
    double jitter_arcsec = tan_pixel_scale(&_solver->best_match.wcstan) * _solver->verify_pix;
    jitter_arcsec = hypot(jitter_arcsec, _solver->index->meta.index_jitter);
    int inverse_order = order;
    int iterations = 5;        //blind.c:628 uses 5
    bool weighted = true;
    int skip_shift = true;

    sip_t *sip = tweak_just_do_it(&_solver->best_match.wcstan, _starxy, 
                                  NULL,
                                  NULL, NULL,
                                  radec,
                                  nstars, jitter_arcsec, 
                                  order, inverse_order, iterations,
                                  weighted, skip_shift);
    free(radec);

    //Check that tweaking worked.
    if (sip == NULL) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"Tweaking failed"));
    }

    ///Astro.net conforms with wcslib in assuming that images are 1-indexed (i.e the bottom left-most pixel
    ///is (1,1). LSST is zero indexed, so we add 1 to the crpix values returned by _solver to convert
    lsst::afw::image::PointD crpix(sip->wcstan.crpix[0]+1,
                                   sip->wcstan.crpix[1]+1);
    lsst::afw::image::PointD crval(sip->wcstan.crval[0],
                                   sip->wcstan.crval[1]);

    //Linear conversion matrix
    int naxis = 2;   //This is hardcoded into the sip_t structure
    Eigen::Matrix2d CD;
    for (int i=0; i<naxis; ++i) {
        for (int j=0; j<naxis; ++j) {
            CD(i,j) = sip->wcstan.cd[i][j];
        }
    }


    //Forward distortion terms. In the SIP notation, these matrices are referred to
    //as A and B. I can find no documentation that insists that these matrices be
    //the same size, so I assume they aren't.
    int aSize = sip->a_order+1;
    Eigen::MatrixXd sipA(aSize, aSize);
    for (int i=0; i<aSize; ++i){
        for (int j=0; j<aSize; ++j){
            sipA(i, j) = sip->a[i][j];
        }
    }

    //Repeat for B
    int bSize = sip->b_order+1;
    Eigen::MatrixXd sipB(bSize, bSize);
    for (int i=0; i<bSize; ++i){
        for (int j=0; j<bSize; ++j){
            sipB(i, j) = sip->b[i][j];
        }
    }

    //Repeat for Ap, for the reverse transform
    int apSize = sip->ap_order+1;
    Eigen::MatrixXd sipAp(apSize, apSize);
    for (int i=0; i<apSize; ++i){
        for (int j=0; j<apSize; ++j){
            sipAp(i, j) = sip->ap[i][j];
        }
    }

    //And finally, Bp, also part of the reverse transform
    int bpSize = sip->bp_order+1;
    Eigen::MatrixXd sipBp(bpSize, bpSize);
    for (int i=0; i<bpSize; ++i){
        for (int j=0; j<bpSize; ++j){
            sipBp(i, j) = sip->bp[i][j];
        }
    }    
        
    
    lsst::afw::image::Wcs::Ptr wcsPtr(new lsst::afw::image::Wcs(crval, crpix, CD, sipA, sipB, sipAp, sipBp, _equinox, _raDecSys));

    
    sip_free(sip);
    return wcsPtr;
}    


///\brief Return a list of the stars used to solve the image.
///After solving an image, use this function to return the set of objects that was used to determine a solution.
///Typically this list will be about 4 or 5 objects long. The ra dec of each object is accessed using
///src.getRa() and src.getDec(), while the chip coords are accessed with  getXAstrom(), getYAstrom()
lsst::afw::detection::SourceSet GlobalAstrometrySolution::getMatchedSources(){
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    
    lsst::afw::detection::SourceSet set;

    for(int i =0; i< _solver->best_match.dimquads; ++i) {
        lsst::afw::detection::Source::Ptr ptr(new lsst::afw::detection::Source());

        ptr->setXAstrom(_solver->best_match.quadpix[2*i]);
        ptr->setYAstrom(_solver->best_match.quadpix[2*i + 1]);

        //This function is defined in astrometry.net. It converts the position of the star
        //as a three dimensional unit vector to ra,dec.
        double ra, dec;
        xyzarr2radecdeg(&_solver->best_match.quadxyz[i*3], &ra, &dec);
        ptr->setRa(ra);
        ptr->setDec(dec);

        set.push_back(ptr);
    }

    return set;
}


///Plate scale of solution in arcsec/pixel. Note this is different than getMin(Max)ImageScale()
///which return the intial guesses of platescale.
double GlobalAstrometrySolution::getSolvedImageScale(){
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    return(_solver->best_match.scale);
} 



///Returns a sourceSet of objects that are nearby in an raDec sense to the best match solution
///This function is still under construction and doesn't do what I want it to yet.
lsst::afw::detection::SourceSet GlobalAstrometrySolution::getCatalogue(double radiusInArcsec) {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    //Don't free this pointer. It points to memory that will still exist
    //when this function goes out of scope.
    double *center = _solver->best_match.center;


    //Generate an array of radec of positions in the field

    //radius of bounding circle of healpix of best match
    double radius2 = arcsec2distsq(radiusInArcsec);
    
    Det::SourceSet out;
    

    //For each index that we've search, pulled out stars that are close in an radec sense    
    for(int i=0; i< pl_size(_solver->indexes); ++i) {
        double *radec=NULL;
        int nstars=0;

        index_t* index = (index_t *) pl_get(_solver->indexes, i);
        startree_search(index->starkd, center, radius2, NULL, &radec, &nstars);
        printf("%i %i %i\n", index->meta.indexid, index->meta.healpix, nstars);

        //Create a  source for every position stored
        for(int j=0; j<nstars; ++j) {
            Det::Source::Ptr ptr(new lsst::afw::detection::Source());

            ptr->setRa(radec[2*j]);
            ptr->setDec(radec[2*j+1]);

            out.push_back(ptr);
        }

        free(radec);    
    }
    return out;

}

///Reset the object so it's ready to match another field.
void GlobalAstrometrySolution::reset() {

    if(_solver != NULL) {
        solver_free(_solver);
        _solver = solver_new();
    }

    
    if(_starxy != NULL) {
        starxy_free(_starxy);
        _starxy= NULL;
    }

    _numBrightObjects = -1;
        
    //I should probably be smarter than this and remember the actual values of
    //the settings instead of just resetting the defaults
    setDefaultValues();

}
                                                         
}}}}
