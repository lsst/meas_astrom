#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"

using namespace std;
namespace Except = lsst::pex::exceptions;

namespace lsst { namespace meas { namespace astrom { namespace net {

//
// Constructors
//

GlobalAstrometrySolution::GlobalAstrometrySolution(const std::string policyPath) :
    _backend(NULL), _solver(NULL), _starxy(NULL), _numBrightObjects(-1)
{
 
    _backend  = backend_new();
    _solver   = solver_new();

    setDefaultValues();
    
    
    lsst::pex::policy::Policy pol(policyPath);
    _equinox = pol.getDouble("equinox");
    _raDecSys = pol.getString("raDecSys");
    std::string pkgDir = lsst::utils::eups::productDir("ASTROMETRY_NET_DATA");

    std::vector<std::string> indexArray = pol.getStringArray("indexFile");
    for(unsigned int i=0; i<indexArray.size(); ++i){
        cout << "Adding index file: " << pkgDir+"/"+indexArray[i] << endl;
        addIndexFile(pkgDir+"/"+indexArray[i]);
    }
}

    
//
// Destructor
//
GlobalAstrometrySolution::~GlobalAstrometrySolution() {

    if( _backend != NULL) {
        backend_free(_backend);
        _backend = NULL;
    }
    
    if(_solver != NULL){
        solver_free(_solver);
        _solver = NULL;
    }
}



//
// Initialisation (Mutators)
//
///    
///Adds a single index file to the backend.
void GlobalAstrometrySolution::addIndexFile(const std::string path ///< Path of index file
                                           ){
    //Copy a constant string into a non-const C style string
    //so it can be passed into a C function without complaint.
    int len = (int) path.length(); 
    char *fn = (char *) malloc((len+1)*sizeof(char));
    strncpy(fn, path.c_str(), len+1);

    //May have to add some error checking. add_index always returns zero regardless
    //of success or failure
    backend_add_index(_backend, fn);
    free(fn);
}

///Read a configuration file that contains a list of indices
///from a stream and load them into memory. Useful for debugging code, but not for production
///\param Name of backend configuration
int GlobalAstrometrySolution::parseConfigFile(const std::string filename) {
    //Copy a constant string into a non-const C style string
    //so it can be passed into a C function without complaint.
    int len = filename.length();
    char *fn = (char *) malloc((len+1)*sizeof(char));
    strncpy(fn, filename.c_str(), len);

    //Tells the backend to load all index files that it finds. Due to a 
    //bug in the astrometry.net 0.24 code, this line is required to load
    //*any* index files
    _backend->inparallel = true;

    return(backend_parse_config_file(_backend, fn));
    free(fn);
}

///
///Same as parseConfigFile except accepts a C FILE* pointer
int GlobalAstrometrySolution::parseConfigStream(FILE* fconf) {
    _backend->inparallel = true;
    return(backend_parse_config_file_stream(_backend, fconf));
}



///Set the image to be solved. The image is abstracted as a list of positions in pixel space
///Only sources with psfFlux > minFlux are added, and the number of sources added is returned.
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
}
    
//
// Accessors
//
///    
///How good must a match be to be accepted.
double GlobalAstrometrySolution::getMatchThreshold(){
    return _solver->logratio_record_threshold;
}


///Returns true is image is flipped, i.e the wcs solution has east pointed in the opposite sense (relative
///to north, than the astronomical convenction (which is north is up, east is left)
///Precondition: Image is already solved    
bool GlobalAstrometrySolution::isFlipped() {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    if(_solver->parity == FLIPPED_PARITY) {
        return(true);
    }
    return(false);
}

 
///Plate scale of solution in arcsec/pixel. Note this is different than getMin(Max)ImageScale()
///which return the intial guesses of platescale.
double GlobalAstrometrySolution::getSolvedImageScale(){
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    return(_solver->best_match.scale);
} 

///\brief Return a list of the stars used to solve the image.
///After solving an image, use this function to return the set of objects that was used to determine a
/// solution.
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
    bool isWeighted = true;
    int skip_shift = true;

    sip_t *sip = tweak_just_do_it(&_solver->best_match.wcstan, _starxy, 
                                  NULL,
                                  NULL, NULL,
                                  radec,
                                  nstars, jitter_arcsec, 
                                  order, inverse_order, iterations,
                                  isWeighted, skip_shift);


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
        
    
    lsst::afw::image::Wcs::Ptr wcsPtr(new lsst::afw::image::Wcs(crval, crpix, CD, 
        sipA, sipB, sipAp, sipBp, _equinox, _raDecSys));

    
    sip_free(sip);
    return wcsPtr;
}    

            

///              
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
///Convert ra dec to pixel coordinates assuming a linear Wcs.
///    
///In general, it is better to create a Wcs object using getDistortedWcs() and use it
///to perform conversions    
lsst::afw::image::PointD GlobalAstrometrySolution::raDecToXY(double ra, double dec) {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    double x, y;
    int flag = tan_radec2pixelxy( &_solver->best_match.wcstan, ra, dec, &x, &y);

    //I don't think this conversion can ever fail
    assert(flag==true);

    lsst::afw::image::PointD ret;
    ret[0] = x;
    ret[1] = y;
    return ret;
}
    
    
///    
///Convert pixels to right ascension, declination
///
///In general, it is better to create a Wcs object using getDistortedWcs() and use it
///to perform conversions  
lsst::afw::image::PointD GlobalAstrometrySolution::xyToRaDec(double x, double y)  {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    double ra, dec;
    tan_pixelxy2radec( &_solver->best_match.wcstan, x, y, &ra, &dec);
    
    lsst::afw::image::PointD ret;
    ret[0] = ra;
    ret[1] = dec;
    
    return ret ;
}

///    
/// Return the number of indices loaded from disk thus far   
int GlobalAstrometrySolution::getNumIndices() {
    return(pl_size(_backend->indexes));
}


///    
/// Return a list of the filenames of the Index files loaded into memory    
vector<string> GlobalAstrometrySolution::getIndexPaths() {
    vector<string> strList;
    int size=getNumIndices();

     index_meta_t *meta;
     meta = (index_meta_t*) malloc(sizeof(index_meta_t));
     for(int i=0; i<size; ++i){
         bl_get(_backend->indexmetas, i, meta);
         string x(meta->indexname);
         strList.push_back(x);
     }
     free(meta);

    return strList;
}
        

/// Print the coordinates and fluxes of all objects used to make a match. Primarily
/// a debugging function
void GlobalAstrometrySolution::printStarlist() {
    if ( ! _solver->fieldxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    
    int max;
    max = starxy_n(_solver->fieldxy);
    max = (_solver->endobj == 0) ? max : _solver->endobj;
    
    for(int i=0; i<max; ++i){
        const double x = starxy_getx(_solver->fieldxy, i);
        const double y = starxy_gety(_solver->fieldxy, i);
        cout << x << " " << y << " " << _solver->fieldxy->flux[i] << endl;
    }
}

    

//
// Mutators
//
///Inform the solver that the image may suffer from some distortion. Turn this on if a purely linear
///wcs solution will get the postions of some stars wrong by more than 1 pixel. Turned on by default
void GlobalAstrometrySolution::allowDistortion(bool distort) {
    _solver->distance_from_quad_bonus = (distort) ? true : false;
}

///Reset the objet so it's ready to match another field, but leave the backend alone because
///it doesn't need to be reset, and it is expensive to do so.
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

    //This must be set to true to ensure data is loaded from the index files
    _backend->inparallel = true;

    //From blind_run() in blind.c
    _solver->numtries=0;
    _solver->nummatches=0;
    _solver->numscaleok=0;
    _solver->num_cxdx_skipped=0;
    _solver->num_verified=0;
    _solver->quit_now = FALSE;

    setParity(UNKNOWN_PARITY);
}

///\brief Only find a solution using the brightest N objects
///Reducing the number of sources in the solution list can reduce the time taken to solve the image
///The truncated list is used by solve() to get the linear wcs, but all input sources are used when
///calculating the distortion terms of the SIP matrix.    
void GlobalAstrometrySolution::setNumBrightObjects(const int N) {

    if (_starxy == NULL) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No field set yet."));
    }

    if (N <= 0) {
        throw(LSST_EXCEPT(Except::RangeErrorException, "Illegal request. N must be greater than zero"));
    }
        
    if (N > starxy_n(_starxy) ) {
        throw(LSST_EXCEPT(Except::RangeErrorException, "Request exceeds number of available sources"));
    }

    _numBrightObjects = N;
}
    
       
///Set the plate scale of the image in arcsec per pixel
void GlobalAstrometrySolution::setImageScaleArcsecPerPixel(double scale ///< Plate scale of image
                                                          ) {
    //Note that the solver will fail if min==max, so we make them different by a small amount.    
    setMinimumImageScale(.99*scale);
    setMaximumImageScale(1.01*scale);
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
///North up and East right (or some rotation thereof) is parity==NORMAL_PARITY, the opposite is
/// parity==FLIPPED_PARITY.
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
// Solve and verify functions
//

///Having correctly set up your GAS object, find the ra dec of your 
///image with no prior information. Returns true if a match was found,
///false otherwise
bool GlobalAstrometrySolution::solve() {

    //Throw exceptions if setup is incorrect
    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }
    
    if(pl_size(_backend->indexes) == 0) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No indices have been set yet"));
    }
    
    if(_solver->best_match_solves){
        string msg="Solver indicated that a match has already been found. Do you need to reset?";
        throw(LSST_EXCEPT(Except::RuntimeErrorException, msg));
    }

    if( _solver->funits_lower >= _solver->funits_upper) {
        string msg = "Minimum image scale must be strictly less than max scale";
        throw(LSST_EXCEPT(Except::DomainErrorException, msg));
    }

    //Move all the indices from the backend structure to the solver structure.
    int size=pl_size(_backend->indexes);
    for(int i=0; i<size; ++i){
        index_t *index = (index_t *) pl_get(_backend->indexes, i);
        solver_add_index(_solver, index);
    }

    //Copies the brightest N objects from _starxy into _solver
    solverSetField();
    
    solver_run( _solver);

    return(_solver->best_match_solves);
}



///Find a solution with an initial guess at the position    
/// \param raDec Right ascension/declination in decimal degrees
bool GlobalAstrometrySolution::solve(const afw::image::PointD raDec   
                                          ) {
    return solve(raDec[0], raDec[1]);
}
    
///Find a solution with an initial guess at the position.
/// \param ra Right ascencsion in decimal degrees
/// \param dec Declination in decimal degrees
//Although the interface is the same as solve() it's actually quite a different function
bool GlobalAstrometrySolution::solve(const double ra, const double dec)  {    
    //Throw exceptions if setup is incorrect
    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }
    
    if(pl_size(_backend->indexes) == 0) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No indices have been set yet"));
    }
    
    if(_solver->best_match_solves){
        string msg="Solver indicated that a match has already been found. Do you need to reset?";
        throw(LSST_EXCEPT(Except::RuntimeErrorException, msg));
    }

    if( _solver->funits_lower >= _solver->funits_upper) {
        string msg = "Minimum image scale must be strictly less than max scale";
        throw(LSST_EXCEPT(Except::DomainErrorException, msg));
    }
    
    //Create a unit vector from the postion to be passed to loadNearbyIndices
    double xyz[3];
    radecdeg2xyzarr(ra, dec, xyz);
    vector<double> unitVector(3);
    for(int i=0; i<3; ++i){
        unitVector[i] = xyz[i];
    }

    solverSetField();      
    loadNearbyIndices(unitVector);
    solver_run( _solver);

    if(_solver->best_match_solves){
        logmsg("Position (%.7f %.7f) verified\n", ra, dec);
    }
    else {
        logmsg("Failed to verify position (%.7f %.7f)\n", ra, dec);
    }
    
    return(_solver->best_match_solves);
}


bool GlobalAstrometrySolution::solve(const lsst::afw::image::Wcs::Ptr wcsPtr, 
                                     const double imageScaleUncertaintyPercent) {

    //Rename the variable to something shorter to make the code easier to read
    double unc = imageScaleUncertaintyPercent/100.;

    //This test is strictly unecessary as solverSetField throws the same exception
    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    solverSetField();
    
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


///An "all-in-one" function. Take a starlist, and an initial guess at a Wcs, and do everything necessary
///to solve the field using the Wcs as an initial guess. The reset the solver so it's ready to be
///used again and return a pointer to the new Wcs
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::solve(const lsst::afw::detection::SourceSet vec, 
                                                           const lsst::afw::image::Wcs::Ptr wcsPtr) {

    setStarlist(vec);

    bool isSuccess = solve(wcsPtr);

    if(! isSuccess) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No solution found"));
    }

    lsst::afw::image::Wcs::Ptr wcsOutPtr = getDistortedWcs();

    //reset so we're ready to do the next match
    reset();

    return(wcsOutPtr);
}
    
    

//
// Private functions
//    

///Astrometry.net defines these in a weird place, and it seems easier to redefine them here
#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a))
#endif   

    
/// Astro.net's Wcs object consists of a very small subset of the 
/// elements of a wcsprm object, plus some fields that define a SIP
/// transformation (used to describe distortions in the image, see Shupe et al.
/// 2005). We ignore the SIP in this function, and fill in the overlapping
/// information as best we can.
sip_t* GlobalAstrometrySolution::convertWcsToSipt(const lsst::afw::image::Wcs::Ptr wcsPtr)  {
    sip_t *sip = sip_create();

    lsst::afw::image::PointD radec = (*wcsPtr).getOriginRaDec();
    sip->wcstan.crval[0] = radec.getX();
    sip->wcstan.crval[1] = radec.getY();

    lsst::afw::image::PointD rowCol = (*wcsPtr).getOriginXY();
    sip->wcstan.crpix[0] = rowCol.getY();    //Check this.
    sip->wcstan.crpix[1] = rowCol.getX();    //Check this.


    Eigen::Matrix2d CD = (*wcsPtr).getLinearTransformMatrix();
    //Test that a field has been assigned
    if (CD.rows() != CD.cols()) {
        throw(LSST_EXCEPT(Except::LengthErrorException, "Linear Transformation Matrix must be square"));
    }

    for (int i=0; i< CD.rows(); ++i) {
        for (int j=0; j< CD.cols(); ++j) {
            sip->wcstan.cd[i][j] = CD(i,j);
        }
    }

    return(sip);
}

///For a given position in the sky, figure out which of the (large) index
///fits files need to be read in from disk, and load them into the object
    void GlobalAstrometrySolution::loadNearbyIndices(std::vector<double> unitVector) {
    
    if(unitVector.size() != 3){
        throw(LSST_EXCEPT(Except::LengthErrorException, "Input vector should have exactly 3 elements") );
    }

    //This function should be called after solverSetField()
    if ( ! _solver->fieldxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "_solver->fieldxy hasn't been set yet"));
    }
    
    
    //Translate the vector into a C array that astrometry.net can understand
    double xyz[3];
    xyz[0] = unitVector[0];
    xyz[1] = unitVector[1];
    xyz[2] = unitVector[2];
    
    //Don't look at very small quads    
    double xSize = _solver->field_maxx - _solver->field_minx;
    double ySize = _solver->field_maxy - _solver->field_miny;
    setMinQuadScale(0.1*MIN(xSize, ySize));

    //What is the largest size we'll want to look at.
    //This function is defined in starutil.h
    double hprange = arcsec2dist(_solver->funits_upper*hypot(xSize,ySize)/2.);
    
    //Each index file contains blocks of healpixes, each containing 
    //a list of star positions of a given region of sky. Different index
    //files contain healpixes of different sized solid angles (i.e large 
    //regions , of sky, or small regions of sky). 
    //In any one index file there should be one healpix covering the position of
    //interest, and 8 surrounding ones. The area of interest will cover
    //one or more of those healpixes.
    int numIndex = pl_size(_backend->indexes);
    if (numIndex < 1) {   //Paranoia...
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No indices in backend"));
    }
        
    for(int i=0; i<numIndex; ++i){
        index_t* index = (index_t*) pl_get(_backend->indexes, i);
        index_meta_t* meta = &(index->meta);
        int healpixes[9];
        il* hplist = il_new(9);
        
        int nhp = healpix_get_neighbours_within_range(xyz, hprange, healpixes, meta->hpnside);

        //Decide which ones are of interest, and add those to the solver
        il_append_array(hplist, healpixes, nhp);
        if(il_contains(hplist, meta->healpix)) {
            logmsg("Adding index %s\n", meta->indexname);
            solver_add_index(_solver, index);
        }
        il_remove_all(hplist);
    }

    //Check that at least one index was added
    if ( pl_size(_solver->indexes) < 1) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No matching indices found for this field"));
    }
        
}    


void GlobalAstrometrySolution::solverSetField() {

    if ( ! _starxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    if( starxy_n(_starxy) == 0){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist has zero elements"));
    }
        

    //The default value, -1, indicates that all objects should be used
    if (_numBrightObjects == -1) {
        _numBrightObjects = starxy_n(_starxy);
    }

    int num = _numBrightObjects;  //Because I'm a lazy typist
    starxy_t *shortlist = starxy_new(num, true, true);
    for(int i=0; i<num; ++i) {
        double x = starxy_getx(_starxy, i);
        double y = starxy_gety(_starxy, i);
        double f = _starxy->flux[i];   //flux has no accessor function

        starxy_setx(shortlist, i, x);
        starxy_sety(shortlist, i, y);
        shortlist->flux[i] = f;
    }

    //Free the old structure stored in the solver, if necessary
    starxy_free(_solver->fieldxy);

    //Add the new, smaller, struture
    solver_set_field(_solver, shortlist);

    //Find field boundaries and precompute kdtree
    solver_preprocess_field(_solver);

}   

}}}}
