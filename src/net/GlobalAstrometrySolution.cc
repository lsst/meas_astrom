

using namespace std;

#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"

namespace Except = lsst::pex::exceptions;

namespace lsst { namespace meas { namespace astrom { namespace net {

//
// Constructors
//

GlobalAstrometrySolution::GlobalAstrometrySolution() {
 
    _backend  = backend_new();
    _solver   = solver_new();
    _starlist = NULL;

    _hprange = 0;

    setDefaultValues();
}

    
GlobalAstrometrySolution::GlobalAstrometrySolution(lsst::afw::detection::SourceSet vec ///< Points indicating pixel
                                                                                          ///<coords of detected
                                                                                          /// objects
                                                  ){
    _backend  = backend_new();
    _solver   = solver_new();
    _starlist = NULL;

    _hprange = 0;

    setDefaultValues();
    
    setStarlist(vec);
}



//
// Destructor
//
GlobalAstrometrySolution::~GlobalAstrometrySolution() {
    backend_free(_backend);
    if(_solver != NULL){
        solver_free(_solver);
    }

    if(_starlist != NULL){
        starxy_free(_starlist);
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
int GlobalAstrometrySolution::parseConfigFile(const std::string filename ///< Name of backend configuration file
                                             ) {
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
///
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
    
    //If this value was previously set, free the memory before reassigning.
    if(_starlist != NULL) {
        starxy_free(_starlist);
    }

    int const size = vec.size();
    _starlist = starxy_new(size, true, false);   

    //Need to add flux information to _starlist, and sort by it.
    int i=0;
    for (lsst::afw::detection::SourceSet::iterator ptr = vec.begin(); ptr != vec.end(); ++ptr) {
        
        double const x = (*ptr)->getXAstrom();
        double const y = (*ptr)->getYAstrom();
        double const flux= (*ptr)->getPsfFlux();
        
        starxy_set(_starlist, i, x, y);
        //There's no function to set the flux, so do it explicitly.
        //This would be a good improvement for the code
        _starlist->flux[i] = flux;
        ++i;
    }
    
    //Now sort the list and add to the solver
    //starxy_sort_by_flux(_starlist);
    solver_set_field(_solver, _starlist);

    //Find field boundaries and precompute kdtree
    solver_preprocess_field(_solver);
}


//
// Accessors
//
///    
///How good must a match be to be accepted.
double GlobalAstrometrySolution::getMatchThreshold(){
    return _solver->logratio_record_threshold;
}

 
///Plate scale of solution in arcsec/pixel. Note this is different than getMin(Max)ImageScale()
///which return the intial guesses of platescale.
double GlobalAstrometrySolution::getSolvedImageScale(){
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }

    return(_solver->best_match.scale);
} 
        

///
///After solving, return a full Wcs including SIP distortion matrics
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getDistortedWcs(int order)  {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
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
    double jitter_arcsec=1.0;
    int inverse_order = order;
    int iterations = 5;        //Taken from blind.c:628
    bool weighted = true;
    int skip_shift = true;

    sip_t *sip = tweak_just_do_it(&_solver->best_match.wcstan, _starlist, 
                                  NULL,
                                  NULL, NULL,
                                  radec,
                                  nstars, jitter_arcsec, 
                                  order, inverse_order, iterations,
                                  weighted, skip_shift);


    //Check that tweaking worked.
    if (sip == NULL) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"Tweaking failed"));
    }

    
    lsst::afw::image::PointD crpix(sip->wcstan.crpix[0],
                                   sip->wcstan.crpix[1]);
    lsst::afw::image::PointD crval(sip->wcstan.crval[0],
                                   sip->wcstan.crval[1]);

    //Linear conversion matrix
    int naxis = 2;   //This is hardcoded into the sip_t structure
    boost::numeric::ublas::matrix<double> CD(2,2);
    for (int i=0; i<naxis; ++i) {
        for (int j=0; j<naxis; ++j) {
            CD.insert_element(i, j, sip->wcstan.cd[i][j]);
        }
    }

    //Forward distortion terms. In the SIP notation, these matrices are referred to
    //as A and B. I can find no documentation that insists that these matrices be
    //the same size, so I assume they aren't.
    int aSize = sip->a_order;
    boost::numeric::ublas::matrix<double> sipA(aSize, aSize);
    for (int i=0; i<aSize; ++i){
        for (int j=0; j<aSize; ++j){
            sipA.insert_element(i, j, sip->a[i][j]);
        }
    }

    //Repeat for B
    int bSize = sip->b_order;
    boost::numeric::ublas::matrix<double> sipB(bSize, bSize);
    for (int i=0; i<bSize; ++i){
        for (int j=0; j<bSize; ++j){
            sipB.insert_element(i, j, sip->b[i][j]);
        }
    }

    //Repeat for Ap, for the reverse transform
    int apSize = sip->ap_order;
    boost::numeric::ublas::matrix<double> sipAp(apSize, apSize);
    for (int i=0; i<apSize; ++i){
        for (int j=0; j<apSize; ++j){
            sipAp.insert_element(i, j, sip->ap[i][j]);
        }
    }

    //And finally, Bp, also part of the reverse transform
    int bpSize = sip->bp_order;
    boost::numeric::ublas::matrix<double> sipBp(bpSize, bpSize);
    for (int i=0; i<bpSize; ++i){
        for (int j=0; j<bpSize; ++j){
            sipBp.insert_element(i, j, sip->bp[i][j]);
        }
    }    
        
    
    lsst::afw::image::Wcs::Ptr wcsPtr = lsst::afw::image::Wcs::Ptr( new(lsst::afw::image::Wcs));
    *wcsPtr = lsst::afw::image::Wcs(crval, crpix, CD, sipA, sipB, sipAp, sipBp);
    
    
    sip_free(sip);
    return wcsPtr;
}    

            

///              
///After solving, return a linear Wcs (i.e without distortion terms)
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getWcs()  {
    if (! _solver->best_match_solves) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException,"No solution found yet. Did you run solve()?"));
    }
    
    lsst::afw::image::PointD crpix(_solver->best_match.wcstan.crpix[0],
                                   _solver->best_match.wcstan.crpix[1]);
    lsst::afw::image::PointD crval(_solver->best_match.wcstan.crval[0],
                                   _solver->best_match.wcstan.crval[1]);
    
    int naxis = 2;   //This is hardcoded into the sip_t structure
    boost::numeric::ublas::matrix<double> CD(2,2);
    for (int i=0; i<naxis; ++i) {
        for (int j=0; j<naxis; ++j) {
            CD.insert_element(i, j, _solver->best_match.wcstan.cd[i][j]);
        }
    }

    lsst::afw::image::Wcs::Ptr wcsPtr = lsst::afw::image::Wcs::Ptr( new(lsst::afw::image::Wcs));
    *wcsPtr = lsst::afw::image::Wcs(crval, crpix, CD);
    
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
    bool flag = tan_radec2pixelxy( &_solver->best_match.wcstan, ra, dec, &x, &y);

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
        cout << x << " " << y << endl;
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
    starxy_free(_starlist);
    _starlist=NULL;

    solver_free(_solver);
    _solver = solver_new();
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

    //A typical image will have thousands of stars, but a match will proceed satisfactorily
    //with only the brightest subset. The number below is arbitary, but seems to work.
    _solver->startobj=0;
    setNumberStars(50);

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
///North up and East right (or some rotation thereof) is parity==0, the opposite is parity==1. If you
///don't know in advance, set parity==2 (the default)
void GlobalAstrometrySolution::setParity(const int parity){

    //Insist on legal values.
    switch (parity){
      case PARITY_BOTH:
      case PARITY_FLIP:
      case PARITY_NORMAL:
        _solver->parity = parity;
      default:
        throw(LSST_EXCEPT(Except::DomainErrorException, "Illegal parity setting\n"));
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
    if ( ! _solver->fieldxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }
    
    if(pl_size(_backend->indexes) == 0) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No indices have been set yet"));
    }
    
    if(_solver->best_match_solves){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Solver indicated that a match has already been found. Do you need to reset?"));
    }

    if( _solver->funits_lower >= _solver->funits_upper) {
        throw(LSST_EXCEPT(Except::DomainErrorException, "Minimum image scale must be strictly less than max scale"));
    }
    
    //Move all the indices from the backend structure to the solver structure.
    int size=pl_size(_backend->indexes);
    for(int i=0; i<size; ++i){
        index_t *index = (index_t *) pl_get(_backend->indexes, i);
        solver_add_index(_solver, index);
    }

    solver_run( _solver);

    return(_solver->best_match_solves);
}



///Find a solution with an initial guess at the position    
bool GlobalAstrometrySolution::solve(const afw::image::PointD raDec   ///<Right ascension/declination
                                               ///in decimal degrees
                                          ) {
    return solve(raDec[0], raDec[1]);
}
    
///Find a solution with an initial guess at the position    
bool GlobalAstrometrySolution::solve(const double ra, const double dec   ///<Right ascension/declination
                                               ///in decimal degrees
                                          )  {    
    //Throw exceptions if setup is incorrect
    if ( ! _solver->fieldxy) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Starlist hasn't been set yet"));
    }
    
    if(pl_size(_backend->indexes) == 0) {
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "No indices have been set yet"));
    }
    
    if(_solver->best_match_solves){
        throw(LSST_EXCEPT(Except::RuntimeErrorException, "Solver indicated that a match has already been found. Do you need to reset?"));
    }

    if( _solver->funits_lower >= _solver->funits_upper) {
        throw(LSST_EXCEPT(Except::DomainErrorException, "Minimum image scale must be strictly less than max scale"));
    }
    
    //Create a unit vector from the postion to be passed to loadNearbyIndices
    double xyz[3];
    radecdeg2xyzarr(ra, dec, xyz);
    vector<double> unitVector(3);
    for(int i=0; i<3; ++i){
        unitVector[i] = xyz[i];
    }

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

#if 0    
//This function isn't working yet    
//Verify that a previously calculated WCS (world coordinate solution) 
//matches the positions of the point sources in the image.
//Returns true if the Wcs matches
//This code is experimental, and probably won't compile
bool GlobalAstrometrySolution::verifyWcs(const lsst::afw::image::Wcs::Ptr wcsPtr)  {
    //Test that a field has been assigned
    if (! _solver-> fieldxy) {
        throw logic_error("No field has been set yet");
    }

    assert(false);  //Make sure the function can't be called because it isn't working yet
    //Extract information from the WCS object and stuff it into
    //a tan_t structure from astrometry.net
    sip_t *sip = convertWcsToSipt(wcsPtr);

     //Find all indices that are near the position indicated in the Wcs
    double ra = sip->wcstan.crval[0];
    double dec = sip->wcstan.crval[1];
    findNearbyIndices(ra, dec);

    solver_verify_sip_wcs(_solver, sip);
    sip_free(sip);

    return(_solver->best_match_solves);
}

    #endif


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


    boost::numeric::ublas::matrix<double> CD = (*wcsPtr).getLinearTransformMatrix();
    //Test that a field has been assigned
    if (CD.size1() != CD.size2()) {
        throw(LSST_EXCEPT(Except::LengthErrorException, "Linear Transformation Matrix must be square"));
    }

    for (unsigned int i=0; i< CD.size1(); ++i) {
        for (unsigned int j=0; j< CD.size2(); ++j) {
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
    
    //Translate the vector into a C array that astrometry.net can understand
    double xyz[3];
    xyz[0] = unitVector[0];
    xyz[1] = unitVector[1];
    xyz[2] = unitVector[2];
    
    //Don't look at very small quads    
    //FIXME: use _solver->field_max[x,y] instead of recomputing
    double x1,x2, y1, y2;
    x1 = x2 = _solver->fieldxy->x[0];
    y1 = y2 = _solver->fieldxy->y[0];    
    for(int i=0; i< _solver->fieldxy->N; ++i){
        x1 = MIN(x1, _solver->fieldxy->x[i]);
        x2 = MAX(x2, _solver->fieldxy->x[i]);

        y1 = MIN(y1, _solver->fieldxy->y[i]);
        y2 = MAX(y2, _solver->fieldxy->y[i]);
    }
    double xSize = x2-x1;
    double ySize = y2-y1;
    _solver->quadsize_min = 0.1*MIN(xSize, ySize);

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
    int N = pl_size(_backend->indexes);
    
    for(int i=0; i<N; ++i){
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
}    




}}}}
