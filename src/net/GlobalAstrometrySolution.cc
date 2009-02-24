

using namespace std;

#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"

//namespace net = lsst::meas::astrom::net;
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

    
GlobalAstrometrySolution::GlobalAstrometrySolution(
       lsst::afw::detection::SourceVector vec) { ///< Points indicating pixel coords of detected objects
    
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
void GlobalAstrometrySolution::addIndexFile(const std::string path) { ///< Path of index file
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
int GlobalAstrometrySolution::parseConfigFile(const std::string filename) { ///< Name of backend configuration file
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
void GlobalAstrometrySolution::setStarlist(lsst::afw::detection::SourceVector vec) ///<List of Sources
        throw(std::domain_error) {

    if (vec.empty()) {
        throw(domain_error("Src list contains no objects"));
    }

    //This number is conservative. A bare minimum of 4 objects is needed, although the
    //search probably won't be unique with that few objects
    if (vec.size() < 20) {
        throw(domain_error("Src list should contain at least 20 objects"));
    }
    
    int const size = vec.size(); 
    _starlist = starxy_new(size, true, false);   

    //Need to add flux information to _starlist, and sort by it.
    int i=0;
    for (lsst::afw::detection::SourceVector::iterator ptr = vec.begin(); ptr != vec.end(); ++ptr) {
        
        double const x = (*ptr)->getXAstrom();
        double const y = (*ptr)->getYAstrom();
        double const flux= (*ptr)->getPsfMag();
        
        starxy_set(_starlist, i, x, y);
        //There's no function to set the flux, so do it explicitly.
        //This would be a good improvement for the code
        _starlist->flux[i] = flux;
        ++i;
    }
    
    //Now sort the list and add to the solver
    starxy_sort_by_flux(_starlist);
    solver_set_field(_solver, _starlist);    
}


//
// Accessors
//
///    
///How good must a match be to be accepted.
double GlobalAstrometrySolution::getMatchThreshold(){
    return _solver->logratio_record_threshold;
}


double GlobalAstrometrySolution::getSolvedImageScale(){
    if (_solver == NULL) {
        throw(logic_error("No solution found yet. Did you run blindSolve()?"));
    }

    return(_solver->best_match.scale);
} 
        

///
///After solving, return a full Wcs including SIP distortion matrics
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getDistortedWcs() throw(logic_error){
    if (_solver == NULL) {
        throw(logic_error("No solution found yet. Did you run blindSolve()?"));
    }

    lsst::afw::image::PointD crpix(_solver->best_match.wcstan.crpix[0],
                                   _solver->best_match.wcstan.crpix[1]);
    lsst::afw::image::PointD crval(_solver->best_match.wcstan.crval[0],
                                   _solver->best_match.wcstan.crval[1]);

    //Linear conversion matrix
    int naxis = 2;   //This is hardcoded into the sip_t structure
    boost::numeric::ublas::matrix<double> CD(2,2);
    for (int i=0; i<naxis; ++i) {
        for (int j=0; j<naxis; ++j) {
            CD.insert_element(i, j, _solver->best_match.wcstan.cd[i][j]);
        }
    }

    //Forward distortion terms. In the SIP notation, these matrices are referred to
    //as A and B. I can find no documentation that insists that these matrices be
    //the same size, so I assume they aren't.
    int aSize = _solver->best_match.sip->a_order;
    boost::numeric::ublas::matrix<double> sipA(aSize, aSize);
    for (int i=0; i<aSize; ++i){
        for (int j=0; j<aSize; ++j){
            sipA.insert_element(i, j, _solver->best_match.sip->a[i][j]);
        }
    }

    //Repeat for B
    int bSize = _solver->best_match.sip->b_order;
    boost::numeric::ublas::matrix<double> sipB(bSize, bSize);
    for (int i=0; i<bSize; ++i){
        for (int j=0; j<bSize; ++j){
            sipB.insert_element(i, j, _solver->best_match.sip->b[i][j]);
        }
    }

    //Repeat for Ap, for the reverse transform
    int apSize = _solver->best_match.sip->ap_order;
    boost::numeric::ublas::matrix<double> sipAp(apSize, apSize);
    for (int i=0; i<apSize; ++i){
        for (int j=0; j<apSize; ++j){
            sipAp.insert_element(i, j, _solver->best_match.sip->ap[i][j]);
        }
    }

    //And finally, Bp, also part of the reverse transform
    int bpSize = _solver->best_match.sip->bp_order;
    boost::numeric::ublas::matrix<double> sipBp(bpSize, bpSize);
    for (int i=0; i<bpSize; ++i){
        for (int j=0; j<bpSize; ++j){
            sipBp.insert_element(i, j, _solver->best_match.sip->bp[i][j]);
        }
    }    
        
    
    lsst::afw::image::Wcs::Ptr wcsPtr = lsst::afw::image::Wcs::Ptr( new(lsst::afw::image::Wcs));
    *wcsPtr = lsst::afw::image::Wcs(crval, crpix, CD, sipA, sipB, sipAp, sipBp);
    
    return wcsPtr;   
}    

            

///              
///After solving, return a linear Wcs (i.e without distortion terms)
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getWcs() throw(logic_error) {
    if (_solver == NULL) {
        throw(logic_error("No solution found yet. Did you run blindSolve()?"));
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
///Convert ra dec to pixel coordinates    
pair<double, double> GlobalAstrometrySolution::raDec2Xy(double ra, double dec)  throw(std::logic_error) {
    if (! _solver->best_match_solves) {
        throw( logic_error("No solution found yet. Did you run blindSolve()?") );
    }

    double x, y;
    bool flag = tan_radec2pixelxy( &_solver->best_match.wcstan, ra, dec, &x, &y);

    //I don't think this conversion can ever fail
    assert(flag==true);

    pair<double, double> ret;
    ret.first = x;
    ret.second=y;
    return ret;
}
    
    
///    
///Convert pixels to right ascension declination    
pair<double, double> GlobalAstrometrySolution::xy2RaDec(double x, double y)  throw(std::logic_error) {
    if (! _solver->best_match_solves) {
        throw( logic_error("No solution found yet. Did you run blindSolve()?") );
    }

    double ra, dec;
    tan_pixelxy2radec( &_solver->best_match.wcstan, x, y, &ra, &dec);
    
    pair<double, double> ret;
    ret.first = ra;
    ret.second=dec;
    
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
        throw(logic_error("Starlist hasn't been set yet"));
    }
    
    int max;
    max = starxy_n(_solver->fieldxy);
    cout << "endobj is " << _solver->endobj << endl;
    max = (_solver->endobj == 0) ? max : _solver->endobj;
    
    for(int i=0; i<max; ++i){
        const double x = starxy_getx(_solver->fieldxy, i);
        const double y = starxy_gety(_solver->fieldxy, i);
        cout << x << ", " << y << endl;
    }
}

    

//
// Mutators
//
///Inform the solver that the image may suffer from some distortion. Turn this on if a purely linear
///wcs solution will get the postions of some stars wrong by more than 1 pixel
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
///attempted match. These values are more reasonable. Note that this function does not call
///solver_set_default_values() in solver.h, which is called automatically when you create a new
///solver object.
void GlobalAstrometrySolution::setDefaultValues() {

    //Set image scale boundaries (in arcseconds per pixel) to non-zero and non-infinity.
    //These values still exceed anything you'll find in a real image
    setMinimumImageScale(1e-6);
    setMaximumImageScale(3600*360);  //2pi radians per pixel

    //A typical image will have thousands of stars, but a match will proceed satisfactorily
    //with only the brightest subset. The number below is arbitary, but seems to work.
    setNumberStars(50);

    //Do we allow the solver to assume the image may have some distortion in it?
    allowDistortion(true);

    //How good must a match be to be considered good enough? Chosen by referring to
    //control-program.c
    setMatchThreshold(30);

    //This must be set to true to ensure data is loaded from the index files
    _backend->inparallel = true;
        
}
        
       
///Set the size of the image. At some point, when I know how to determine what my inputs are
///I may include this in the default constructor.
///Note that the solver will fail if min==max, so we make them different by a small amount.    
void GlobalAstrometrySolution::setImageScaleArcsecPerPixel(double scale) { ///< Plate scale of image
    setMinimumImageScale(.99*scale);
    setMaximumImageScale(1.01*scale);
}


///Set the verbosity level for astrometry.net. The higher the level the more information is returned.
///1 and 2 are typically good values to use.
void GlobalAstrometrySolution::setLogLevel(const int level) {
    if (level < 0 || level > 4) {
        throw( domain_error("Logging level must be between 0 and 4") );
    }
    
    log_init((enum log_level) level);
}


///    
///How good does a match need to be to be accepted. Typical value is log(1e12) approximately 27
void GlobalAstrometrySolution::setMatchThreshold(const double threshold) {
    _solver->logratio_record_threshold = threshold;
}
        

    

//
// Solve and verify functions
//

///Having correctly set up your GAS object, find the ra dec of your 
///image with no prior information. Returns true if a match was found,
///false otherwise
int GlobalAstrometrySolution::blindSolve() {

    //Throw exceptions if setup is correct
    if (! _solver->fieldxy) {
        throw logic_error("No field has been set yet");
    }

    //Check that minimum image scale != max scale, or the solver will do all the work and
    //then fail
    if( _solver->funits_lower >= _solver->funits_upper) {
        throw domain_error("Minimum image scale must be strictly less than max scale");
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


#if 0    
///    
///This function isn't working yet because I need to set _hprange first
bool GlobalAstrometrySolution::verifyRaDec(const double ra,   ///<Right ascension in decimal degrees
                                           const double dec   ///<Declination in decimal degrees
                                          ) throw(logic_error)  {
    
    assert(false);  //Make sure the function can't be called yet because it isn't working
    //Test that a field has been assigned
    if (! _solver-> fieldxy) {
        throw logic_error("No field has been set yet");
    }

    findNearbyIndices(ra, dec);
    solver_run( _solver);

    return(_solver->best_match_solves);
}


///This function isn't working yet    
///Verify that a previously calculated WCS (world coordinate solution) 
///matches the positions of the point sources in the image.
///Returns true if the Wcs matches
///This code is experimental, and probably won't compile
bool GlobalAstrometrySolution::verifyWcs(const lsst::afw::image::Wcs::Ptr wcsPtr) throw(logic_error) {
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

/// Astro.net's Wcs object consists of a very small subset of the 
/// elements of a wcsprm object, plus some fields that define a SIP
/// transformation (used to describe distortions in the image, see Shupe et al.
/// 2005). We ignore the SIP in this function, and fill in the overlapping
/// information as best we can.
sip_t* GlobalAstrometrySolution::convertWcsToSipt(const lsst::afw::image::Wcs::Ptr wcsPtr) throw(logic_error) {
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
        throw logic_error("Linear Transformation Matrix must be square");
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
void GlobalAstrometrySolution::findNearbyIndices(double ra,    ///<Right ascension in decimal degrees
                                                 double dec) { ///<Declination in decimal degrees
    //Leave this assertion in for the moment until I decide what the best way of setting hprange is.
    //Currently you have to call a separate function, which is unecessary work for the programmer
    assert( _hprange != 0);

    //Convert ra,dec to a 3 dimensional cartesian unit vector
    double centerxyz[3]; 
    radecdeg2xyzarr(ra, dec, centerxyz);

    //Each healpix represents a patch of sky, so we load the one pointed to by radec, and each
    //one around it
    int healpixes[9];

    //Linked list of arrays
    il* hplist = il_new(4);

    //Number of index files listed in the configuration file
    int nIndex = getNumIndices();
    for(int i=0; i<nIndex; ++i){
        index_t* index = (index_t *) pl_get( _backend->indexes, i);
        int nHealPixSides = index->meta.hpnside;

        int nHealpix = healpix_get_neighbours_within_range(centerxyz, _hprange, healpixes, nHealPixSides);

        il_append_array(hplist, healpixes, nHealpix);

        //If the index is nearby, add it
        if(il_contains(hplist, index->meta.healpix)){
            solver_add_index(_solver, index);
        }

        il_remove_all(hplist);
    }

    il_free(hplist);


}


}}}}
