

using namespace std;

#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"

//namespace net = lsst::meas::astrom::net;
namespace lsst { namespace meas { namespace astrom { namespace net {

//
// Constructors
//


GlobalAstrometrySolution::GlobalAstrometrySolution(const std::string filename) { ///< Name of backend configuration file

    _backend  = backend_new();
    _solver   = solver_new();
    _starlist = NULL;
    _sip = NULL;

    _hprange = 0;

    parseConfigFile(filename);
}



GlobalAstrometrySolution::GlobalAstrometrySolution(const std::string filename,  ///< Name of backend configuration file
                                                   std::vector<lsst::afw::detection::Source::Ptr> src) { ///< Points
    ///<indicating pixel coords of detected objects
    _backend  = backend_new();
    _solver   = solver_new();
    _starlist = NULL;
    _sip = NULL;

    _hprange = 0;

    parseConfigFile(filename);
    setStarlist(src);

}



//
// Destructor
//
GlobalAstrometrySolution::~GlobalAstrometrySolution() {
    backend_free(_backend);
    solver_free(_solver);

    if(_starlist != NULL){
        starxy_free(_starlist);
    }
}



//
// Initialisation (Mutators)
//

///Read a configuration file that contains a list of indices
///from a stream and load them into memory
int GlobalAstrometrySolution::parseConfigFile(const std::string filename) { ///< Name of backend configuration file
    //Copy a constant string into a non-const C style string
    //so it can be passed into a C function without complaint.
    int len = filename.length();
    char *fn = (char *) malloc((len+1)*sizeof(char));
    strncpy(fn, filename.c_str(), len);
    cout << "Converted " << filename << "  to " << fn << endl;

    //Tells the backend to load all index files that it finds. Due to a 
    //bug in the astrometry.net code, this line is required to load
    //*any* index files
    _backend->inparallel = true;

    return(backend_parse_config_file(_backend, fn));
    free(fn);
}


///Same as parseConfigFile except accepts a C FILE* pointer
int GlobalAstrometrySolution::parseConfigStream(FILE* fconf) {
    _backend->inparallel = true;
    return(backend_parse_config_file_stream(_backend, fconf));
}



// ********************************************************************
///Init an object using a vector of sources
void GlobalAstrometrySolution::setStarlist(std::vector<lsst::afw::detection::Source::Ptr> src) throw(domain_error) {
    if (src.empty()) {
        throw(domain_error("Src list contains no objects"));
    }

    //This number is conservative. A bare minimum of 4 objects is needed, although the
    //search probably won't be unique with that few objects
    if (src.size() < 20) {
        throw(domain_error("Src list should contain at least 20 objects"));
    }

    int const size = src.size(); 
    _starlist = starxy_new(size, true, false);   

    int i=0;
    for (std::vector<lsst::afw::detection::Source::Ptr>::iterator ptr = src.begin(); ptr != src.end(); ++ptr) {
        lsst::afw::detection::Source::Ptr obj = *ptr;
        double const col = obj->getColc();
        double const row = obj->getRowc();
        starxy_set(_starlist, i++, col, row);
    }

    solver_set_field(_solver, _starlist);
}


void GlobalAstrometrySolution::tmpStarlist(std::vector<lsst::afw::detection::Source> src) {
    cout << "Cp" << endl;

}


//
// Accessors
//
///After solving, return a linear Wcs (i.e without distortion terms)
lsst::afw::image::Wcs GlobalAstrometrySolution::getWcs() throw(domain_error) {
    if (_sip == NULL) {
        throw(domain_error("Wcs not yet calculated"));
    }

    lsst::afw::image::PointD crval(_sip->wcstan.crval[0], _sip->wcstan.crval[1]);
    lsst::afw::image::PointD crpix(_sip->wcstan.crpix[1], _sip->wcstan.crpix[0]);

    int naxis = 2;   //This is hardcoded into the sip_t structure
    boost::numeric::ublas::matrix<double> CD(2,2);
    for (int i=0; i<naxis; ++i) {
        for (int j=0; j<naxis; ++j) {
            CD.insert_element(i, j, _sip->wcstan.cd[i][j]);
        }
    }

    lsst::afw::image::Wcs::Wcs wcs(crval, crpix, CD);
    return wcs;
}


/// Return the number of indices loaded from disk thus far   
int GlobalAstrometrySolution::getNumIndices() {
    return(pl_size(_backend->indexes));
}


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
        

    

//
// Mutators
//

//This function is commented out because it used the add_index() routine from astrometry.net
//but that function is not declared in a header file. Until if figure out what to do, this
//code won't compile. This is fixed in version 0.25, but 0.25 not installed yet.

///Adds a single index file to the backend
void GlobalAstrometrySolution::addIndexFile(const std::string path) { ///< Path of index file
    //Copy a constant string into a non-const C style string
    //so it can be passed into a C function without complaint.
    int len = (int) path.length(); 
    char *fn = (char *) malloc((len+1)*sizeof(char));
    strncpy(fn, path.c_str(), len);

    //May have to add some error checking. add_index always returns zero regardless
    //of success or failure
#if 0                                   // needs to be declare publicly visible before this will compile
    add_index(_backend, path);
#endif
    free(fn);
}



///Set the size of the image. At some point, when I know how to determine what my inputs are
///I may include this in the default constructor
void GlobalAstrometrySolution::setImageScaleArcsecPerPixel(double scale) { ///< Plate scale of image
    setMinimumImageScale(scale);
    setMaximumImageScale(scale);
}


///Inform the solver that the image may suffer from some distortion. Turn this on if a purely linear
///wcs solution will get the postions of some stars wrong by more than 1 pixel
void GlobalAstrometrySolution::allowDistortion(bool distort) {
    _solver->distance_from_quad_bonus = (distort) ? true : false;
}
    
        

//
// Solve and verify functions
//

///Having correctly set up your GAS object, find the ra dec of your 
///image with no prior information. Returns true if a match was found,
///false otherwise
int GlobalAstrometrySolution::blindSolve() {

    //Move all the indices from the backend structure to the solver structure.
    int size=pl_size(_backend->indexes);
    for(int i=0; i<size; ++i){
        index_t *index = (index_t *) pl_get(_backend->indexes, i);
        solver_add_index(_solver, index);
    }

    //I don't know what this setting does, but see control-program.c
    _solver->distance_from_quad_bonus = true;

    solver_run( _solver);

    return(_solver->best_match_solves);
}


bool GlobalAstrometrySolution::verifyRaDec(const double ra,   ///<Right ascension in decimal degrees
                                           const double dec   ///<Declination in decimal degrees
                                          ) throw(logic_error)  {
    //Test that a field has been assigned
    if (! _solver-> fieldxy) {
        throw logic_error("No field has been set yet");
    }

    findNearbyIndices(ra, dec);
    solver_run( _solver);

    return(_solver->best_match_solves);
}


///Verify that a previously calculated WCS (world coordinate solution) 
///matches the positions of the point sources in the image.
///Returns true if the Wcs matches
///This code is experimental, and probably won't compile
bool GlobalAstrometrySolution::verifyWcs(const lsst::afw::image::Wcs::Ptr wcsPtr) throw(logic_error) {
    //Test that a field has been assigned
    if (! _solver-> fieldxy) {
        throw logic_error("No field has been set yet");
    }


    //Extract information from the WCS object and stuff it into
    //a tan_t structure from astrometry.net
    _sip = convertWcsToSipt(wcsPtr);

     //Find all indices that are near the position indicated in the Wcs
    double ra = _sip->wcstan.crval[0];
    double dec = _sip->wcstan.crval[1];
    findNearbyIndices(ra, dec);

    solver_verify_sip_wcs(_solver, _sip);
    sip_free(_sip);

    return(_solver->best_match_solves);
}



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

    lsst::afw::image::PointD radec = (*wcsPtr).getRaDecCenter();
    sip->wcstan.crval[0] = radec.getX();
    sip->wcstan.crval[1] = radec.getY();

    lsst::afw::image::PointD rowCol = (*wcsPtr).getColRowCenter();
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
