
#ifndef GLOBAL_ASTROMETRY_SOLUTION_H
#define GLOBAL_ASTROMETRY_SOLUTION_H

#include <iostream>
#include <vector>

#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "wcslib/wcs.h"

namespace lsst { namespace meas { namespace astrom { namespace net {

extern "C" {
#include "blind/backend.h"
#include "blind/solver.h"
#include "util/healpix.h"
#include "util/bl.h"
}



// Function declarations marked //o are not yet implemented

//!\class 
//!\brief Stores the indices of the star catalogue for the 
//!\brief astrometry.net code
///
class GlobalAstrometrySolution {
public:
    typedef boost::shared_ptr<GlobalAstrometrySolution> Ptr;
    typedef boost::shared_ptr<GlobalAstrometrySolution const> ConstPtr;

    //Constructors
    GlobalAstrometrySolution(const std::string filename); 
    GlobalAstrometrySolution(const std::string filename,  std::vector<lsst::afw::detection::Source::Ptr> src);
    
    //Destructor
    ~GlobalAstrometrySolution();

    //Initialisation routines, for those who prefer fine grained control
    int parseConfigStream(FILE* fconf);                     
    int parseConfigFile(const std::string filename);        
    void setStarlist(std::vector<lsst::afw::detection::Source::Ptr> src); //o

    //Accessors
    lsst::afw::image::Wcs getWcs();

    //The following accessors are mostly for debugging, and let you know what's going on inside
    //the object
    int getNumIndices();
    //vector<string> getIndexPaths();//o Postponed until I learn how to implement vectors in swig
    string getIndexPath(int i);

    //Mutators, mostly for tweaking parameters
    void addIndexFile(const std::string path);
    inline void setHpRange(const double range) { _hprange=range;}
    
    //Solve and verify functions.
    int blindSolve();    
    bool verifyRaDec(const double ra, const double dec);
    bool verifyWcs(const lsst::afw::image::Wcs::Ptr wcsPtr);
    
private:
    backend_t *_backend;
    solver_t *_solver;
    starxy_t *_starlist;
    sip_t *_sip;
    
    double _hprange;
    
    sip_t *convertWcsToSipt(const lsst::afw::image::Wcs::Ptr);
    void findNearbyIndices(const double ra, const double dec);
};

}}}}

#endif


