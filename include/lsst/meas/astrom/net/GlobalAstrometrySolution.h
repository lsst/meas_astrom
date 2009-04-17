

#ifndef GLOBAL_ASTROMETRY_SOLUTION_H
#define GLOBAL_ASTROMETRY_SOLUTION_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/utils/Utils.h"

#include "wcslib/wcs.h"

namespace lsst { namespace meas { namespace astrom { namespace net {

extern "C" {
#include "backend.h"
#include "solver.h"
#include "index.h"
#include "tweak.h"
#include "healpix.h"
#include "bl.h"
#include "log.h"
}


///Define variables to indicate tha parity of the image (i.e whether, when the image is rotated so that North is
///up, East is to the left (Normal), or to the right, (FLIPPED)). The constants PARITY_* are defined
///in solver.h, but not visible to Python. The default setting is UNKNOWN_PARITY. Setting the parity correctly
///reduces the time taken to solve an image by about a factor of two.
enum {
    NORMAL_PARITY = PARITY_NORMAL,
    FLIPPED_PARITY = PARITY_FLIP,
    UNKNOWN_PARITY = PARITY_BOTH
};
    


//!\brief Solve for WCS based only on the positions of stars in an image 
///
///See the examples/ directory for an example of how to use
class GlobalAstrometrySolution {
public:
    typedef boost::shared_ptr<GlobalAstrometrySolution> Ptr;
    typedef boost::shared_ptr<GlobalAstrometrySolution const> ConstPtr;

    //Constructors
    GlobalAstrometrySolution(const std::string policyPath);
    
    //Destructor
    ~GlobalAstrometrySolution();

    //Initialisation routines, for those who prefer fine grained control.
    void addIndexFile(const std::string path);        
    int parseConfigFile(const std::string filename);        
    int parseConfigStream(FILE* fconf);                     
    void setStarlist(lsst::afw::detection::SourceSet vec);


    //Accessors
    double getMatchThreshold();
    inline double getMinimumImageScale() {    return _solver->funits_lower; }
    inline double getMaximumImageScale() {    return _solver->funits_upper; }
    inline double getMinQuadScale(){    return _solver->quadsize_min;}
    bool isFlipped(); 
    
    double getSolvedImageScale();
    lsst::afw::image::Wcs::Ptr getDistortedWcs(int order=3);
    lsst::afw::image::Wcs::Ptr getWcs();
    lsst::afw::image::PointD raDecToXY(double ra, double dec);
    lsst::afw::image::PointD xyToRaDec(double x, double y);

    //The following accessors are mostly for debugging, and let you know what's going on inside
    //the object
    int getNumIndices();
    std::vector<std::string> getIndexPaths();
    void printStarlist();

    //Mutators, mostly for tweaking parameters
    void allowDistortion(bool distort);
    void reset();
    void setDefaultValues();
    void setImageScaleArcsecPerPixel(double scale);
    void setLogLevel(const int level);
    void setMatchThreshold(const double threshold);
    ///Note than minimum image scale should be strictly less than Maximum scale
    inline void setMinimumImageScale(double scale){   _solver->funits_lower=scale;}
    inline void setMaximumImageScale(double scale){   _solver->funits_upper=scale;}
    
    ///Set the maximum number of stars to pass to the solver. Fewer stars means a faster
    ///solution, more stars means a greater chance of finding a match. In practice, 50
    ///is a good number to choose
    inline void setNumberStars(const int num)  {    _solver->endobj = num;}

    ///Set the scale (in pixels) of the smallest quad (group of 4 stars) to match
    ///against the database. You don't ususually need to use this function, but it
    ///can be useful in debugging
    inline void setMinQuadScale(const double scale){    _solver->quadsize_min = scale;}
    void setParity(const int parity);
    
    //Solve and verify functions.
    bool solve();
    bool solve(const afw::image::PointD raDec);
    bool solve(const double ra, const double dec);
    bool solve(const lsst::afw::image::Wcs::Ptr wcsPtr, const double imageScaleUncertaintyPercent=20);
    lsst::afw::image::Wcs::Ptr solve(const lsst::afw::detection::SourceSet vec, const lsst::afw::image::Wcs::Ptr wcsPtr);
    
    //Not implemented yet
    #if 0
    bool verifyWcs(const lsst::afw::image::Wcs::Ptr wcsPtr);
    #endif
    
private:
    backend_t *_backend;
    solver_t *_solver;
    starxy_t *_starxy;   ///List of sources to solve for
    int _numBrightObjects;   //Only use the brightest objects in solve.
    
    //Variables indicating the coordinate system of the solution
    double _equinox;
    std::string _raDecSys;
    
    sip_t *convertWcsToSipt(const lsst::afw::image::Wcs::Ptr);
    void loadNearbyIndices(std::vector<double> unitVector);
    void solverSetField();
};

}}}}

#endif


