// -*- LSST-C++ -*-

#ifndef GLOBAL_ASTROMETRY_SOLUTION_H
#define GLOBAL_ASTROMETRY_SOLUTION_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/format.hpp"
#include "Eigen/Core.h"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/utils/Utils.h"

#include "wcslib/wcs.h"

extern "C" {
#include "solver.h"
#include "index.h"
#include "tweak.h"
#include "healpix.h"
#include "bl.h"
#include "log.h"
}

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace net {


///Define variables to indicate tha parity of the image (i.e whether, when the image is rotated so that 
///North is up, East is to the left (Normal), or to the right, (FLIPPED)). The constants PARITY_* are defined
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
    explicit GlobalAstrometrySolution(const std::string policyPath);
    
    //Destructor
    ~GlobalAstrometrySolution();

    //Tune the solver
    void setDefaultValues();
    void setStarlist(lsst::afw::detection::SourceSet vec);
    void setNumBrightObjects(int N);
    inline void setMinimumImageScale(double scale){   _solver->funits_lower = scale; }
    inline void setMaximumImageScale(double scale){   _solver->funits_upper = scale; }
    void setImageScaleArcsecPerPixel(double scale);
    void allowDistortion(bool hasDistortion);
    void setLogLevel(int level);
    void setMatchThreshold(double threshold);
    void setParity(int parity);

    //Solve for a wcs solution
    bool solve();
    bool solve(const afw::image::PointD raDec);
    bool solve(double ra, double dec);
    bool solve(const lsst::afw::image::Wcs::Ptr wcsPtr, double imageScaleUncertaintyPercent = 20);

    //Return the solution
    lsst::afw::image::Wcs::Ptr getWcs();
    lsst::afw::image::Wcs::Ptr getDistortedWcs(int order = 3);
    lsst::afw::detection::SourceSet getMatchedSources();
    double getSolvedImageScale();
    lsst::afw::detection::SourceSet getCatalogue(double ra, double dec, double radiusInArcsec);    
    lsst::afw::detection::SourceSet getCatalogue(double radiusInArcsec);

    //Call this before performing a new match
    void reset();



private:
    lsst::pex::logging::Log _mylog;

    std::vector<index_t*> _indexList;

    solver_t *_solver;
    starxy_t *_starxy;   ///List of sources to solve for
    int _numBrightObjects;   //Only use the brightest objects in solve.
        
    //Variables indicating the coordinate system of the solution
    double _equinox;
    std::string _raDecSys;

    index_t *_loadIndexMeta(std::string filename);

    void _solverSetField();
    int _addSuitableIndicesToSolver(double minImageSizeArcsec, double maxImageSizeArcsec, \
        double ra=-360, double dec=-360);    

};

}}}}
#endif
