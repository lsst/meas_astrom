

#ifndef GLOBAL_ASTROMETRY_SOLUTION_H
#define GLOBAL_ASTROMETRY_SOLUTION_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/utils/Utils.h"

#include "wcslib/wcs.h"

namespace lsst { namespace meas { namespace astrom { namespace net {

using namespace std;

extern "C" {
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

    //Tune the solver
    void setDefaultValues();
    void setStarlist(lsst::afw::detection::SourceSet vec);
    void setNumBrightObjects(const int N);
    inline void setMinimumImageScale(double scale){   _solver->funits_lower=scale;}
    inline void setMaximumImageScale(double scale){   _solver->funits_upper=scale;}
    void setImageScaleArcsecPerPixel(double scale);
    void allowDistortion(bool distort);
    void setLogLevel(const int level);
    void setMatchThreshold(const double threshold);
    void setParity(const int parity);

    //Solve for a wcs solution
    bool solve();
    bool solve(const afw::image::PointD raDec);
    bool solve(const double ra, const double dec);
    bool solve(const lsst::afw::image::Wcs::Ptr wcsPtr, const double imageScaleUncertaintyPercent=20);

    //Return the solution
    lsst::afw::image::Wcs::Ptr getWcs();
    lsst::afw::image::Wcs::Ptr getDistortedWcs(int order=3);
    lsst::afw::detection::SourceSet getMatchedSources();
    double getSolvedImageScale();

    //Call this before performing a new match
    void reset();



private:
    pl *_indexList;
    pl *_metaList;
    solver_t *_solver;
    starxy_t *_starxy;   ///List of sources to solve for
    int _numBrightObjects;   //Only use the brightest objects in solve.
    
    //Variables indicating the coordinate system of the solution
    double _equinox;
    std::string _raDecSys;

    index_meta_t *_loadIndexMeta(string filename);

    void _solverSetField();
    bool _isIndexMetaPossibleMatch(index_meta_t *meta, double ra, double dec);
    bool _isMetaNearby(index_meta_t *meta, double ra, double dec); 
    bool _isMetaSuitableScale(index_meta_t *meta);                               

};

}}}}
#endif
