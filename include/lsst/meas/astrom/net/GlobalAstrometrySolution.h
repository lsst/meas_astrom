

#ifndef GLOBAL_ASTROMETRY_SOLUTION_H
#define GLOBAL_ASTROMETRY_SOLUTION_H

#include <iostream>
#include <vector>
#include <string>

#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "wcslib/wcs.h"

namespace lsst { namespace meas { namespace astrom { namespace net {

extern "C" {
#include "backend.h"
#include "solver.h"
#include "index.h"
#include "healpix.h"
#include "bl.h"
#include "log.h"
}




//!\brief Solve for WCS based only on the positions of stars in an image 
///
///See the tests/ directory for an example of how to use
class GlobalAstrometrySolution {
public:
    typedef boost::shared_ptr<GlobalAstrometrySolution> Ptr;
    typedef boost::shared_ptr<GlobalAstrometrySolution const> ConstPtr;

    //Constructors
    GlobalAstrometrySolution();
    GlobalAstrometrySolution(lsst::afw::detection::SourceVector vec);

    
    //Destructor
    ~GlobalAstrometrySolution();

    //Initialisation routines, for those who prefer fine grained control.
    void addIndexFile(const std::string path);        
    int parseConfigFile(const std::string filename);        
    int parseConfigStream(FILE* fconf);                     
    void setStarlist(lsst::afw::detection::SourceVector vec) throw (std::domain_error);


    //Accessors
    double getMatchThreshold();
    inline double getMinimumImageScale() {    return _solver->funits_lower; }
    inline double getMaximumImageScale() {    return _solver->funits_upper; }
    lsst::afw::image::Wcs::Ptr getWcs() throw(std::logic_error);
    std::pair<double, double> raDec2Xy(double ra, double dec) throw(std::logic_error);
    std::pair<double, double> xy2RaDec(double x, double y) throw(std::logic_error);

    //The following accessors are mostly for debugging, and let you know what's going on inside
    //the object
    int getNumIndices();
    std::vector<std::string> getIndexPaths();
    void printStarlist();

    //Mutators, mostly for tweaking parameters
    void allowDistortion(bool distort);
    void reset();
    void setDefaultValues();
    void setHpRange(const double range) { _hprange=range;}
    void setImageScaleArcsecPerPixel(double scale);
    void setLogLevel(const int level);
    void setMatchThreshold(const double threshold);
    ///Note than minimum image scale should be strictly less than Maximum scale
    inline void setMinimumImageScale(double scale){   _solver->funits_lower=scale;}
    inline void setMaximumImageScale(double scale){   _solver->funits_upper=scale;}
    inline void setNumberStars(const int num)  { _solver->endobj = num;}
    
    //Solve and verify functions.
    int blindSolve();    
    bool verifyRaDec(const double ra, const double dec) throw(std::logic_error);
    bool verifyWcs(const lsst::afw::image::Wcs::Ptr wcsPtr) throw(std::logic_error);
    
private:
    backend_t *_backend;
    solver_t *_solver;
    starxy_t *_starlist;
    
    double _hprange;
    
    sip_t *convertWcsToSipt(const lsst::afw::image::Wcs::Ptr) throw(std::logic_error);
    void findNearbyIndices(const double ra, const double dec);
};

}}}}

#endif


