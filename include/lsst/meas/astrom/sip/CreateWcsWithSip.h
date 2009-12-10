// -*- LSST-C++ -*-

#ifndef CREATE_WCS_WITH_SIP
#define CREATE_WCS_WITH_SIP

#include <cstdio>
#include <vector>
#include <algorithm>

#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"

#include "lsst/pex/logging/Log.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/detection/SourceMatch.h"
#include "lsst/afw/detection/Source.h"

#include "lsst/meas/astrom/sip/LeastSqFitter2d.h"


namespace except = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwImg = lsst::afw::image;
namespace det = lsst::afw::detection;


namespace lsst { 
namespace meas { 
namespace astrom { 
namespace sip {


class CreateWcsWithSip {
public:

    typedef boost::shared_ptr<CreateWcsWithSip> Ptr;
    typedef boost::shared_ptr<CreateWcsWithSip const> ConstPtr;

    CreateWcsWithSip(const std::vector<det::SourceMatch> match,
                     const afwImg::Wcs &linearWcs,
                     int order);

    CreateWcsWithSip(const std::vector<det::SourceMatch> match,
                          const afwImg::Wcs &linearWcs,
                          double maxScatterInArcsec,
                          int maxOrder);

    afwImg::Wcs getNewWcs();
    double getScatterInPixels();
    double getScatterInArcsec();
    inline int getOrder() { return  _sipA.rows(); }

private:
    
    const std::vector<det::SourceMatch> _matchList;
    const afwImg::Wcs _linearWcs;
    afwImg::Wcs _newWcs;
    Eigen::MatrixXd _sipA, _sipB;
    Eigen::MatrixXd _sipAp, _sipBp;
    
    void _createWcs(int order);

    Eigen::MatrixXd _calculateSip(const std::vector<double> u, const std::vector<double> v, 
                  const std::vector<double> z, const std::vector<double> s, int order);

    void _scaleVector(std::vector<double>& v);
    
    Eigen::MatrixXd _convertChebyToSip(Eigen::MatrixXd cheby);
                  
};    
    

}}}}

#endif
