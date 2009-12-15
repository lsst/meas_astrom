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


///\brief Measure the distortions in an image plane and express them a SIP polynomials 
///
/// Given a list of matching sources between a catalogue and an image,
/// and a linear Wcs that describes the mapping from pixel space in the image
/// and ra/dec space in the catalogue, calculate discrepancies between the two
/// and compute SIP distortion polynomials to describe the discrepancy
///
///SIP polynomials are defined in Shupe at al. (2005) ASPC 347 491.
///
/// Note that the SIP standard insists (although it is only mentioned obliquly
/// between Eqns 3 and 4) that the lowest three terms in the distortion
/// polynomials be zero (A00, A10, A01, B00, etc.). To achieve this, we need to 
/// adjust the values of CD and CRPIX from the input wcs. This may not be the 
/// behaviour you expect.
/// 
/// A Wcs may be created in a variety of ways (e.g. lsst::meas::astrom::net::GlobalAstrometrySolution ), 
/// and the
/// list of matched sources (match) can be generated with MatchSrcToCatalogue)
/// 
/// \code
/// #Example usage 
/// match = MatchSrcToCatalogue(catSet, srcSet)
/// wcs = getWcsFromSomewhere()
/// 
/// maxScatter=0.1
/// maxOrder= 10
/// sipObject = CreateWcsWithSip(match, wcs, maxScatter, maxOrder)
/// wcs = sipObject.getNewWcs()
/// \endcode
/// 
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
    ///Get the number of terms in the SIP matrix
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
