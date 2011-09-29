// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 

#ifndef CREATE_WCS_WITH_SIP
#define CREATE_WCS_WITH_SIP

#include <cstdio>
#include <vector>
#include <algorithm>

#include "boost/shared_ptr.hpp"
#include "lsst/base.h"
#include "Eigen/Core.h"
#include "Eigen/SVD"
#include "Eigen/Cholesky"
#include "Eigen/LU"


#include "lsst/pex/logging/Log.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/detection/SourceMatch.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"



namespace lsst { 
    namespace afw {
        namespace image {
            class Wcs;
            class TanWcs;
        }
    }
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

    CreateWcsWithSip(const std::vector<lsst::afw::detection::SourceMatch> match,
                     const lsst::afw::image::Wcs::Ptr linearWcs,
                     int order);


    PTR(lsst::afw::image::TanWcs) getNewWcs();
    double getScatterInPixels();
    double getScatterInArcsec();

    ///Get the number of terms in the SIP matrix
    inline int getOrder() { return  _sipA.rows(); }

    inline int getNPoints() { return _size; }

private:
    
    const std::vector<lsst::afw::detection::SourceMatch> _matchList;
    CONST_PTR(lsst::afw::image::Wcs) _linearWcs;
    CONST_PTR(lsst::afw::image::Wcs) _origWcs;
    //size is number of input points. _sipOrder is polynomial order for forward transform.
    //_reverseSipOrder is order for reverse transform, not necessarily the same.
    const int _sipOrder, _reverseSipOrder, _size;      

    Eigen::MatrixXd _sipA, _sipB;
    Eigen::MatrixXd _sipAp, _sipBp;

    PTR(lsst::afw::image::TanWcs) _newWcs;    
    
    void _calculateForwardMatrices();
    void _calculateReverseMatrices();
    
    Eigen::MatrixXd _calculateCMatrix(Eigen::VectorXd axis1, Eigen::VectorXd axis2, int order);
    Eigen::VectorXd _leastSquaresSolve(Eigen::VectorXd b, Eigen::MatrixXd A);
    
    lsst::afw::geom::Point2D _getCrvalAsGeomPoint();
    
    int getUIndex(int j, int order);
    int getVIndex(int j, int order);
    
    void debug();
    
};    
    

}}}}

#endif
