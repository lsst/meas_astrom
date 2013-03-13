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

#include <vector>

#include "boost/shared_ptr.hpp"
#include "lsst/base.h"
#include "Eigen/Core"

#include "lsst/afw/table/Match.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/pex/logging/Log.h"

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


///
/// \brief Measure the distortions in an image plane and express them a SIP polynomials 
///
/// Given a list of matching sources between a catalogue and an image,
/// and a linear Wcs that describes the mapping from pixel space in the image
/// and ra/dec space in the catalogue, calculate discrepancies between the two
/// and compute SIP distortion polynomials to describe the discrepancy
///
/// SIP polynomials are defined in Shupe at al. (2005) ASPC 347 491.
///
/// Note that the SIP standard insists (although it is only mentioned obliquly
/// between Eqns 3 and 4) that the lowest three terms in the distortion
/// polynomials be zero (A00, A10, A01, B00, etc.). To achieve this, we need to 
/// adjust the values of CD and CRPIX from the input wcs. This may not be the 
/// behaviour you expect.
/// 
/// A Wcs may be created in a variety of ways (e.g. lsst::meas::astrom::net::GlobalAstrometrySolution ), 
/// and the list of matched sources (matches) can be generated with the matchRaDec function.
/// 
/// \code
/// #Example usage 
/// matches = matchRaDec(catSet, srcSet, 1.0*afwGeom.arcseconds, true)
/// wcs = getWcsFromSomewhere()
/// 
/// maxScatter=0.1
/// maxOrder= 10
/// sipObject = CreateWcsWithSip(matches, wcs, maxScatter, maxOrder)
/// wcs = sipObject.getNewWcs()
/// \endcode
/// 
/// Note that the matches must be one-to-one; this is ensured by passing closest=true to matchRaDec.
///
template<class MatchT>
class CreateWcsWithSip {
public:

    typedef boost::shared_ptr<CreateWcsWithSip> Ptr;
    typedef boost::shared_ptr<CreateWcsWithSip const> ConstPtr;

    CreateWcsWithSip(
        std::vector<MatchT> const & matches,
        CONST_PTR(afw::image::Wcs) linearWcs,
        int const order,
        afw::geom::Box2I const& bbox = afw::geom::Box2I(),
        int const ngrid=0
    );

    PTR(afw::image::TanWcs) getNewWcs() { return _newWcs; }

    /*
     Returns the median separation between points in this object's match list,
     projecting reference sources from RA,Dec to pixels using the SIP WCS, and
     comparing with the matched source pixel positions.
     */
    double getScatterInPixels();

    /*
     Returns the median separation between points in this object's match list,
     projecting sources from pixel space to RA,Dec using the SIP WCS, and
     comparing with the reference source RA,Dec positions.
     */
    afw::geom::Angle getScatterOnSky();

    /*
     Returns the median separation between points in this object's match list,
     projecting reference sources from RA,Dec to pixels using the input TAN
     (linear) WCS, and comparing with the matched source pixel positions.
     */
    double getLinearScatterInPixels();

    /*
     Returns the median separation between points in this object's match list,
     projecting sources from pixel space to RA,Dec using the input TAN (linear)
     WCS, and comparing with the reference source RA,Dec positions.
     */
    afw::geom::Angle getLinearScatterOnSky();

    /// Get the number of terms in the SIP matrix
    int getOrder() const { return  _sipA.rows(); }
    /// Return the number of points in the catalogue
    int getNPoints() const { return _matches.size(); }
    /// Return the number of grid points (on each axis) used in inverse SIP transform
    int getNGrid() const { return _ngrid; }

private:

    lsst::pex::logging::Log _log;
    
    std::vector<MatchT> const _matches;
    afw::geom::Box2I mutable _bbox;
    int _ngrid;                         // grid size to calculate inverse SIP coefficients (1-D)
    CONST_PTR(afw::image::Wcs) _linearWcs;
    // _sipOrder is polynomial order for forward transform.
    // _reverseSipOrder is order for reverse transform, not necessarily the same.
    int const _sipOrder, _reverseSipOrder;      

    Eigen::MatrixXd _sipA, _sipB;
    Eigen::MatrixXd _sipAp, _sipBp;

    PTR(afw::image::TanWcs) _newWcs;

    double _getScatterPixels(afw::image::Wcs const& wcs,
                             std::vector<MatchT> const & matches);
    afw::geom::Angle _getScatterSky(afw::image::Wcs const& wcs,
                                    std::vector<MatchT> const & matches);
    
    void _calculateForwardMatrices();
    void _calculateReverseMatrices();
    
    afw::geom::Point2D _getCrvalAsGeomPoint();
};    

/// Factory function for CreateWcsWithSip
template<class MatchT>
CreateWcsWithSip<MatchT> makeCreateWcsWithSip(
    std::vector<MatchT> const & matches,
    CONST_PTR(afw::image::Wcs) linearWcs,
    int const order,
    afw::geom::Box2I const& bbox = afw::geom::Box2I(),
    int const ngrid=0
) {
    return CreateWcsWithSip<MatchT>(matches, linearWcs, order, bbox, ngrid);
}

}}}}

#endif
