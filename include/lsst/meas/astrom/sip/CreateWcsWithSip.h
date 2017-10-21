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

#include <memory>
#include <vector>

#include "lsst/base.h"
#include "Eigen/Core"

#include "lsst/afw/table/Match.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Angle.h"

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace sip {


/**
 \brief Measure the distortions in an image plane and express them a SIP polynomials 

 Given a list of matching sources between a catalogue and an image,
 and a linear Wcs that describes the mapping from pixel space in the image
 and ra/dec space in the catalogue, calculate discrepancies between the two
 and compute SIP distortion polynomials to describe the discrepancy

 SIP polynomials are defined in Shupe at al. (2005) ASPC 347 491.

 Note that the SIP standard insists (although it is only mentioned obliquly
 between Eqns 3 and 4) that the lowest three terms in the distortion
 polynomials be zero (A00, A10, A01, B00, etc.). To achieve this, we need to 
 adjust the values of CD and CRPIX from the input wcs. This may not be the 
 behaviour you expect.

 A Wcs may be created in a variety of ways (e.g. lsst::meas::astrom::net::GlobalAstrometrySolution ), 
 and the list of matched sources (matches) can be generated with the matchRaDec function.

 \code
 #Example usage 
 matches = matchRaDec(catSet, srcSet, 1.0*afwGeom.arcseconds, true)
 wcs = getWcsFromSomewhere()

 maxScatter=0.1
 maxOrder= 10
 sipObject = CreateWcsWithSip(matches, wcs, maxScatter, maxOrder)
 wcs = sipObject.getNewWcs()
 \endcode

 Note that the matches must be one-to-one; this is ensured by passing closest=true to matchRaDec.
 */
template<class MatchT>
class CreateWcsWithSip {
public:

    typedef std::shared_ptr<CreateWcsWithSip> Ptr;
    typedef std::shared_ptr<CreateWcsWithSip const> ConstPtr;

    /**
     Construct a CreateWcsWithSip

     \param[in] matches  list of matches
     \param[in] linearWcs  initial WCS, typically pure TAN but need not be
     \param[in] order  SIP order for fit WCS
     \param[in] bbox  bounding box over which to compute the reverse SIP transform.
                         If empty then a bounding box is computed based on the matches,
                         extended a bit to allow for the fact that the sources
                         will not necessarily reach to each edge of the image.
                         Specifially the box is grown by dimensions/sqrt(number of matches).
     \param[in] ngrid  number of points along x or y for the grid of points on which
                         the reverse SIP transform is computed
     */
    CreateWcsWithSip(
        std::vector<MatchT> const & matches,
        afw::geom::SkyWcs const & linearWcs,
        int const order,
        afw::geom::Box2I const& bbox = afw::geom::Box2I(),
        int const ngrid=0
    );

    std::shared_ptr<afw::geom::SkyWcs> getNewWcs() { return _newWcs; }

    /**
     Compute the median separation, in pixels, between items in this object's match list

     For each match, project the reference object coord to pixels using the fit TAN-SIP WCS,
     and measure the radial separation to the source centroid
     */
    double getScatterInPixels() const;

    /**
     Compute the median on-sky separation between items in this object's match list

     For each match, project the source centroid to RA,Dec using the fit TAN-SIP WCS,
     and measure the on-sky angular separation to the reference source coord.
     */
    afw::geom::Angle getScatterOnSky() const;

    /**
     Compute the median radial separation between items in this object's match list

     For each match, project the reference object coord to pixels using the initial "linearWcs" WCS,
     and measure the radial separation to the source centroid.
     */
    double getLinearScatterInPixels() const;

    /**
     Compute the median on-sky separation between items in this object's match list,

     For each match, project the source centroid to RA,Dec using the initial "linearWcs" WCS,
     and measure the on-sky angular separation to the reference source coord.
     */
    afw::geom::Angle getLinearScatterOnSky() const;

    /// Return the number of terms in the SIP matrix
    int getOrder() const { return  _sipA.rows(); }
    /// Return the number of points in the catalogue
    int getNPoints() const { return _matches.size(); }
    /// Return the number of grid points (on each axis) used in inverse SIP transform
    int getNGrid() const { return _ngrid; }

    // Return the SIP A matrix
    Eigen::MatrixXd const getSipA() { return _sipA; }
    // Return the SIP B matrix
    Eigen::MatrixXd const getSipB() { return _sipB; }
    // Return the SIP Ap matrix
    Eigen::MatrixXd const getSipAp() { return _sipAp; }
    // Return the SIP Bp matrix
    Eigen::MatrixXd const getSipBp() { return _sipBp; }

private:

    std::vector<MatchT> const _matches;
    afw::geom::Box2I mutable _bbox;
    int _ngrid;                         // grid size to calculate inverse SIP coefficients (1-D)
    std::shared_ptr<const afw::geom::SkyWcs> _linearWcs;
    // _sipOrder is polynomial order for forward transform.
    // _reverseSipOrder is order for reverse transform, not necessarily the same.
    int const _sipOrder, _reverseSipOrder;      

    Eigen::MatrixXd _sipA, _sipB;
    Eigen::MatrixXd _sipAp, _sipBp;

    std::shared_ptr<afw::geom::SkyWcs> _newWcs;

    void _calculateForwardMatrices();
    void _calculateReverseMatrices();
};    

/// Factory function for CreateWcsWithSip
template<class MatchT>
CreateWcsWithSip<MatchT> makeCreateWcsWithSip(
    std::vector<MatchT> const & matches,
    afw::geom::SkyWcs const& linearWcs,
    int const order,
    afw::geom::Box2I const& bbox = afw::geom::Box2I(),
    int const ngrid=0
) {
    return CreateWcsWithSip<MatchT>(matches, linearWcs, order, bbox, ngrid);
}

}}}}

#endif
