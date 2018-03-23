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
#include <cmath> 

#include "Eigen/SVD"
#include "Eigen/Cholesky"
#include "Eigen/LU"

#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/meas/astrom/sip/CreateWcsWithSip.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/log/Log.h"
#include "lsst/meas/astrom/makeMatchStatistics.h"

namespace {
LOG_LOGGER _log = LOG_GET("meas.astrom.sip");
}

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace sip {

namespace except   = lsst::pex::exceptions;
namespace afwGeom  = lsst::afw::geom;
namespace afwImg   = lsst::afw::image;
namespace afwDet   = lsst::afw::detection;
namespace afwMath  = lsst::afw::math;
namespace afwTable = lsst::afw::table;

namespace {

double const MAX_DISTANCE_CRPIX_TO_BBOXCTR = 1000;

/*
 * Given an index and a SIP order, calculate p and q for the index'th term u^p v^q
 * (Cf. Eqn 2 in http://fits.gsfc.nasa.gov/registry/sip/SIP_distortion_v1_0.pdf)
 */
std::pair<int, int>
indexToPQ(int const index, int const order)
{
    int p = 0, q = index;
    
    for (int decrement = order; q >= decrement && decrement > 0; --decrement) {
        q -= decrement;
        p++;
    }
    
    return std::make_pair(p, q);
}

Eigen::MatrixXd
calculateCMatrix(Eigen::VectorXd const& axis1, Eigen::VectorXd const& axis2, int const order)
{
    int nTerms = 0.5*order*(order+1);

    int const n = axis1.size();
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n, nTerms);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nTerms; ++j) {
            std::pair<int, int> pq = indexToPQ(j, order);
            int p = pq.first, q = pq.second;

            assert(p + q < order);
            C(i, j) = ::pow(axis1[i], p)*::pow(axis2[i], q);
        }
    }
    
    return C;
}
    
///Given a vector b and a matrix A, solve b - Ax = 0
/// b is an m x 1 vector, A is an n x m matrix, and x, the output is a 
///\param b An m x 1 vector, where m is the number of parameters in the fit
///\param A An n x m vecotr, where n is the number of equations in the solution
///
///\returns x, an m x 1 vector of best fit params
Eigen::VectorXd
leastSquaresSolve(Eigen::VectorXd b, Eigen::MatrixXd A) {
    assert(A.rows() == b.rows());
    Eigen::VectorXd par = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    return par;
}

} // anonymous namespace

///Constructor
template<class MatchT>
CreateWcsWithSip<MatchT>::CreateWcsWithSip(
    std::vector<MatchT> const & matches,
    afw::geom::SkyWcs const& linearWcs,
    int const order,
    afwGeom::Box2I const& bbox,
    int const ngrid
):
    _matches(matches),
    _bbox(bbox),
    _ngrid(ngrid),
    _linearWcs(std::make_shared<afw::geom::SkyWcs>(linearWcs)),
    _sipOrder(order+1),
    _reverseSipOrder(order+2), //Higher order for reverse transform
    _sipA(Eigen::MatrixXd::Zero(_sipOrder, _sipOrder)),
    _sipB(Eigen::MatrixXd::Zero(_sipOrder, _sipOrder)),
    _sipAp(Eigen::MatrixXd::Zero(_reverseSipOrder, _reverseSipOrder)),
    _sipBp(Eigen::MatrixXd::Zero(_reverseSipOrder, _reverseSipOrder)),
    _newWcs()
{
    if  (order < 2) {
        throw LSST_EXCEPT(except::OutOfRangeError, "SIP must be at least 2nd order");        
    }
    if (_sipOrder > 9) {
        throw LSST_EXCEPT(except::OutOfRangeError,
                          str(boost::format("SIP forward order %d exceeds the convention limit of 9") %
                              _sipOrder));
    }
    if (_reverseSipOrder > 9) {
        throw LSST_EXCEPT(except::OutOfRangeError,
                          str(boost::format("SIP reverse order %d exceeds the convention limit of 9") %
                              _reverseSipOrder));
    }
    
    if (_matches.size() < std::size_t(_sipOrder)) {
        throw LSST_EXCEPT(except::LengthError, "Number of matches less than requested sip order");
    }

    if (_ngrid <= 0) {
        _ngrid = 5*_sipOrder;           // should be plenty
    }

    /*
     * We need a bounding box to define the region over which:
     *    The forward transformation should be valid
     *    We calculate the reverse transformartion
     * If no BBox is provided, guess one from the input points (extrapolated a bit to allow for fact
     * that a finite number of points won't reach to the edge of the image)
     */
    if (_bbox.isEmpty() && !_matches.empty() > 0) {
        for (
            typename std::vector<MatchT>::const_iterator ptr = _matches.begin();
            ptr != _matches.end();
            ++ptr
        ) {
            afwTable::SourceRecord const & src = *ptr->second;
            _bbox.include(afwGeom::PointI(src.getX(), src.getY()));
        }
        float const borderFrac = 1/::sqrt(_matches.size()); // fractional border to add to exact BBox
        afwGeom::Extent2I border(borderFrac*_bbox.getWidth(), borderFrac*_bbox.getHeight());

        _bbox.grow(border);
    }

    // If crpix is too far from the center of the fit bbox, move it to the center to improve the fit
    auto const initialCrpix = _linearWcs->getPixelOrigin();
    auto const bboxCenter = afw::geom::Box2D(_bbox).getCenter();
    if (std::hypot(initialCrpix[0] - bboxCenter[0], initialCrpix[1] - bboxCenter[1]) >
        MAX_DISTANCE_CRPIX_TO_BBOXCTR) {
        LOGL_DEBUG(_log,
                   "_linearWcs crpix = %d, %d farther than %f from bbox center; shifting to center %d, %d",
                   initialCrpix[0], initialCrpix[1], MAX_DISTANCE_CRPIX_TO_BBOXCTR, bboxCenter[0],
                   bboxCenter[1]);
        _linearWcs = _linearWcs->getTanWcs(bboxCenter);
    }

    // Calculate the forward part of the SIP distortion
    _calculateForwardMatrices();

    // Build a new wcs incorporating the forward SIP matrix, it's all we know so far
    auto const crval = _linearWcs->getSkyOrigin();
    auto const crpix = _linearWcs->getPixelOrigin();
    Eigen::MatrixXd cdMatrix = _linearWcs->getCdMatrix();
    _newWcs = afw::geom::makeTanSipWcs(crpix, crval, cdMatrix, _sipA, _sipB);

    // Use _newWcs to calculate the forward transformation on a grid, and derive the back transformation
    _calculateReverseMatrices();

    // Build a new wcs incorporating all SIP matrices
    _newWcs = afw::geom::makeTanSipWcs(crpix, crval, cdMatrix, _sipA, _sipB, _sipAp, _sipBp);
}

template<class MatchT>
void
CreateWcsWithSip<MatchT>::_calculateForwardMatrices()
{
    // Assumes FITS (1-indexed) coordinates.
    afwGeom::Point2D crpix = _linearWcs->getPixelOrigin();

    // Calculate u, v and intermediate world coordinates
    int const nPoints = _matches.size();
    Eigen::VectorXd u(nPoints), v(nPoints), iwc1(nPoints), iwc2(nPoints);

    int i = 0;
    auto linearIwcToSky = getIntermediateWorldCoordsToSky(*_linearWcs);
    for (
        typename std::vector<MatchT>::const_iterator ptr = _matches.begin();
        ptr != _matches.end();
        ++ptr, ++i
    ) {
        afwTable::ReferenceMatch const & match = *ptr;

        // iwc: intermediate world coordinate positions of catalogue objects
        auto c = match.first->getCoord();
        afwGeom::Point2D p = linearIwcToSky->applyInverse(c);
        iwc1[i] = p[0];
        iwc2[i] = p[1];
        // u and v are intermediate pixel coordinates of observed (distorted) positions
        u[i] = match.second->getX() - crpix[0];
        v[i] = match.second->getY() - crpix[1];
    }
    // Scale u and v down to [-1,,+1] in order to avoid too large numbers in the polynomials
    double uMax = u.cwiseAbs().maxCoeff();
    double vMax = v.cwiseAbs().maxCoeff();
    double norm = std::max(uMax, vMax);
    u = u/norm;
    v = v/norm;
   
    // Forward transform
    int ord = _sipOrder;
    Eigen::MatrixXd forwardC = calculateCMatrix(u, v, ord);
    Eigen::VectorXd mu = leastSquaresSolve(iwc1, forwardC);
    Eigen::VectorXd nu = leastSquaresSolve(iwc2, forwardC);

    // Use mu and nu to refine CD

    // Given the implementation of indexToPQ(), the refined values
    // of the elements of the CD matrices are in elements 1 and "_sipOrder" of mu and nu
    // If the implementation of indexToPQ() changes, these assertions
    // will catch that change.
    assert ((indexToPQ(0,   ord) == std::pair<int, int>(0, 0)));
    assert ((indexToPQ(1,   ord) == std::pair<int, int>(0, 1)));
    assert ((indexToPQ(ord, ord) == std::pair<int, int>(1, 0)));

    // Scale back CD matrix
    Eigen::Matrix2d CD;
    CD(1,0) = nu[ord]/norm;
    CD(1,1) = nu[1]/norm;
    CD(0,0) = mu[ord]/norm;
    CD(0,1) = mu[1]/norm;

    Eigen::Matrix2d CDinv = CD.inverse();   //Direct inverse OK for 2x2 matrix in Eigen

    // The zeroth elements correspond to a shift in crpix
    crpix[0] -= mu[0]*CDinv(0,0) + nu[0]*CDinv(0,1); 
    crpix[1] -= mu[0]*CDinv(1,0) + nu[0]*CDinv(1,1);

    auto const crval = _linearWcs->getSkyOrigin();

    _linearWcs = afw::geom::makeSkyWcs(crpix, crval, CD);

    //Get Sip terms
    
    //The rest of the elements correspond to
    //mu[i] == CD11*Apq + CD12*Bpq and
    //nu[i] == CD21*Apq + CD22*Bpq and
    //
    //We solve for Apq and Bpq with the equation
    // (Apq)  = (CD11 CD12)-1  * (mu[i])  
    // (Bpq)    (CD21 CD22)      (nu[i])
    
    for(int i=1; i< mu.rows(); ++i) {
        std::pair<int, int> pq = indexToPQ(i, ord);
        int p = pq.first, q = pq.second;

        if(p + q > 1 && p + q < ord) {    
            Eigen::Vector2d munu(2,1);
            munu(0) = mu(i);
            munu(1) = nu(i);
            Eigen::Vector2d AB = CDinv*munu;
            // Scale back sip coefficients
            _sipA(p,q) = AB[0]/::pow(norm,p+q);
            _sipB(p,q) = AB[1]/::pow(norm,p+q);
        }
    }
}

template<class MatchT>
void CreateWcsWithSip<MatchT>::_calculateReverseMatrices() {
    int const ngrid2 = _ngrid*_ngrid;

    Eigen::VectorXd U(ngrid2), V(ngrid2);
    Eigen::VectorXd delta1(ngrid2), delta2(ngrid2);
    
    int const x0 = _bbox.getMinX();
    double const dx = _bbox.getWidth()/(double)(_ngrid - 1);
    int const y0 = _bbox.getMinY();
    double const dy = _bbox.getHeight()/(double)(_ngrid - 1);

    // wcs->getPixelOrigin() returns LSST-style (0-indexed) pixel coords.
    afwGeom::Point2D crpix = _newWcs->getPixelOrigin();

    LOGL_DEBUG(_log, "_calcReverseMatrices: x0,y0 %i,%i, W,H %i,%i, ngrid %i, dx,dy %g,%g, CRPIX %g,%g",
               x0, y0, _bbox.getWidth(), _bbox.getHeight(), _ngrid, dx, dy, crpix[0], crpix[1]);

    auto tanWcs = _newWcs->getTanWcs(_newWcs->getPixelOrigin());
    auto applySipAB = afwGeom::makeWcsPairTransform(*_newWcs, *tanWcs);
    int k = 0;
    for (int i = 0; i < _ngrid; ++i) {
        double const y = y0 + i*dy;
        std::vector<afwGeom::Point2D> beforeSipABPoints;
        beforeSipABPoints.reserve(_ngrid);
        for (int j = 0; j < _ngrid; ++j) {
            double const x = x0 + j*dx;
            beforeSipABPoints.emplace_back(afwGeom::Point2D(x, y));
        }
        auto const afterSipABPoints = applySipAB->applyForward(beforeSipABPoints);
        for (int j = 0; j < _ngrid; ++j, ++k) {
            double const x = x0 + j*dx;
            double u,v;
            // u and v are intermediate pixel coordinates on a grid of positions
            u = x - crpix[0];
            v = y - crpix[1];

            // U and V are the result of applying the "forward" (A,B) SIP coefficients
            afwGeom::Point2D const xy = afterSipABPoints[j];
            U[k] = xy[0] - crpix[0];
            V[k] = xy[1] - crpix[1];

            if ((i == 0 || i == (_ngrid-1) || i == (_ngrid/2)) &&
                (j == 0 || j == (_ngrid-1) || j == (_ngrid/2))) {
                LOGL_DEBUG(_log, "  x,y (%.1f, %.1f), u,v (%.1f, %.1f), U,V (%.1f, %.1f)", x, y, u, v, U[k], V[k]);
            }

            delta1[k] = u - U[k];
            delta2[k] = v - V[k];
        }
    }
    
    // Scale down U and V in order to avoid too large numbers in the polynomials
    double UMax = U.cwiseAbs().maxCoeff();
    double VMax = V.cwiseAbs().maxCoeff();
    double norm = (UMax > VMax) ? UMax : VMax;
    U = U/norm;
    V = V/norm;

    // Reverse transform
    int const ord = _reverseSipOrder;
    Eigen::MatrixXd reverseC = calculateCMatrix(U, V, ord); 
    Eigen::VectorXd tmpA = leastSquaresSolve(delta1, reverseC);
    Eigen::VectorXd tmpB = leastSquaresSolve(delta2, reverseC);
    
    assert(tmpA.rows() == tmpB.rows());
    for(int j=0; j< tmpA.rows(); ++j) {
        std::pair<int, int> pq = indexToPQ(j, ord);
        int p = pq.first, q = pq.second;
        // Scale back sip coefficients
        _sipAp(p, q) = tmpA[j]/::pow(norm,p+q);
        _sipBp(p, q) = tmpB[j]/::pow(norm,p+q);   
    } 
}

template<class MatchT>
double CreateWcsWithSip<MatchT>::getScatterInPixels() const {
    assert(_newWcs.get());
    return makeMatchStatisticsInPixels(*_newWcs, _matches, afw::math::MEDIAN).getValue();
}

template<class MatchT>
double CreateWcsWithSip<MatchT>::getLinearScatterInPixels() const {
    assert(_linearWcs.get());
    return makeMatchStatisticsInPixels(*_linearWcs, _matches, afw::math::MEDIAN).getValue();
}

template<class MatchT>
afwGeom::Angle CreateWcsWithSip<MatchT>::getScatterOnSky() const {
    assert(_newWcs.get());
    return makeMatchStatisticsInRadians(
        *_newWcs, _matches, afw::math::MEDIAN).getValue()*afw::geom::radians;
}

template<class MatchT>
afwGeom::Angle CreateWcsWithSip<MatchT>::getLinearScatterOnSky() const {
    assert(_linearWcs.get());
    return makeMatchStatisticsInRadians(
        *_linearWcs, _matches, afw::math::MEDIAN).getValue()*afw::geom::radians;
}


#define INSTANTIATE(MATCH) \
    template class CreateWcsWithSip<MATCH>;

INSTANTIATE(afwTable::ReferenceMatch);
INSTANTIATE(afwTable::SourceMatch);

}}}}


