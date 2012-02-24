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
 

#include "Eigen/SVD"
#include "Eigen/Cholesky"
#include "Eigen/LU"

#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/meas/astrom/sip/CreateWcsWithSip.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/image/TanWcs.h"

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace sip {

using namespace std;

namespace except = lsst::pex::exceptions;
namespace afwCoord = lsst::afw::coord;
namespace afwGeom = lsst::afw::geom;
namespace afwImg = lsst::afw::image;
namespace afwDet = lsst::afw::detection;
namespace afwMath = lsst::afw::math;

namespace {
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
    int nTerms = 0;
    for (int i = 1; i <= order; ++i) {
        nTerms += i;
    }

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
    if (A.rows() != b.rows()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "vector b of wrong size");        
    }
    Eigen::VectorXd par = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    return par;
}
}

///Constructor
CreateWcsWithSip::CreateWcsWithSip(
    afw::table::ReferenceMatchVector const & matches,
    CONST_PTR(lsst::afw::image::Wcs) linearWcs,
    int const order,
    lsst::afw::geom::Box2I const& bbox,
    int const ngrid
):
    _matches(matches),
    _bbox(bbox),
    _ngrid(ngrid),
    _linearWcs(linearWcs->clone()),
    _sipOrder(order+1),
    _reverseSipOrder(order+2), //Higher order for reverse transform
    _sipA(Eigen::MatrixXd::Zero(_sipOrder, _sipOrder)),
    _sipB(Eigen::MatrixXd::Zero(_sipOrder, _sipOrder)),
    _sipAp(Eigen::MatrixXd::Zero(_reverseSipOrder, _reverseSipOrder)),
    _sipBp(Eigen::MatrixXd::Zero(_reverseSipOrder, _reverseSipOrder)),
    _newWcs()
{
    if  (order < 2) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Sip matrices are at least 2nd order");        
    }

    if (_sipOrder > 9) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                          str(boost::format("SIP forward order %d exceeds the IAU limit of 9") % _sipOrder));
    }
    if (_reverseSipOrder > 9) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                          str(boost::format("SIP reverse order %d exceeds the IAU limit of 9") %
                              _reverseSipOrder));
    }
    
    if (_matches.size() < std::size_t(_sipOrder)) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Number of matches less than requested sip order");
    }

    if (_ngrid <= 0) {
        _ngrid = 3*_sipOrder;           // should be plenty
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
            afw::table::ReferenceMatchVector::const_iterator ptr = _matches.begin();
            ptr != _matches.end();
            ++ptr
        ) {
            lsst::afw::table::SourceRecord const & src = *ptr->second;
            _bbox.include(afwGeom::PointI(src.getX(), src.getY()));
        }
        float const borderFrac = 1/::sqrt(_matches.size()); // fractional border to add to exact BBox
        afwGeom::Extent2I border(borderFrac*_bbox.getWidth(), borderFrac*_bbox.getHeight());

        _bbox.grow(border);
    }
    // Calculate the forward part of the SIP distortion
    _calculateForwardMatrices();

    //Build a new wcs incorporating the forward sip matrix, it's all we know so far
    afwGeom::Point2D crval = _getCrvalAsGeomPoint();
    afwGeom::Point2D crpix = _linearWcs->getPixelOrigin();
    Eigen::MatrixXd CD = _linearWcs->getCDMatrix();
    
    _newWcs = PTR(afwImg::TanWcs)(new afwImg::TanWcs(crval, crpix, CD, _sipA, _sipB, _sipAp, _sipBp));
    // Use _newWcs to calculate the forward transformation on a grid, and derive the back transformation
    _calculateReverseMatrices();

    //Build a new wcs incorporating both of the sip matrices
    _newWcs = PTR(afwImg::TanWcs)(new afwImg::TanWcs(crval, crpix, CD, _sipA, _sipB, _sipAp, _sipBp));

}

void
CreateWcsWithSip::_calculateForwardMatrices()
{
    // Assumes FITS (1-indexed) coordinates.
    afwGeom::Point2D crpix = _linearWcs->getPixelOrigin();

    // Calculate u, v and intermediate world coordinates
    int const nPoints = _matches.size();
    Eigen::VectorXd u(nPoints), v(nPoints), iwc1(nPoints), iwc2(nPoints);

    int i = 0;
    for (
        afw::table::ReferenceMatchVector::const_iterator ptr = _matches.begin();
        ptr != _matches.end();
        ++ptr, ++i
    ) {
        afw::table::ReferenceMatch const & match = *ptr;

        // iwc's store the intermediate world coordinate positions of catalogue objects
        afwCoord::IcrsCoord c = match.first->getCoord();

        afwGeom::Point2D p = _linearWcs->skyToIntermediateWorldCoord(c);
        iwc1[i] = p[0];
        iwc2[i] = p[1];
        // u and v are intermediate pixel coordinates of observed (distorted) positions
        u[i] = match.second->getX() - crpix[0];
        v[i] = match.second->getY() - crpix[1];
    }
    
   
    // Forward transform
    int ord = _sipOrder;
    Eigen::MatrixXd forwardC = calculateCMatrix(u, v, ord);
    Eigen::VectorXd mu = leastSquaresSolve(iwc1, forwardC);
    Eigen::VectorXd nu = leastSquaresSolve(iwc2, forwardC);

    // Use mu and nu to refine CD

    // Given the implementation of indexToPQ(), the refined values
    // of the elements of the CD matrices are in elements 1 and "_sipOrder" of mu and nu.
    // If the implementation of indexToPQ() changes, these assertions
    // will catch that change.
    assert ((indexToPQ(0,   ord) == std::pair<int, int>(0, 0)));
    assert ((indexToPQ(1,   ord) == std::pair<int, int>(0, 1)));
    assert ((indexToPQ(ord, ord) == std::pair<int, int>(1, 0)));

    Eigen::Matrix2d CD;
    CD(1,0) = nu[ord];
    CD(1,1) = nu[1];
    CD(0,0) = mu[ord];
    CD(0,1) = mu[1];

    Eigen::Matrix2d CDinv = CD.inverse();   //Direct inverse OK for 2x2 matrix in Eigen

    // The zeroth elements correspond to a shift in crpix
    crpix[0] -= mu[0]*CDinv(0,0) + nu[0]*CDinv(0,1); 
    crpix[1] -= mu[0]*CDinv(1,0) + nu[0]*CDinv(1,1);

    afwGeom::Point2D crval = _getCrvalAsGeomPoint();
    _linearWcs = afwImg::Wcs::Ptr( new afwImg::Wcs(crval, crpix, CD));

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
            _sipA(p,q) = AB[0];
            _sipB(p,q) = AB[1];
        }
    }
}

void CreateWcsWithSip::_calculateReverseMatrices() {
    // Assumes FITS (1-indexed) coordinates.
    int const ngrid2 = _ngrid*_ngrid;

    Eigen::VectorXd u(ngrid2), v(ngrid2);
    Eigen::VectorXd U(ngrid2), V(ngrid2);
    Eigen::VectorXd delta1(ngrid2), delta2(ngrid2);
    
    int const x0 = _bbox.getMinX();
    float const dx = _bbox.getWidth()/(_ngrid - 1);
    int const y0 = _bbox.getMinY();
    float const dy = _bbox.getHeight()/(_ngrid - 1);

    afwGeom::Point2D crpix = _linearWcs->getPixelOrigin();
    int k = 0;
    for (int i = 0; i < _ngrid; ++i) {
        float const y = y0 + i*dy;
        for (int j = 0; j < _ngrid; ++j, ++k) {
            float const x = x0 + j*dx;
            // u and v are intermediate pixel coordinates on a grid of positions
            u[k] = x - crpix[0];
            v[k] = y - crpix[1];
            // U and V are the true, undistorted intermediate pixel positions as calculated
            // using the new Tan-Sip forward coefficients (to sky) and the linear Wcs (back to pixels)
            afwCoord::Coord::ConstPtr c = _newWcs->pixelToSky(x, y);
            afwGeom::Point2D p = _linearWcs->skyToPixel(*c);
            U[k] = p[0] - crpix[0];
            V[k] = p[1] - crpix[1];
            delta1[k] = u[k] - U[k];
            delta2[k] = v[k] - V[k];
        }
    }

    // Reverse transform
    int const ord = _reverseSipOrder;
    Eigen::MatrixXd reverseC = calculateCMatrix(U, V, ord); 
    Eigen::VectorXd tmpA = leastSquaresSolve(delta1, reverseC);
    Eigen::VectorXd tmpB = leastSquaresSolve(delta2, reverseC);
    
    assert(tmpA.rows() == tmpB.rows());
    for(int j=0; j< tmpA.rows(); ++j) {
        std::pair<int, int> pq = indexToPQ(j, ord);
        int p = pq.first, q = pq.second;

        _sipAp(p, q) = tmpA[j];
        _sipBp(p, q) = tmpB[j];   
    } 
}

///Get the scatter in position in pixel space 
double CreateWcsWithSip::getScatterInPixels() {
    vector<double> val;
    val.reserve(_matches.size());
    
    for (
        afw::table::ReferenceMatchVector::const_iterator ptr = _matches.begin();
        ptr != _matches.end();
        ++ptr
    ) {
        afw::table::ReferenceMatch const & match = *ptr;

        PTR(afw::table::SimpleRecord) catRec = match.first;
        PTR(afw::table::SourceRecord) srcRec = match.second;
        
        double imgX = srcRec->getX();
        double imgY = srcRec->getY();
        
        afwGeom::Point2D xy = _newWcs->skyToPixel(catRec->getCoord());
        double catX = xy[0];
        double catY = xy[1];
        
        val.push_back(::hypot(imgX - catX, imgY - catY));
   }
    
    return afwMath::makeStatistics(val, afwMath::MEDIAN).getValue();
}
    

///Get the scatter in (celestial) position
afwGeom::Angle CreateWcsWithSip::getScatterOnSky() {
    vector<double> val;
    val.reserve(_matches.size());

    for (
        afw::table::ReferenceMatchVector::const_iterator ptr = _matches.begin();
        ptr != _matches.end();
        ++ptr
    ) {
        afw::table::ReferenceMatch const & match = *ptr;

        PTR(afw::table::SimpleRecord) catRec = match.first;
        PTR(afw::table::SourceRecord) srcRec = match.second;
        afwCoord::IcrsCoord catRadec = catRec->getCoord();
        CONST_PTR(afwCoord::Coord) imgRadec = _newWcs->pixelToSky(srcRec->getCentroid());
        afwGeom::Angle sep = catRadec.angularSeparation(*imgRadec);
        val.push_back(sep.asDegrees());

    }
    assert(val.size() > 0);

    return afwMath::makeStatistics(val, afwMath::MEDIAN).getValue()*afwGeom::degrees;
}


afwGeom::Point2D CreateWcsWithSip::_getCrvalAsGeomPoint() {
    afwCoord::Fk5Coord coo = _linearWcs->getSkyOrigin()->toFk5();
    return coo.getPosition(afwGeom::degrees);
}

        
}}}}
