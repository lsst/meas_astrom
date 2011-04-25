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
 

#include "lsst/meas/astrom/sip/CreateWcsWithSip.h"

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace sip {

using namespace std;

namespace except = lsst::pex::exceptions;
namespace pexLog = lsst::pex::logging;
namespace afwCoord = lsst::afw::coord;
namespace afwGeom = lsst::afw::geom;
namespace afwImg = lsst::afw::image;
namespace det = lsst::afw::detection;
namespace math = lsst::afw::math;




///Constructor
CreateWcsWithSip::CreateWcsWithSip(const std::vector<lsst::afw::detection::SourceMatch> match,
                     const lsst::afw::image::Wcs::Ptr linearWcs,
                     int order):
                     _matchList(match), 
                     _linearWcs(linearWcs->clone()),
                     _sipOrder(order+1),
                     _reverseSipOrder(order+2), //Higher order for reverse transform
                     _size(match.size()),
                     _sipA(Eigen::MatrixXd::Zero(_sipOrder, _sipOrder)),
                     _sipB(Eigen::MatrixXd::Zero(_sipOrder, _sipOrder)),
                     _sipAp(Eigen::MatrixXd::Zero(_reverseSipOrder, _reverseSipOrder)),
                     _sipBp(Eigen::MatrixXd::Zero(_reverseSipOrder, _reverseSipOrder)),
                     _newWcs() {

    if  (order < 2) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Sip matrices are at least 2nd order");        
    }
    
    if (_size < _sipOrder) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Number of matches less than requested sip order");
    }


    _calculateForwardMatrices();
    _calculateReverseMatrices();

    //Build a new wcs incorporating the sip matrices
    afwGeom::Point2D crval = _getCrvalAsGeomPoint();
    afwGeom::Point2D crpix = _linearWcs->getPixelOrigin();
    Eigen::MatrixXd CD = _linearWcs->getCDMatrix();
    
    _newWcs = afwImg::TanWcs::Ptr(new afwImg::TanWcs(crval, crpix, CD, _sipA, _sipB, _sipAp, _sipBp));

}


afwImg::TanWcs::Ptr CreateWcsWithSip::getNewWcs() {
    return _newWcs;
}


int CreateWcsWithSip::getUIndex(int j, int order) {
    int decrement = order;
    int k=0;
    
    while ((j>= decrement) && (decrement > 0)) {
        j -= decrement;
        decrement--;
        k++;
    }
    
    return k;
}

int CreateWcsWithSip::getVIndex(int j, int order) {
    int decrement = order;
    
    while ((j>= decrement) && (decrement > 0)) {
        j -= decrement;
        decrement--;
    }
    return j;
    
}


void CreateWcsWithSip::_calculateForwardMatrices() {
    // Assumes FITS (1-indexed) coordinates.
    afwGeom::Point2D crpix = _linearWcs->getPixelOrigin();

    // Calculate u, v and intermediate world coordinates
    Eigen::VectorXd u(_size), v(_size), iwc1(_size), iwc2(_size);
    
    for(int i=0; i < _size; ++i) {
        // iwc's store the intermediate world coordinate positions of catalogue objects
        afwCoord::Coord::Ptr c = _matchList[i].first->getRaDec();
        afwGeom::Point2D p = _linearWcs->skyToIntermediateWorldCoord(c);
        iwc1[i] = p[0];
        iwc2[i] = p[1];
        // u and v are intermediate pixel coordinates of observed (distorted) positions
        u[i] = _matchList[i].second->getXAstrom() - crpix[0];
        v[i] = _matchList[i].second->getYAstrom() - crpix[1];
    }
    
   
    // Forward transform
    int ord = _sipOrder;
    Eigen::MatrixXd forwardC = _calculateCMatrix(u, v, ord);
    Eigen::VectorXd mu = _leastSquaresSolve(iwc1, forwardC);
    Eigen::VectorXd nu = _leastSquaresSolve(iwc2, forwardC);

    // Use mu and nu to refine CD

    // Given the implementation of getUIndex() and getVIndex(), the refined values
    // of the elements of the CD matrices are in elements 1 and "_sipOrder" of mu and nu.
    // If the implementation of getUIndex() and getVIndex() change, these assertions
    // will catch that change.
    assert(getUIndex(0, ord) == 0 && getVIndex(0, ord) == 0);
    assert(getUIndex(1, ord) == 0 && getVIndex(1, ord) == 1);
    assert(getUIndex(ord, ord) == 1 && getVIndex(ord, ord) == 0);

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
        int p = getUIndex(i, ord);
        int q = getVIndex(i, ord);
        if( (p+q > 1) && (p+q < ord)) {    
            Eigen::Vector2d munu(2,1);
            munu(0) = mu(i);
            munu(1) = nu(i);
            Eigen::Vector2d AB = CDinv * munu;
            _sipA(p,q) = AB[0];
            _sipB(p,q) = AB[1];
        }
    }
}




    
void CreateWcsWithSip::_calculateReverseMatrices() {
    // Assumes FITS (1-indexed) coordinates.
    afwGeom::Point2D crpix = _linearWcs->getPixelOrigin();

    Eigen::VectorXd u(_size), v(_size);
    Eigen::VectorXd U(_size), V(_size);
    Eigen::VectorXd delta1(_size), delta2(_size);
    
    for(int i=0; i < _size; ++i) {
        // u and v are intermediate pixel coordinates of observed (distorted) positions
        u[i] = _matchList[i].second->getXAstrom() - crpix[0];
        v[i] = _matchList[i].second->getYAstrom() - crpix[1];
        // U and V are the true, undistorted intermediate pixel positions as calculated
        // from the catalogue ra and decs and the (updated) linear wcs
        afwCoord::Coord::Ptr c = _matchList[i].first->getRaDec();
        afwGeom::Point2D p = _linearWcs->skyToPixel(c);
        U[i] = p[0] - crpix[0];
        V[i] = p[1] - crpix[1];
        delta1[i] = u[i] - U[i];
        delta2[i] = v[i] - V[i];
    }
    
    
    // Reverse transform
    int ord = _reverseSipOrder;
    Eigen::MatrixXd reverseC = _calculateCMatrix(U, V, ord); 
    Eigen::VectorXd tmpA = _leastSquaresSolve(delta1, reverseC);
    Eigen::VectorXd tmpB = _leastSquaresSolve(delta2, reverseC);
    
    assert(tmpA.rows() == tmpB.rows());
    for(int j=0; j< tmpA.rows(); ++j) {
        int p = getUIndex(j, ord);
        int q = getVIndex(j, ord);
        _sipAp(p,q) = tmpA[j];
        _sipBp(p,q) = tmpB[j];   
    } 
}



///Get the scatter in position in pixel space 
double CreateWcsWithSip::getScatterInPixels() {
    unsigned int size = _matchList.size();
    
    vector<double> val;
    val.reserve(size);
    
    for (unsigned int i = 0; i< size; ++i) {
        det::Source::Ptr catSrc = _matchList[i].first;
        det::Source::Ptr imgSrc = _matchList[i].second;
        
        double imgX = imgSrc->getXAstrom();
        double imgY = imgSrc->getYAstrom();
        
        afwGeom::Point2D xy = _newWcs->skyToPixel(catSrc->getRaDec());
        double catX = xy[0];
        double catY = xy[1];
        
        val.push_back(hypot(imgX - catX, imgY - catY));
   }
    
    math::Statistics stat = math::makeStatistics(val, math::MEDIAN);
    double scatter = stat.getValue(math::MEDIAN);
    
    return scatter;
}
    

///Get the scatter in position in arcseconds
double CreateWcsWithSip::getScatterInArcsec() {
    unsigned int size = _matchList.size();
    vector<double> val;
    val.reserve(size);
    for (unsigned int i = 0; i < size; ++i) {
        det::Source::Ptr catSrc = _matchList[i].first;
        det::Source::Ptr imgSrc = _matchList[i].second;
        afwCoord::Coord::ConstPtr catRadec = catSrc->getRaDec();
        afwCoord::Coord::ConstPtr imgRadec = _newWcs->pixelToSky(imgSrc->getXAstrom(), imgSrc->getYAstrom());
        val.push_back(catRadec->angularSeparation(*imgRadec, afwCoord::DEGREES));
        /*
         printf("ras: (%.3f, %.3f), decs: (%.3f, %.3f), dist: %.3f deg / %.3f deg.\n",
         catRadec->toFk5().getRa(afwCoord::DEGREES),
         imgRadec->toFk5().getRa(afwCoord::DEGREES),
         catRadec->toFk5().getDec(afwCoord::DEGREES),
         imgRadec->toFk5().getDec(afwCoord::DEGREES),
         catRadec->angularSeparation(*imgRadec, afwCoord::DEGREES),
         imgRadec->angularSeparation(*catRadec, afwCoord::DEGREES));
         */
    }
    assert(val.size() > 0);
    math::Statistics stat = math::makeStatistics(val, math::MEDIAN);
    double scatter = stat.getValue(math::MEDIAN);
    return scatter*3600;    //Convert to arcsec
}


Eigen::MatrixXd CreateWcsWithSip::_calculateCMatrix(Eigen::VectorXd axis1, Eigen::VectorXd axis2, int order) {

    int nTerms = 0;
    for(int i =1; i<= order; ++i) {
        nTerms += i;
    }
    
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(_size, nTerms);
    for(int i=0; i< _size; ++i) {
        for(int j=0; j < nTerms; ++j) {
            int p = getUIndex(j, order);
            int q = getVIndex(j, order);
            assert(p+q < order);
            
            C(i,j) = pow(axis1[i], p) * pow(axis2[i], q);
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
Eigen::VectorXd CreateWcsWithSip::_leastSquaresSolve(Eigen::VectorXd b, Eigen::MatrixXd A) {

    int rowsA = A.rows();
    int colsA = A.cols();
    int rowsB = b.rows();
    
    if  (rowsA != rowsB) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "vector b of wrong size");        
    }

    Eigen::MatrixXd Atr = A.transpose();
    Eigen::MatrixXd AtrA = Atr * A;  //A transpose x A
    Eigen::MatrixXd Atrb = Atr * b; //A transpose x b

    //Try three different methods of solving the linear equation
    Eigen::VectorXd par(colsA);
    if (! AtrA.ldlt().solve(Atrb, &par)) {
         pexLog::TTrace<5>("lsst.meas.astrom.sip.LeastSquaresSolve",
                           "Unable fit data with Cholesky LDL^t");

        if (! AtrA.llt().solve(Atrb, &par)) {
             pexLog::TTrace<5>("lsst.meas.astrom.sip.LeastSquaresSolve",
                           "Unable fit data with Cholesky LL^t either");
                        
            if (! AtrA.lu().solve(Atrb, &par)) {
                 pexLog::TTrace<5>("lsst.meas.astrom.sip.LeastSquaresSolve",
                               "Unable fit data with LU decomposition either");

                 throw LSST_EXCEPT(pexExcept::Exception,
                     "Unable to solve least squares equation in LeastSquaresSolve()");
            }
        }
    }
    
    return par;
}



afwGeom::Point2D CreateWcsWithSip::_getCrvalAsGeomPoint() {
    afwCoord::Fk5Coord coo = _linearWcs->getSkyOrigin()->toFk5();
    // NOTE, these are in degrees to agree with Wcs / TanWcs
    double ra =coo.getRa(afwCoord::DEGREES);
    double dec=coo.getDec(afwCoord::DEGREES);
    return afwGeom::Point2D(ra,dec);
}

        
}}}}
