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
 

#ifndef LEAST_SQ_FITTER_1D
#define LEAST_SQ_FITTER_1D

#include <cstdio>
#include <vector>

#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"
#include "Eigen/SVD"
#include "Eigen/Cholesky"
#include "Eigen/LU"

#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/math/FunctionLibrary.h"

namespace except = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;


namespace lsst { 
namespace meas { 
namespace astrom { 
namespace sip {

/**
\brief Fit an lsst::afw::math::Function1 object to a set of data points in one dimension

The class is templated over the kind of object to fit. 

Input is a list of x ordinates for a set of points, the y coordinate, and the uncertainties, s.
order is order of the polynomial to fit (e.g if the templated function is 
lsst::afw::math::PolynomialFunction1, then order=3 => fit a function of the form \f$ax^2+bx+c\f$

\tparam FittingFunc The 1d function to fit in both dimensions. Must inherit from 
lsst::afw::math::Function1 

\param x Ordinate of points to fit
\param y Co-ordinate of pionts to fit
\param s 1\f$\sigma\f$ uncertainties in z
\param order Polynomial order to fit

\sa LeastSqFitter1d
*/

template <class FittingFunc>class LeastSqFitter1d {
public:
    LeastSqFitter1d(const std::vector<double> &x, const std::vector<double> &y, 
                    const std::vector<double> &s,  
                    unsigned int order);

    Eigen::VectorXd getParams();
    Eigen::VectorXd getErrors();
    FittingFunc getBestFitFunction();
    double valueAt(double x);
    std::vector<double> residuals();
    
    double getChiSq();
    double getReducedChiSq();

    
private:
    void initFunctions();
             
    void calculateA();
    void calculateBeta();
    double func1d(double value, int exponent);
    

    std::vector<double> _x, _y, _s;
    unsigned int _order;   //Degree of polynomial to fit, e.g 4=> cubic
    unsigned int _nData;   //Number of data points, == _x.size()
    
    Eigen::MatrixXd _A;
    Eigen::VectorXd _beta;
    Eigen::VectorXd _par;

    std::vector<boost::shared_ptr<FittingFunc> > _funcArray;
};
    



//The .cc part
namespace sip = lsst::meas::astrom::sip;
namespace math = lsst::afw::math;


///Fit a 1d polynomial to a set of data points z(x, y)
///
///\tparam FittingFunc  The type of function to fit. This function extends the base
///       class of lsst::afw::math::Function1
///\param x vector of x positions of data
///\param y vector of y positions of data
///\param s Vector of measured uncertainties in the values of z
///\param order Order of 2d function to fit
template<class FittingFunc> LeastSqFitter1d<FittingFunc>::LeastSqFitter1d(const std::vector<double> &x, 
    const std::vector<double> &y, const std::vector<double> &s, unsigned int order) :
    _x(x), _y(y), _s(s), _order(order), _A(1,1), _beta(1) {
    
    if (order == 0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Fit order must be >= 1");        
    }
    
    _nData = _x.size();
    if  (_nData != _y.size()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "x and y vectors of different lengths");        
    }
    if (_nData != _s.size()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "x and s vectors of different lengths");        
    }

    if (_nData < _order) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Fewer data points than parameters");        
    }


    initFunctions();
    calculateBeta();
    calculateA();


    //Try three different methods of solving the linear equation
    _par = Eigen::VectorXd(_order);
    Eigen::FullPivLU<Eigen::MatrixXd> lu(_A);
    _par = _A.ldlt().solve(_beta);
}
        


///Return the best fit parameters as an Eigen::Matrix
template<class FittingFunc> Eigen::VectorXd LeastSqFitter1d<FittingFunc>::getParams() {

    Eigen::VectorXd vec = Eigen::VectorXd::Zero(_order);
    for (unsigned int i = 0; i< _order; ++i) {
        vec(i) = _par(i);
    }
    return vec;
}


///Return the 1 sigma uncertainties in the best fit parameters as an Eigen::Matrix
template<class FittingFunc> Eigen::VectorXd LeastSqFitter1d<FittingFunc>::getErrors() {

    Eigen::MatrixXd Ainv = _A.partialPivLU().inverse();
    Eigen::VectorXd vec = Eigen::VectorXd::Zero(_order);
    for (unsigned int i = 0; i< _order; ++i) {
        vec(i) = std::sqrt(Ainv(i,i));
    }
    return vec;
}


///Return the best fit polynomial as a lsst::afw::math::Function1 object
template<class FittingFunc> FittingFunc LeastSqFitter1d<FittingFunc>::getBestFitFunction() {

    //FittingFunc and LeastSqFitter disagree on the definition of order of a function.
    //LSF says that a linear function is order 2 (two coefficients), FF says only 1
    FittingFunc func(_order - 1);

    for (unsigned int i = 0; i< _order; ++i) {
        func.setParameter(i, _par(i));
    }
    return func;
}


///Calculate the value of the function at a given point
template<class FittingFunc> double LeastSqFitter1d<FittingFunc>::valueAt(double x) {
        FittingFunc f = getBestFitFunction();
        
        return f(x);
}



///Return a vector of residuals of the fit (i.e the difference between the input y values, and 
///the value of the fitting function at that point.
template<class FittingFunc>  std::vector<double> LeastSqFitter1d<FittingFunc>::residuals(){
    std::vector<double> out;
    out.reserve(_nData);
    
    FittingFunc f = getBestFitFunction();
    
    for (unsigned int i = 0; i< _nData; ++i) {
        out.push_back(_y[i] - f(_x[i]));
    }
    
    return out;
}


/// \brief Return a measure of the goodness of fit. 
/// \f[ \chi_r^2 = \sum \left( \frac{y_i - f(x_i)}{s} \right)^2 \f]
/// 
template<class FittingFunc> double LeastSqFitter1d<FittingFunc>::getChiSq() {
    FittingFunc f = getBestFitFunction();

    double chisq = 0;
    for (unsigned int i = 0; i < _nData; ++i) {
        double val = _y[i] - f(_x[i]);
        val /= _s[i];
        chisq += pow(val, 2);
    }

    return chisq;
}


/// \brief Return a measure of the goodness of fit. 
/// \f[ \chi_r^2 = \sum \left( \frac{y_i - f(x_i)}{s} \right)^2 \div (N-p) \f]
/// Where \f$ N \f$ is the number of data points, and \f$ p \f$ is the number 
/// of parameters in the fit
/// 
template<class FittingFunc> double LeastSqFitter1d<FittingFunc>::getReducedChiSq() {
    return getChiSq()/(double) (_nData - _order);
}



/// Initialise the array of functions. _funcArray[i] is a object of type math::Function1 of order 
///_norder
template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::initFunctions() {

    _funcArray.reserve(_order);
    
    std::vector<double> coeff;
    coeff.reserve( _order);
    
    coeff.push_back(1.0);
    for ( unsigned int i = 0; i< _order; ++i) {
        boost::shared_ptr<FittingFunc> p(new FittingFunc(coeff));
        _funcArray.push_back(p);
        coeff[i] = 0.0;
        coeff.push_back(1.0);  //coeff now looks like [0,0,...,0,1]
    }
}
    

template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::calculateA() {

    _A = Eigen::MatrixXd( static_cast<int> (_order), static_cast<int> (_order));

    unsigned int i, j;
    for (i = 0; i< _order; ++i) {
        for (j = 0; j< _order; ++j) {
            double val = 0;
            for (unsigned int k = 0; k < _nData; ++k) {
                val += func1d(_x[k], i) * func1d(_x[k], j)/( _s[k]*_s[k]);
            }
            _A(i, j) = val;
        }
        
    }        

}

    
template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::calculateBeta() {

    _beta = Eigen::VectorXd(_order);

    double val;
    unsigned int j;
    for (unsigned int i = 0; i < _order; ++i) {
        _beta(i) = 0;
        for (j = 0; j< _nData; ++j) {
            val = _y[j]*func1d(_x[j], i)/ (_s[i]*_s[i]);
            _beta(i) += val;
        }
    }
}
 
 

template<class FittingFunc> double LeastSqFitter1d<FittingFunc>::func1d(double value, int exponent) {
    return (*_funcArray[exponent])(value);
}


}}}}

#endif
