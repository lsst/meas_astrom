// -*- LSST-C++ -*-

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

template <class FittingFunc>class LeastSqFitter1d {
public:
    LeastSqFitter1d(const std::vector<double> &x, const std::vector<double> &y, 
                    const std::vector<double> &s,  
                    unsigned int order);

    Eigen::VectorXd getParams();
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
    

    std::vector<double> _x, _y, _z, _s;
    unsigned int _order;   //Degree of polynomial to fit, e.g 4=> cubic
    unsigned int _nData;   //Number of data points, == _x.size()
    
    Eigen::MatrixXd _A;
    Eigen::VectorXd _beta;
    Eigen::MatrixXd _Ainv;
    Eigen::VectorXd _par;

    std::vector<boost::shared_ptr<FittingFunc> > _funcArray;
};
    



//The .cc part
namespace sip = lsst::meas::astrom::sip;
namespace math = lsst::afw::math;


///Fit a 1d polynomial to a set of data points z(x, y)
///
///\tparam FittingFunc  The type of function to fit. This function should extend the base
///       class of lsst::afw::math::Function1. 
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
    if (! _A.ldlt().solve(_beta, &_par)) {
         pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter1d",
                           "Unable fit data with Cholesky LDL^t");

        if (! _A.llt().solve(_beta, &_par)) {
             pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter1d",
                           "Unable fit data with Cholesky LL^t either");
                        
            if (! _A.lu().solve(_beta, &_par)) {
                 pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter1d",
                               "Unable fit data with LU decomposition either");

                 throw LSST_EXCEPT(pexExcept::Exception,
                     "Unable to solve least squares equation in LeastSqFitter1d()");
            }
        }
    }
}
        


///Return the best fit paramets as an Eigen::Matrix
template<class FittingFunc> Eigen::VectorXd LeastSqFitter1d<FittingFunc>::getParams() {

    Eigen::VectorXd vec = Eigen::VectorXd::Zero(_order);
    for (unsigned int i = 0; i< _order; ++i) {
        vec(i) = _par(i);
    }
    return vec;
}


///Return the best fit polynomial
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
    
    for (unsigned int i; i< _nData; ++i) {
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
    for (unsigned int i = 0; i< _nData; ++i) {
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
            for (unsigned int k = 0; k< _nData; ++k) {
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
    for (unsigned int i = 0; i< _order; ++i) {
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
