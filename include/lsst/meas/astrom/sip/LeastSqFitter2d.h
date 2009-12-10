// -*- LSST-C++ -*-

#ifndef LEAST_SQ_FITTER_2D
#define LEAST_SQ_FITTER_2D


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

///TODO: Implement method to return a afw::math::Function2d object
template <class FittingFunc>class LeastSqFitter2d {
public:
    LeastSqFitter2d(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, 
                    const std::vector<double> &s, int order);

    Eigen::MatrixXd getParams();
    Eigen::MatrixXd getErrors();
    double valueAt(double x, double y);
    std::vector<double> residuals();
    
    void printParams();
    double getChiSq();
    double getReducedChiSq();


private:
    void initFunctions();
             
    void calculateA();
    void calculateBeta();
    double func2d(double x, double y, int term);
    double func1d(double value, int exponent);
    
    std::vector<double> _x, _y, _z, _s;
    unsigned int _order;   //Degree of polynomial to fit, e.g 4=> cubic
    unsigned int _nPar;    //Number of parameters in fitting eqn, e.g x^2, xy, y^2, x^3, 
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

///Fit a 2d polynomial to a set of data points z(x, y)
///
///\tparam FittingFunc  The type of function to fit in each dimension. This function should extend the base
///       class of lsst::afw::math::Function1. Note that this is a 1d function
///\param x vector of x positions of data
///\param y vector of y positions of data
///\param z Value of data for a given x,y. z = z_i = z_i(x_i, y_i)
///\param s Vector of measured uncertainties in the values of z
///\param order Order of 2d function to fit
template<class FittingFunc> LeastSqFitter2d<FittingFunc>::LeastSqFitter2d(const std::vector<double> &x, 
    const std::vector<double> &y, const std::vector<double> &z, const std::vector<double> &s, int order) :
    _x(x), _y(y), _z(z), _s(s), _order(order), _nPar(0), _A(1,1), _beta(1), _par(1) {
    
    //_nPar, the number of terms to fix (x^2, xy, y^2 etc.) is \Sigma^(order+1) 1
    _nPar = 0;
    for (int i = 0; i<order; ++i) {
        _nPar += i + 1;
    }
    
    //Check input vectors are the same size
    _nData = _x.size();
    if (_nData != _y.size()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "x and y vectors of different lengths");        
    }
    if (_nData != _s.size()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "x and s vectors of different lengths");        
    }
    if (_nData != _z.size()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "x and z vectors of different lengths");        
    }

    for (unsigned int i = 0; i<_nData; ++i) {
        if ( _s[i] == 0 ) {
            std::string msg = "Illegal zero value for fit weight encountered.";
            throw LSST_EXCEPT(except::RuntimeErrorException, msg);        
        }
    }

    if (_nData < _order) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Fewer data points than parameters");        
    }
    
    
    initFunctions();
    calculateBeta();
    calculateA();

    
    //Try three different methods of solving the linear equation
    _par = Eigen::VectorXd(_nPar);
    if (! _A.ldlt().solve(_beta, &_par)) {
         pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter2d",
                           "Unable fit data with Cholesky LDL^t");

        if (! _A.llt().solve(_beta, &_par)) {
             pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter2d",
                           "Unable fit data with Cholesky LL^t either");
                        
            if (! _A.lu().solve(_beta, &_par)) {
                 pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter2d",
                               "Unable fit data with LU decomposition either");

                 throw LSST_EXCEPT(pexExcept::Exception,
                     "Unable to solve least squared equation in LeastSqFitter2d()");
            }
        }
    }
    
}

        
///Build up a triangular matrix of the parameters. The shape of the matrix is
///such that the values correspond to the coefficients of the following polynomials
/// 1   y    y^2  y^3
/// x   xy   xy^2 0
/// x^2 x^2y 0    0
/// x^3 0    0    0   (order==4)
///
/// where row*column < _order
template<class FittingFunc> Eigen::MatrixXd LeastSqFitter2d<FittingFunc>::getParams() {

    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(_order, _order);  //Should be a shared ptr?

    int count = 0;
    for (unsigned int i = 0; i< _order; ++i) {
        for (unsigned int j = 0; j< _order - i; ++j) {
            out(i, j) = _par(count++);
        }
    }
    return out;
}


///Print the parameters of the fit matrix. This function is helpful in debugging 
///because we haven't wrapped the Eigen::Matrix class into Python, so we can't 
///extract values easily.
template<class FittingFunc> void LeastSqFitter2d<FittingFunc>::printParams() {

    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(_order, _order);  //Should be a shared ptr?

    int count = 0;
    for (unsigned int i = 0; i< _order; ++i) {
        for (unsigned int j = 0; j< _order - i; ++j) {
            printf("%7e  ", _par(count++));
        }
        printf("\n");
    }
}


/// \brief Return a measure of the goodness of fit. 
/// \f[ \chi^2 = \sum \left( \frac{z_i - f(x_i, y_i)}{s_i} \right)^2  \f]
template<class FittingFunc> double LeastSqFitter2d<FittingFunc>::getChiSq() {
    
    double chisq = 0;
    for (unsigned int i = 0; i< _nData; ++i) {
        double val = _z[i] - valueAt(_x[i], _y[i]);
        val /= _s[i];
        chisq += pow(val, 2);
    }
    
    return chisq;
}


/// \brief Return a measure of the goodness of fit. 
/// \f[ \chi_r^2 = \sum \left( \frac{z_i - f(x_i, y_i)}{s} \right)^2 \div (N-p) \f]
/// Where \f$ N \f$ is the number of data points, and \f$ p \f$ is the number 
/// of parameters in the fit
/// 
template<class FittingFunc> double LeastSqFitter2d<FittingFunc>::getReducedChiSq() {
    return getChiSq()/(double) (_nData - _nPar);
}



///Return the value of the best fit function at a given position (x,y)
template<class FittingFunc>  double LeastSqFitter2d<FittingFunc>::valueAt(double x, double y){
    double val = 0;
    
    //Sum the values of the different orders to get the value of the fitted function
    for (unsigned int i = 0; i< _nPar; ++i) {
        val += _par[i] * func2d(x, y, i);
    }
    return val;
}

///Return a vector of residuals of the fit (i.e the difference between the input z values, and 
///the value of the fitting function at that point.
template<class FittingFunc>  std::vector<double> LeastSqFitter2d<FittingFunc>::residuals(){

    std::vector<double> out;
    out.reserve(_nData);
    
    for (unsigned int i; i< _nData; ++i) {
        out.push_back(_z[i] - valueAt(_x[i], _y[i]));
    }
    
    return out;
}



template<class FittingFunc> Eigen::MatrixXd LeastSqFitter2d<FittingFunc>::getErrors() {
    Eigen::MatrixXd V = _A.svd().matrixV();
    Eigen::VectorXd w = _A.svd().singularValues();

    Eigen::VectorXd err(_nPar);
    for (unsigned int i = 0; i<_nPar; ++i) {
        err(i) = 0;
        for (unsigned int j = 0; j<_nPar; ++j) {
            double val =  V(i, j)/w(j);
            err(i) += val*val;
        }
    }
    return err;
}

    
template<class FittingFunc> void LeastSqFitter2d<FittingFunc>::initFunctions() {
    //Initialise the array of functions. _funcArray[i] is a object of type math::Function1 of order i
    _funcArray.reserve(_order);
    
    std::vector<double> coeff;
    coeff.reserve( _order);
    
    coeff.push_back(1);
    for (unsigned int i = 0; i< _order; ++i) {
        boost::shared_ptr<FittingFunc> p(new FittingFunc(coeff));
        _funcArray.push_back(p);
        
        coeff[i] = 0;
        coeff.push_back(1);  //coeff now looks like [0,0,...,0,1]
    }
}
    

template<class FittingFunc> void LeastSqFitter2d<FittingFunc>::calculateA() {

    assert(_nPar != 0);
    _A = Eigen::MatrixXd(static_cast<int> (_nPar), static_cast<int> (_nPar));

    for (unsigned int i = 0; i< _nPar; ++i) {
        for (unsigned int j = 0; j< _nPar; ++j) {
            double val = 0;
            for (unsigned int k = 0; k< _nData; ++k) {
                val += func2d(_x[k], _y[k], i) * func2d(_x[k], _y[k], j)/( _s[k]*_s[k]);
            }
            _A(i, j) = val;
        }
        
    }        
}

    
template<class FittingFunc> void LeastSqFitter2d<FittingFunc>::calculateBeta() {

    assert(_nPar != 0);
    _beta = Eigen::VectorXd(_nPar);

    double val;
    for (unsigned int i = 0; i< _nPar; ++i) {
        _beta(i) = 0;
        for (unsigned int j = 0; j< _nData; ++j) {
            val = _z[j]*func2d(_x[j], _y[j], i)/ (_s[i]*_s[i]);
            _beta(i) += val;
        }
    }
}
 
 

//The ith term in the fitting polynomial is of the form x^a * y^b. This function figures
//out the value of a and b, then calculates the value of the ith term at the given x and y
template<class FittingFunc> double LeastSqFitter2d<FittingFunc>::func2d(double x, double y, int term) {

    int yexp = 0;  //y exponent
    int xexp = 0;  //x exponent

    for (int i = 0; i<term; ++i) {
        yexp = (yexp + 1) % (_order - xexp);
        if ( yexp == 0){
            xexp++;
        }
    }

    double xcomp = func1d(x, xexp);  //x component of polynomial
    double ycomp = func1d(y, yexp);  //y component

    #if 0   //A useful debugging printf statement
        printf("The %i(th) function: x^%i * y^%i = %.1f * %.1f\n", term, xexp, yexp, xcomp, ycomp);
    #endif
    
    return xcomp*ycomp;
}

template<class FittingFunc> double LeastSqFitter2d<FittingFunc>::func1d(double value, int exponent) {
    return (*_funcArray[exponent])(value);
}


    
}}}}

#endif
