
using namespace std;

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


namespace lsst { namespace meas { namespace astrom { namespace sip {

template <class FittingFunc>class LeastSqFitter1d
{
public:
    LeastSqFitter1d(const vector<double> &x, const vector<double> &y, const vector<double> &s,  unsigned int order);

    Eigen::VectorXd getParams();
    FittingFunc getBestFitFunction();
    double valueAt(double x);
    
private:
    void initFunctions();
             
    void calculateA();
    void calculateBeta();
    double func1d(double value, int exponent);
    

    vector<double> _x, _y, _z, _s;
    unsigned int _order;   //Degree of polynomial to fit, e.g 4=> cubic
    unsigned int _nData;   //Number of data points, == _x.size()
    
    Eigen::MatrixXd _A;
    Eigen::VectorXd _beta;
    Eigen::MatrixXd _Ainv;
    Eigen::VectorXd _par;

    vector<boost::shared_ptr<FittingFunc> > _funcArray;
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
///\param z Value of data for a given x,y. z = z_i = z_i(x_i, y_i)
///\param s Vector of measured uncertainties in the values of z
///\param order Order of 2d function to fit
template<class FittingFunc> LeastSqFitter1d<FittingFunc>::LeastSqFitter1d(const vector<double> &x, const vector<double> &y, const vector<double> &s, unsigned int order) :
    _x(x), _y(y), _s(s), _order(order), _A(1,1), _beta(1) {
    
    if(order == 0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Fit order must be >= 1");        
    }
    
    _nData = _x.size();
    if(_nData != _y.size()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "x and y vectors of different lengths");        
    }
    if(_nData != _s.size()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "x and s vectors of different lengths");        
    }

    if(_nData < _order) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Fewer data points than parameters");        
    }

printf("\nInput data\n");
for(int i=0; i<_nData; ++i) {
    printf("%.16f %.16f %.16f\n", _x[i], _y[i], _s[i]);
}

    initFunctions();
    calculateBeta();
    calculateA();

printf("A\n");
std::cout << _A << std::endl;    
printf("beta\n");
std::cout << _beta << std::endl;    

    //Try three different methods of solving the linear equation
    _par = Eigen::VectorXd(_order);
    if(! _A.ldlt().solve(_beta, &_par)) {
         pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter1d",
                           "Unable fit data with Cholesky LDL^t");

        if(! _A.llt().solve(_beta, &_par)) {
             pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter1d",
                           "Unable fit data with Cholesky LL^t either");
                        
            if(! _A.lu().solve(_beta, &_par)) {
                 pexLogging::TTrace<5>("lsst.meas.astrom.sip.LeastSqFitter1d",
                               "Unable fit data with LU decomposition either");

                 throw LSST_EXCEPT(pexExcept::Exception,
                     "Unable to determine kernel solution in LeastSqFitter1d()");
            }
        }
    }
        
    
printf("Par\n");
std::cout << _par << std::endl;
    
}
        


///Return the best fit paramets as an Eigen::Matrix
template<class FittingFunc> Eigen::VectorXd LeastSqFitter1d<FittingFunc>::getParams() {

    Eigen::VectorXd vec = Eigen::VectorXd::Zero(_order);
    for(unsigned int i=0; i< _order; ++i) {
        vec(i) = _par(i);
        printf("%g ", vec(i));
    }
    printf("\n");
    return vec;
}


///Return the best fit polynomial
template<class FittingFunc> FittingFunc LeastSqFitter1d<FittingFunc>::getBestFitFunction() {

    //FittingFunc and LeastSqFitter disagree on the definition of order of a function.
    //LSF says that a linear function is order 2 (two coefficients), FF says only 1
    FittingFunc func(_order-1);

    for(unsigned int i=0; i< _order; ++i) {
        func.setParameter(i, _par(i));
        //printf("bff %i %g\n", i, func.getParameter(i));
    }
    return func;
}


///Calculate the value of the function at a given point
template<class FittingFunc> double LeastSqFitter1d<FittingFunc>::valueAt(double x) {
        FittingFunc f = getBestFitFunction();
        
        return f(x);
}


/// Initialise the array of functions. _funcArray[i] is a object of type math::Function1 of order 
///_norder
template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::initFunctions() {

    _funcArray.reserve(_order);
    
    vector<double> coeff;
    coeff.reserve( _order);
    
    coeff.push_back(1.);
    for( unsigned int i=0; i< _order; ++i) {
        boost::shared_ptr<FittingFunc> p(new FittingFunc(coeff));
        _funcArray.push_back(p);

printf("Init i=%i: Coeff= ", i);
for(int j=0; j<=i; ++j) {
    printf("%.0f ", coeff[j]);
}
printf("\n");

        coeff[i] = 0.;
        coeff.push_back(1.);  //coeff now looks like [0,0,...,0,1]
    }
}
    

template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::calculateA() {

    _A = Eigen::MatrixXd( (int) _order, (int) _order);

    unsigned int i, j;
    for(i=0; i< _order; ++i) {
        for(j=0; j< _order; ++j) {
            double val=0;
            for(unsigned int k=0; k< _nData; ++k) {
                val += func1d(_x[k], i) * func1d(_x[k], j)/( _s[k]*_s[k]);
printf("(i,j,k)=(%i,%i,%i) x=%.1f y=%.1f fi=%.1f fj=%.1f, val=%.1f\n", i, j, k, _x[j], _y[j], func1d(_x[j], i), func1d(_x[j], j), val);
            }
            _A(i,j) = val;
        }
        
    }        
std::cout << _A << std::endl;
}

    
template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::calculateBeta() {

    _beta = Eigen::VectorXd(_order);

printf("CalculateBeta()\n");
    double val;
    unsigned int j;
    for(unsigned int i=0; i< _order; ++i) {
        _beta(i) = 0;
        for(j=0; j< _nData; ++j) {
            val = _y[j]*func1d(_x[j], i)/ (_s[i]*_s[i]);
            _beta(i) += val;
printf("(i,j)=(%i,%i) x=%.1f y=%.1f f=%.1f, val=%.1f\n", i, j, _x[j], _y[j], func1d(_x[j], i), val);
        }
    }
std::cout << _beta << std::endl;
}
 
 

template<class FittingFunc> double LeastSqFitter1d<FittingFunc>::func1d(double value, int exponent) {
    return (*_funcArray[exponent])(value);
}


}}}}
