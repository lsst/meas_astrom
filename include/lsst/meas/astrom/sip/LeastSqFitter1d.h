
using namespace std;


#include <vector>
#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"
#include "Eigen/SVD"
#include "lsst/afw/math/FunctionLibrary.h"

namespace lsst { namespace meas { namespace astrom { namespace sip {

template <class FittingFunc>class LeastSqFitter1d
{
public:
    LeastSqFitter1d(const vector<double> &x, const vector<double> &y, const vector<double> &s,  int order);

    boost::shared_ptr<Eigen::VectorXd> getParams();
    boost::shared_ptr<FittingFunc> getBestFitFunction();
    double valueAt(double x);
    
private:
    void initFunctions();
             
    void calculateA();
    void calculateBeta();
    double func1d(double value, int exponent);
    

    vector<double> _x, _y, _z, _s;
    int _order;   //Degree of polynomial to fit, e.g 4=> cubic
    int _nData;   //Number of data points, == _x.size()
    
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
template<class FittingFunc> LeastSqFitter1d<FittingFunc>::LeastSqFitter1d(const vector<double> &x, const vector<double> &y, const vector<double> &s, int order) :
    _x(x), _y(y), _s(s), _order(order), _A(1,1), _beta(1) {
    
    
    _nData = _x.size();
    //@TODO check _y.size and _z.size == _nData
    
    initFunctions();
    calculateBeta();
    calculateA();

    _par = Eigen::VectorXd(_order);
    _A.svd().solve(_beta, &_par);
}
        


///Return the best fit paramets as an Eigen::Matrix
template<class FittingFunc> boost::shared_ptr<Eigen::VectorXd> LeastSqFitter1d<FittingFunc>::getParams() {

    boost::shared_ptr<Eigen::VectorXd> vPtr(new Eigen::VectorXd(_order)); 
    for(int i=0; i< _order; ++i) {
        (*vPtr)(i) = _par(i);
    }
    return vPtr;
}


///Return the best fit polynomial
template<class FittingFunc> 
boost::shared_ptr<FittingFunc> LeastSqFitter1d<FittingFunc>::getBestFitFunction() {

    //FittingFunc and LeastSqFitter disagree on the definition of order of a function.
    //LSF says that a linear function is order 2 (two coefficients), FF says only 1
    boost::shared_ptr<FittingFunc> funcPtr(new FittingFunc(_order-1) );

    for(int i=0; i< _order; ++i) {
        funcPtr->setParameter(i, _par(i));
    }
    return funcPtr;
}


///Calculate the value of the function at a given point
template<class FittingFunc> double LeastSqFitter1d<FittingFunc>::valueAt(double x) {
        static boost::shared_ptr<FittingFunc> f = getBestFitFunction();
        
        return (*f)(x);
}

/// Initialise the array of functions. _funcArray[i] is a object of type math::Function1 of order 
///_norder
template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::initFunctions() {
    _funcArray.reserve(_order);
    
    vector<double> coeff;
    coeff.reserve( _order);
    
    coeff.push_back(1);
    for( int i=0; i< _order; ++i) {
        boost::shared_ptr<FittingFunc> p(new FittingFunc(coeff));
        _funcArray.push_back(p);
        
        coeff[i] = 0;
        coeff.push_back(1);  //coeff now looks like [0,0,...,0,1]
    }
}
    

template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::calculateA() {

    assert(_order >= 0);
    _A = Eigen::MatrixXd(_order, _order);

    int i, j;
    for(i=0; i< _order; ++i) {
        for(j=0; j< _order; ++j) {
            double val=0;
            for(int k=0; k< _nData; ++k) {
                val += func1d(_x[k], i) * func1d(_x[k], j)/( _s[k]*_s[k]);
            }
            _A(i,j) = val;
        }
        
    }        
}

    
template<class FittingFunc> void LeastSqFitter1d<FittingFunc>::calculateBeta() {

    assert(_order >= 0);
    _beta = Eigen::VectorXd(_order);

    double val;
    int j;
    for(int i=0; i< _order; ++i) {
        _beta(i) = 0;
        for(j=0; j< _nData; ++j) {
            val = _y[j]*func1d(_x[j], i)/ (_s[i]*_s[i]);
            _beta(i) += val;
        }
    }
}
 
 

template<class FittingFunc> double LeastSqFitter1d<FittingFunc>::func1d(double value, int exponent) {
    return (*_funcArray[exponent])(value);
}


}}}}
