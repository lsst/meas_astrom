
using namespace std;


#include <vector>
#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"
#include "Eigen/SVD"
#include "lsst/afw/math/FunctionLibrary.h"

namespace lsst { namespace meas { namespace astrom { namespace sip {

template <class FittingFunc>class LeastSqFitter2d
{
public:
    LeastSqFitter2d(const vector<double> x, const vector<double> y, const vector<double> z, const vector<double> s, \
                    int order);

    Eigen::VectorXd getParams();
    Eigen::VectorXd getErrors();

    void testInit() {
        for(int i=0; i< _order; ++i) {
            double val = (*_funcArray[i])(2.);
            cout << i << " " << val << endl;
        }
    }

private:
    void init();
             
    void calculateA();
    void calculateBeta();
    double func2d(double x, double y, int term);
    double func1d(double value, int exponent);

    vector<double> _x, _y, _z, _s;
    int _order;   //Degree of polynomial to fit, e.g 4=> cubic
    int _nPar;    //Number of parameters in fitting eqn, e.g x^2, xy, ^y^2, x^3, x^2y...
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

///Fit a 2d polynomial to a set of data points z(x, y)
///
///\param ReturnT e.g float or double
///\param FittingFunc  The type of function to fit in each dimension. This function should extend the base
///       class of lsst::afw::math::Function1. Note that this is a 1d function
///\param x vector of x positions of data
///\param y vector of y positions of data
///\param z Value of data for a given x,y. z = z_i = z_i(x_i, y_i)
///\param order Order of 2d function to fit
template<class FittingFunc> LeastSqFitter2d<FittingFunc>::LeastSqFitter2d(const vector<double> x, const vector<double> y, const vector<double> z, const vector<double> s, \
                int order) :
    _x(x), _y(y), _z(z), _order(order), _nPar(0), _A(order, order), _beta(order){
    
    //_nPar, the number of terms to fix (x^2, xy, y^2 etc.) is \Sigma^(order+1) 1
    _nPar=0;
    for(int i=0; i<order; ++i) {
        _nPar+= i+1;
    }
    
    _nData = _x.size();
    //@TODO check _y.size and _z.size == _nData
    
    init();
    calculateBeta();
    calculateA();
    _A.svd().solve(_beta, &_par);
}
        

template<class FittingFunc> Eigen::VectorXd LeastSqFitter2d<FittingFunc>::getParams() {
    return _par;
}


template<class FittingFunc> Eigen::VectorXd LeastSqFitter2d<FittingFunc>::getErrors() {
    Eigen::MatrixXd V = _A.svd().matrixV();
    Eigen::VectorXd w = _A.svd().singularValues();

    Eigen::VectorXd err(_nPar);
    for(int i=0; i<_nPar; ++i) {
        err(i)=0;
        for(int j=0; j<_nPar; ++j) {
            double val =  V(i,j)/w(j);
            err(i) += val*val;
        }
    }
    return err;
}

    
template<class FittingFunc> void LeastSqFitter2d<FittingFunc>::init() {
    //Initialise the array of functions. _funcArray[i] is a object of type math::Function1 of order i
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
    

template<class FittingFunc> void LeastSqFitter2d<FittingFunc>::calculateA() {
    for(int i=0; i< _nPar; ++i) {
        for(int j=0; j< _nPar; ++j) {
            double val=0;
            for(int k=0; k< _nData; ++k) {
                val = 1/(_s[k]* _s[k]);
                val *= func2d(_x[k], _y[k], i) * func2d(_x[k], _y[k], j);
            }
            _A(i,j) = val;
        }
    }        
}

    
template<class FittingFunc> void LeastSqFitter2d<FittingFunc>::calculateBeta() {
    for(int i=0; i< _nPar; ++i) {
        _beta(i) = 0;
        for(int j=0; j< _nData; ++j) {
            _beta(i) += _z[j]*func2d(_x[j], _y[j], i)/ (_s[i]*_s[i]);
        }
    }
}
 
 
 
//The ith term in the fitting polynomial is of the form x^a * y^b. This function figures
//out the value of a and b, then calculates the value of the ith term at the given x and y
template<class FittingFunc> double LeastSqFitter2d<FittingFunc>::func2d(double x, double y, int term) {
    int yexp = term % _order;  //y exponent
    int xexp = (term - yexp) / _order;  //x exponent

    double xcomp = func1d(x, xexp);  //x component of polynomial
    double ycomp = func1d(y, yexp);  //y component

    return xcomp*ycomp;
}

template<class FittingFunc> double LeastSqFitter2d<FittingFunc>::func1d(double value, int exponent) {
    return (*_funcArray[exponent])(value);
}


}}}}
