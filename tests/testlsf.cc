#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testlsf
 
//The boost unit test header
#include "boost/test/unit_test.hpp"


using namespace std;


#include <vector>
#include <cmath>
#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"
#include "lsst/afw/math/FunctionLibrary.h"

#include "lsst/meas/astrom/sip/LeastSqFitter2d.h"

namespace sip = lsst::meas::astrom::sip;
namespace math = lsst::afw::math;


BOOST_AUTO_TEST_CASE(fitLine)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    for(int i=0; i<7; ++i) {
        x.push_back( (double) i);
        y.push_back( (double) i);
        s.push_back( (double) 1);
        z.push_back((double) 2*i);

    }
    
    int order=2;
    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(x, y, z, s, order);
    
    Eigen::MatrixXd par = lsf.getParams();
    
    BOOST_CHECK(par.rows() == order);
    BOOST_CHECK(par.cols() == order);
    BOOST_CHECK_CLOSE( par(0,0)+1, 1., .001); //Check within .1% of zero
    BOOST_CHECK_CLOSE( par(0,1), 1., .001);   //Check is 1
    BOOST_CHECK_CLOSE( par(1,0), 1., .001);   //Check is 1
    BOOST_CHECK_CLOSE( par(1,1)+1, 1., .001); //Check within .1% of zero    
    
}



#include <cstdio>

BOOST_AUTO_TEST_CASE(fitLinearSurfaceX)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    int nData = 7;
    for(int i=0; i<nData*nData; ++i) {
        x.push_back( (double) (i % nData));
        y.push_back( (double) ((i- (int)x[i])/ nData));
        z.push_back(1 + 2*x[i]); 
        s.push_back(1);
        //printf("%.0f %.0f --> %.0f\n", x[i], y[i], z[i]);
    }        
    
    int order=2;
    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(x, y, z, s, order);
    
    Eigen::MatrixXd par = lsf.getParams();
    
    BOOST_CHECK(par.rows() == order);
    BOOST_CHECK(par.cols() == order);

    BOOST_CHECK_CLOSE( par(0,0), 1., .001); 
    BOOST_CHECK_CLOSE( par(1,0), 2., .001);     //linear x term
    BOOST_CHECK_CLOSE( par(0,1)+1, 1., .001); //Check within .1% of zero    
    BOOST_CHECK_CLOSE( par(1,1)+1, 1., .001); //Check within .1% of zero
        
}


BOOST_AUTO_TEST_CASE(fitLinearSurfaceY)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    int nData = 7;
    for(int i=0; i<nData*nData; ++i) {
        x.push_back( (double) (i % nData));
        y.push_back( (double) ((i- (int)x[i])/ nData));
        z.push_back(1 + 2*y[i]); 
        s.push_back(1);
        //printf("%.0f %.0f --> %.0f\n", x[i], y[i], z[i]);
    }        
    
    int order=2;
    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(x, y, z, s, order);
    
    Eigen::MatrixXd par = lsf.getParams();
    
    BOOST_CHECK(par.rows() == order);
    BOOST_CHECK(par.cols() == order);

    BOOST_CHECK_CLOSE( par(0,0), 1., .001); 
    BOOST_CHECK_CLOSE( par(0,1), 2., .001);     //linear y term
    BOOST_CHECK_CLOSE( par(1,0)+1, 1., .001); //Check within .1% of zero    
    BOOST_CHECK_CLOSE( par(1,1)+1, 1., .001); //Check within .1% of zero
        
}



BOOST_AUTO_TEST_CASE(fitLinearSurfaceXY)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    int nData = 7;
    for(int i=0; i<nData*nData; ++i) {
        x.push_back( (double) (i % nData));
        y.push_back( (double) ((i- (int)x[i])/ nData));
        z.push_back(1 + 2*x[i] + 3*y[i]); 
        s.push_back(1);
        //printf("%.0f %.0f --> %.0f\n", x[i], y[i], z[i]);
    }        
    
    int order=2;
    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(x, y, z, s, order);
    
    Eigen::MatrixXd par = lsf.getParams();
    
    BOOST_CHECK(par.rows() == order);
    BOOST_CHECK(par.cols() == order);

    BOOST_CHECK_CLOSE( par(0,0), 1., .001); 
    BOOST_CHECK_CLOSE( par(0,1), 3., .001);     //linear y term
    BOOST_CHECK_CLOSE( par(1,0), 2., .001);     //linear x term    
    BOOST_CHECK_CLOSE( par(1,1)+1, 1., .001); //Check within .1% of zero
        
}


BOOST_AUTO_TEST_CASE(fitQuadraticX)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    int nData = 3;
    for(int i=0; i<nData*nData; ++i) {
        x.push_back( (double) (i % nData));
        y.push_back( (double) ((i- (int)x[i])/ nData));
        z.push_back(1 + 2*x[i] + 3*pow(x[i], 2)); 
        s.push_back(1);
        //printf("%.0f %.0f --> %.0f\n", x[i], y[i], z[i]);
    }        
    
    int order=3;
    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(x, y, z, s, order);
    
    Eigen::MatrixXd par = lsf.getParams();
    
    BOOST_CHECK(par.rows() == order);
    BOOST_CHECK(par.cols() == order);

    BOOST_CHECK_CLOSE( par(0,0), 1., .001); 
    BOOST_CHECK_CLOSE( par(1,0), 2., .001);     //linear x term    
    BOOST_CHECK_CLOSE( par(2,0), 3., .001);     //quadratic x term    

    BOOST_CHECK_CLOSE( par(0,1)+1, 1., .001); //Check within .1% of zero        
    BOOST_CHECK_CLOSE( par(0,2)+1, 1., .001); //Check within .1% of zero       
     
    BOOST_CHECK_CLOSE( par(1,1)+1, 1., .001); //Check within .1% of zero        
    BOOST_CHECK_CLOSE( par(1,2)+1, 1., .001); //Check within .1% of zero        
    
    BOOST_CHECK_CLOSE( par(2,1)+1, 1., .001); //Check within .1% of zero        
    BOOST_CHECK_CLOSE( par(2,2)+1, 1., .001); //Check within .1% of zero        
}


BOOST_AUTO_TEST_CASE(fitQuadraticXY)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    int nData = 3;
    for(int i=0; i<nData*nData; ++i) {
        x.push_back( (double) (i % nData));
        y.push_back( (double) ((i- (int)x[i])/ nData));
        z.push_back(1 + 2*y[i] + 3*pow(y[i], 2) + 4*x[i] + 5*x[i]*y[i] + 6*pow(x[i], 2)); 
        s.push_back(1);
        //printf("%.0f %.0f --> %.0f\n", x[i], y[i], z[i]);
    }        
    
    int order=3;
    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(x, y, z, s, order);
    
    Eigen::MatrixXd par = lsf.getParams();
    
    BOOST_CHECK(par.rows() == order);
    BOOST_CHECK(par.cols() == order);



    BOOST_CHECK_CLOSE( par(0,0), 1., .001); 
    BOOST_CHECK_CLOSE( par(0,1), 2., .001); 
    BOOST_CHECK_CLOSE( par(0,2), 3., .001); 
    
    BOOST_CHECK_CLOSE( par(1,0), 4., .001); 
    BOOST_CHECK_CLOSE( par(1,1), 5., .001); 
    
    BOOST_CHECK_CLOSE( par(2,0), 6., .001); 
}
