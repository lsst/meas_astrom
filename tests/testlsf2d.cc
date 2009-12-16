#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testlsf2d
 
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
        y.push_back( (double) 0.);
        s.push_back( (double) 1.);
        z.push_back((double) 2*i);

    }
    
    int order=2;
    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(x, y, z, s, order);
    
    Eigen::MatrixXd par = lsf.getParams();
    
    BOOST_CHECK(par.rows() == order);
    BOOST_CHECK(par.cols() == order);
    BOOST_CHECK_CLOSE( par(0,0)+1, 1., .001); //Check within .1% of zero
    BOOST_CHECK_CLOSE( par(0,1), 0., .001);   //Check is 1
    BOOST_CHECK_CLOSE( par(1,0), 2., .001);   //Check is 1
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
    BOOST_CHECK_CLOSE( par(0,1)+1, 1., .001); //Check within .1% of zero    
    BOOST_CHECK_CLOSE( par(1,0), 2., .001);     //linear x term
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


BOOST_AUTO_TEST_CASE(fitLinearXSurface2)
{

    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> s;

    x.push_back(599.59899999999993);
    x.push_back(1172.7726619709097);
    x.push_back(512.51199999999994);
    x.push_back(1083.6078436082901);
        
    y.push_back(512.0);
    y.push_back(539.77214401699996);
    y.push_back(541.0);
    y.push_back(562.09371856300004);

    z.push_back(0.5989999999999327);
    z.push_back(1.1716010609097793);
    z.push_back(0.51199999999994361);
    z.push_back(1.0825253182899814);
    
    s.push_back(.1);
    s.push_back(.1);
    s.push_back(.1);
    s.push_back(.1);
    
    int order=2;
    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(x, y, z, s, order);
    Eigen::MatrixXd par = lsf.getParams();

    lsf.printParams();
    for(unsigned int i=0; i<x.size(); ++i)
    {
        printf("%.3f %.3f %.3f %.3f\n", x[i], y[i], z[i], lsf.valueAt(x[i], y[i]));
        BOOST_CHECK_CLOSE(z[i], lsf.valueAt(x[i], y[i]), .001);
    }
}



BOOST_AUTO_TEST_CASE(fitChebyshevX1)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    int nData = 3;
    for(int i=0; i<nData*nData; ++i) {
        x.push_back( (double) (i % nData));
        y.push_back( (double) ((i- (int)x[i])/ nData));
        z.push_back(1 + 2*x[i]); //T0 + 2T1
        s.push_back(1);
    }        
    
    int order=3;
    sip::LeastSqFitter2d<math::Chebyshev1Function1<double> > lsf(x, y, z, s, order);
    Eigen::MatrixXd par = lsf.getParams();

    BOOST_CHECK_CLOSE( par(0,0), 1., .001); 
    BOOST_CHECK_CLOSE( par(1,0), 2., .001); 
    
    BOOST_CHECK_CLOSE( par(0,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(0,2)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(1,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(1,2)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(2,0)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(2,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(2,2)+1 , 1., .001); 

}


BOOST_AUTO_TEST_CASE(fitChebyshevX2)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    //Create a Cheby polynomial
    std::vector<double> cp;
    cp.push_back(1);
    cp.push_back(2);
    cp.push_back(3);
    math::Chebyshev1Function1<double> cheby(cp);
    
    
    int nData = 7;
    for(int i=0; i<nData*nData; ++i) {
        x.push_back( (double) (i % nData));
        y.push_back( (double) ((i- (int)x[i])/ nData));
        z.push_back(cheby(x[i])); 
        s.push_back(1);
    }        
    
    int order=3;
    sip::LeastSqFitter2d<math::Chebyshev1Function1<double> > lsf(x, y, z, s, order);
    Eigen::MatrixXd par = lsf.getParams();
    std::cout << par << std::endl;

    BOOST_CHECK_CLOSE( par(0,0), 1., .001); 
    BOOST_CHECK_CLOSE( par(1,0), 2., .001); 
    BOOST_CHECK_CLOSE( par(2,0), 3., .001); 
        
    BOOST_CHECK_CLOSE( par(0,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(0,2)+1 , 1., .001); 
    
    BOOST_CHECK_CLOSE( par(1,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(1,2)+1 , 1., .001); 

    BOOST_CHECK_CLOSE( par(2,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(2,2)+1 , 1., .001); 

}



BOOST_AUTO_TEST_CASE(fitChebyshevX3)
{
    //A test case for a specific problem I've run into and don't understand
    
    vector<double> x;
    vector<double> y;
    vector<double> s;
    vector<double> z;

    x.push_back(-1);
    y.push_back(0);
    z.push_back(-1e-3);
    s.push_back(1);
    
    x.push_back(-.75);
    y.push_back(0);
    z.push_back(-1.25e-3);
    s.push_back(1);

    x.push_back(-.5);
    y.push_back(0);
    z.push_back(-1.5e-3);
    s.push_back(1);

    x.push_back(-.25);
    y.push_back(0);
    z.push_back(-1.75e-3);
    s.push_back(1);

    x.push_back(0.);
    y.push_back(0);
    z.push_back(-2e-3);
    s.push_back(1);

    x.push_back(.25);
    y.push_back(0);
    z.push_back(-2.25e-3);
    s.push_back(1);

    x.push_back(.5);
    y.push_back(0);
    z.push_back(-2.5e-3);
    s.push_back(1);

    x.push_back(.75);
    y.push_back(0);
    z.push_back(-2.75e-3);
    s.push_back(1);

    x.push_back(1.);
    y.push_back(0);
    z.push_back(-3e-3);
    s.push_back(1);
    
    int order=3;
    sip::LeastSqFitter2d<math::Chebyshev1Function1<double> > lsf(x, y, z, s, order);
    Eigen::MatrixXd par = lsf.getParams();
    std::cout << par << std::endl;

    BOOST_CHECK_CLOSE( par(0,0), -2e-3, .001); 
    BOOST_CHECK_CLOSE( par(1,0), -1e-3, .001); 
        
    BOOST_CHECK_CLOSE( par(0,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(0,2)+1 , 1., .001); 
    
    BOOST_CHECK_CLOSE( par(1,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(1,2)+1 , 1., .001); 

    BOOST_CHECK_CLOSE( par(2,0)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(2,1)+1 , 1., .001); 
    BOOST_CHECK_CLOSE( par(2,2)+1 , 1., .001); 

}



