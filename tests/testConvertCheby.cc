
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestMatchSrcToCatalogue

//The boost unit test header
#include "boost/test/unit_test.hpp"

#include <iostream>
#include <cmath>
#include <vector>

#include "boost/format.hpp"
#include "Eigen/Core.h"

#include "lsst/meas/astrom/sip/CreateWcsWithSip.h"

namespace sip = lsst::meas::astrom::sip;

//This test doesn't work now that I've made a class out of createWcsWithSip

#if 0
///The simplest test I can think of. Create a series of points in a diagonal line starting
/// at (1,1) in pixel space, and a wcs that converts 1 pixel length to 1 arcsecond. Then 
/// create a catalogue which half over laps with input data and check for matches.    
BOOST_AUTO_TEST_CASE( convertChebyTest00)
{
    Eigen::MatrixXd cheby(Eigen::MatrixXd::Zero(2,2));
    
    cheby(0,0) = 1;
    
    Eigen::MatrixXd test = sip::convertChebyToSip(cheby);
    
    BOOST_CHECK_CLOSE(test(0,0), 1., 1e-4);
    BOOST_CHECK_CLOSE(test(0,1)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,1)+1, 1., 1e-4);
}


BOOST_AUTO_TEST_CASE( convertChebyTest01)
{
    Eigen::MatrixXd cheby(Eigen::MatrixXd::Zero(2,2));
    
    cheby(0,1) = 1;
    
    Eigen::MatrixXd test = sip::convertChebyToSip(cheby);
    
    BOOST_CHECK_CLOSE(test(0,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(0,1), 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,1)+1, 1., 1e-4);
}


BOOST_AUTO_TEST_CASE( convertChebyTest10)
{
    Eigen::MatrixXd cheby(Eigen::MatrixXd::Zero(2,2));
    
    cheby(1,0) = 1;
    
    Eigen::MatrixXd test = sip::convertChebyToSip(cheby);
    
    BOOST_CHECK_CLOSE(test(0,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(0,1)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,0), 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,1)+1, 1., 1e-4);
}


BOOST_AUTO_TEST_CASE( convertChebyTest11)
{
    Eigen::MatrixXd cheby(Eigen::MatrixXd::Zero(2,2));
    
    cheby(1,1) = 1;
    
    Eigen::MatrixXd test = sip::convertChebyToSip(cheby);
    
    BOOST_CHECK_CLOSE(test(0,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(0,1)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,1), 1., 1e-4);
}


BOOST_AUTO_TEST_CASE( convertChebyTest00_11)
{
    Eigen::MatrixXd cheby(Eigen::MatrixXd::Zero(2,2));

    cheby(0,0) = 1;    
    cheby(1,1) = 1;
    
    Eigen::MatrixXd test = sip::convertChebyToSip(cheby);
    
    BOOST_CHECK_CLOSE(test(0,0), 1., 1e-4);
    BOOST_CHECK_CLOSE(test(0,1)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,1), 1., 1e-4);
}


//A harder example. T2(x)*T3(y) = 
//8x^2y^3 + 6x^2y - 4y^2 - 3y
BOOST_AUTO_TEST_CASE( convertChebyTest23)
{
    Eigen::MatrixXd cheby(Eigen::MatrixXd::Zero(4,4));

    cheby(2,3) = 1;    

    Eigen::MatrixXd test = sip::convertChebyToSip(cheby);
    
    //Non zero components
    BOOST_CHECK_CLOSE(test(0,1), 3., 1e-4);    
    BOOST_CHECK_CLOSE(test(0,3), -4., 1e-4);
    BOOST_CHECK_CLOSE(test(2,1), -6., 1e-4);
    BOOST_CHECK_CLOSE(test(2,3), 8., 1e-4);
    
    BOOST_CHECK_CLOSE(test(0,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(0,2)+1, 1., 1e-4);    
    
    BOOST_CHECK_CLOSE(test(1,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(1,1)+1, 1., 1e-4);    
    BOOST_CHECK_CLOSE(test(1,2)+1, 1., 1e-4);    
    BOOST_CHECK_CLOSE(test(1,3)+1, 1., 1e-4);    

    BOOST_CHECK_CLOSE(test(2,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(2,2)+1, 1., 1e-4);    

    BOOST_CHECK_CLOSE(test(3,0)+1, 1., 1e-4);
    BOOST_CHECK_CLOSE(test(3,1)+1, 1., 1e-4);    
    BOOST_CHECK_CLOSE(test(3,2)+1, 1., 1e-4);    
    BOOST_CHECK_CLOSE(test(3,3)+1, 1., 1e-4);    

}
#endif
