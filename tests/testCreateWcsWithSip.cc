#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testlsf1d
 
//The boost unit test header
#include "boost/test/unit_test.hpp"


using namespace std;


#include <vector>
#include <cmath>
#include <cstdio>
#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"
#include "lsst/afw/math/FunctionLibrary.h"

#include "lsst/meas/astrom/sip/createWcsSip.h"

namespace sip = lsst::meas::astrom::sip;
namespace math = lsst::afw::math;


BOOST_AUTO_TEST_CASE(reScale1)
{
    vector<double> x;

    int size=7;
    for(int i=0; i<size; ++i) {
        x.push_back( (double) i);
    }
    
    x = sip::scaleVector<double>(x, 10, 20);
        
    BOOST_CHECK_CLOSE( x[0], 10., .001); 
    BOOST_CHECK_CLOSE( x[size-1], 20., .001);  

}


BOOST_AUTO_TEST_CASE(reScale2)
{
    vector<double> x;

    int size=7;
    for(int i=0; i<size; ++i) {
        x.push_back( (double) i);
    }
    
    x = sip::scaleVector<double>(x, 10, 20, 3., 4.);
        
    BOOST_CHECK_CLOSE( x[0], -20., .001); 
    BOOST_CHECK_CLOSE( x[size-1], 40., .001);  

}

