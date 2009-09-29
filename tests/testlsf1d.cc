#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testlsf
 
//The boost unit test header
#include "boost/test/unit_test.hpp"


using namespace std;


#include <vector>
#include <cmath>
#include <cstdio>
#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"
#include "lsst/afw/math/FunctionLibrary.h"

#include "lsst/meas/astrom/sip/LeastSqFitter1d.h"

namespace sip = lsst::meas::astrom::sip;
namespace math = lsst::afw::math;


BOOST_AUTO_TEST_CASE(fitLine)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;

    for(int i=0; i<7; ++i) {
        x.push_back( (double) i);
        y.push_back( (double) 3+ 2*i);
        s.push_back( (double) 1);

    }
    
    int order=2;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);
    
    boost::shared_ptr<Eigen::VectorXd> parPtr = lsf.getParams();
    
    BOOST_CHECK(parPtr->size() == order);
    BOOST_CHECK_CLOSE( (*parPtr)(0), 3., .001); 
    BOOST_CHECK_CLOSE( (*parPtr)(1), 2., .001);  

    
}


BOOST_AUTO_TEST_CASE(fitQuadratic)
{
    vector<double> x;
    vector<double> y;
    vector<double> s;

    for(int i=0; i<7; ++i) {
        x.push_back( (double) i);
        y.push_back( (double) 4 + i*(3+ i*2));
        s.push_back( (double) 1);

    }
    
    int order=3;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);
    
    boost::shared_ptr<math::PolynomialFunction1<double> > f =lsf.getBestFitFunction();

    printf("%i == %i\n", (*f).getNParameters(), order);
    BOOST_CHECK( (int) (*f).getNParameters() == order);
    BOOST_CHECK_CLOSE( (*f).getParameter(0), 4., .001); 
    BOOST_CHECK_CLOSE( (*f).getParameter(1), 3., .001); 
    BOOST_CHECK_CLOSE( (*f).getParameter(2), 2., .001);         

    BOOST_CHECK_CLOSE( lsf.valueAt(1.), 9., .001);

    
}

/*
#include <vector>
#include <cmath>
#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"
#include "lsst/afw/math/FunctionLibrary.h"

#include "lsst/meas/astrom/sip/LeastSqFitter1d.h"

namespace sip = lsst::meas::astrom::sip;
namespace math = lsst::afw::math;


int main()
{

    vector<double> x;
    vector<double> y;
    vector<double> s;

    for(int i=0; i<7; ++i) {
        x.push_back( (double) i);
        y.push_back( (double) 2*i);
        s.push_back( (double) 1);
    }
    
    int order=2;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> >  lsf(x,y,s, order);
    
    boost::shared_ptr<math::PolynomialFunction1<double> > f =lsf.getBestFitFunction();
    std::cout << (*f)(5.) << std::endl;
    
    
}
    
*/
