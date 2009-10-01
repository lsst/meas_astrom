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

/*
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
    
    Eigen::VectorXd par = lsf.getParams();
    
    BOOST_CHECK(par.size() == order);
    BOOST_CHECK_CLOSE( par(0), 3., .001); 
    BOOST_CHECK_CLOSE( par(1), 2., .001);  

    
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
    
    math::PolynomialFunction1<double> f =lsf.getBestFitFunction();

    printf("%i == %i\n", f.getNParameters(), order);
    BOOST_CHECK( (int) f.getNParameters() == order);
    BOOST_CHECK_CLOSE( f.getParameter(0), 4., .001); 
    BOOST_CHECK_CLOSE( f.getParameter(1), 3., .001); 
    BOOST_CHECK_CLOSE( f.getParameter(2), 2., .001);         

    BOOST_CHECK_CLOSE( lsf.valueAt(1.), 9., .001);

    
}



//Check a more extreme case where there's a big difference between the x and y values
BOOST_AUTO_TEST_CASE(fitLinear2)
{

    vector<double> x;
    vector<double> y;
    vector<double> s;

    x.push_back(631.920343622);
    x.push_back(1000.02009955);
    x.push_back(1205.63487689);
    x.push_back(1427.86231185);
    
    y.push_back(0.622249643993);
    y.push_back(1.01887634344);
    y.push_back(1.1982950985);
    y.push_back(1.4294170431);
    
    s.push_back(1);
    s.push_back(1);
    s.push_back(1);
    s.push_back(1);
    
    int order=2;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);
    Eigen::VectorXd par = lsf.getParams();
    

}
*/

BOOST_AUTO_TEST_CASE(fitLinear3)
{

    vector<double> x;
    vector<double> y;
    vector<double> s;

    x.push_back(689.301136505);
    x.push_back(1112.8573687);
    x.push_back(1386.67168477);
    
    y.push_back(0.66911456573);
    y.push_back(1.1147439759);
    y.push_back(1.39597284177);
    
    s.push_back(1);
    s.push_back(1);
    s.push_back(1);
    
    int order=2;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);
    Eigen::VectorXd par = lsf.getParams();
    
    BOOST_CHECK_CLOSE(par(0), -0.0488399, .001);
    BOOST_CHECK_CLOSE(par(1), 0.00104313, .001);

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
