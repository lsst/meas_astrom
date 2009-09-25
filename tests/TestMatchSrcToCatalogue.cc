
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestMatchSrcToCatalogue

//The boost unit test header
#include "boost/test/unit_test.hpp"

#include <iostream>
#include <cmath>
#include <vector>

#include "boost/format.hpp"
#include "Eigen/Core.h"

#include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"

namespace det = lsst::afw::detection;
namespace img = lsst::afw::image;
    
///The simplest test I can think of. Create a series of points in a diagonal line starting at (1,1) in pixel
///space, and a wcs that converts 1 pixel length to 1 arcsecond. Then create a catalogue which half over laps
///with input data and check for matches.    
BOOST_AUTO_TEST_CASE( OutputTest)
{
    
    int size=7;
    det::SourceSet data;
    det::SourceSet cat;

    //Populate a set of data points and catalogue points
    //to make partially overlapping sets
    for(int i=1; i< size; ++i) {
        det::Source::Ptr s1(new det::Source());
        s1->setXAstrom(i);
        s1->setYAstrom(i);
        data.push_back(s1);

        det::Source::Ptr s2(new det::Source());
        s2->setRa(2*i/3600.);
        s2->setDec(2*i/3600.);
        cat.push_back(s2);
    }

    //Create a trivial wcs solution
    img::PointD crval(0,0);
    img::PointD crpix(0,0);
    
    Eigen::Matrix2d CD(2,2);
    CD(0,0) = CD(1,1) = 1/3600.;
    CD(0,1) = CD(1,0) = 0;
    img::Wcs wcs(crval, crpix, CD);

    lsst::meas::astrom::sip::MatchSrcToCatalogue match(cat, data, wcs, 1);
    std::vector<det::SourceMatch> sm = match.getMatches();

    //Debugging code: Print out the values of the matched variables
    std::cout << sm.size() << std::endl;
    BOOST_CHECK_MESSAGE(sm.size() == 3, "Did not find 3 matches as expected");
        
    for(unsigned int i=0; i<sm.size(); ++i) {
        det::SourceMatch s = sm[i];
        
        boost::format fmt("(%.1f, %.1f) -- (%.7f, %.7f) -- (%.7f, %.7f) -- %.2f\n");
        fmt % (boost::tuples::get<0>(s))->getXAstrom();
        fmt % (boost::tuples::get<0>(s))->getYAstrom();
        fmt % (boost::tuples::get<0>(s))->getRa();
        fmt % (boost::tuples::get<0>(s))->getDec();
        fmt % (boost::tuples::get<1>(s))->getRa();
        fmt % (boost::tuples::get<1>(s))->getDec();
        fmt % (boost::tuples::get<2>(s));
        
        std::string text = boost::str(fmt);
        std::cout << text;
        
        BOOST_CHECK_CLOSE((boost::tuples::get<0>(s))->getRa(), (boost::tuples::get<1>(s))->getRa(), .01);
    }
        
        

    
    
}
