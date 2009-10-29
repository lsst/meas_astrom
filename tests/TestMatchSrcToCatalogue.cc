
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

//Create a trivial wcs solution
img::Wcs trivialWcs()
{
    img::PointD crval(0,0);
    img::PointD crpix(0,0);
    
    Eigen::Matrix2d CD(2,2);
    CD(0,0) = CD(1,1) = 1/3600.;
    CD(0,1) = CD(1,0) = 0;
    img::Wcs wcs(crval, crpix, CD);
    return wcs;
}


//Debugging function
void printSourceMatch(det::SourceMatch s) {

    boost::format fmt("(%.1f, %.1f) -- (%.7f, %.7f) -- (%.7f, %.7f) -- %.2f\n");
    fmt % (s.first)->getXAstrom();
    fmt % (s.first)->getYAstrom();
    fmt % (s.first)->getRa();
    fmt % (s.first)->getDec();
    fmt % (s.second)->getRa();
    fmt % (s.second)->getDec();
    fmt %  s.distance;

    std::string text = boost::str(fmt);
    std::cout << text;
}

//Debugging function
void printSourceSet(det::SourceSet s) {
    for(unsigned int i=0; i< s.size(); ++i) {
        boost::format fmt("%i: (%.1f, %.1f) -- (%.7f, %.7f)\n");
        fmt % i;
        fmt % (s[i])->getXAstrom();
        fmt % (s[i])->getYAstrom();
        fmt % (s[i])->getRa();
        fmt % (s[i])->getDec();
        
        std::string text = boost::str(fmt);
        std::cout << text;
    }
}


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

    img::Wcs wcs = trivialWcs();
    
    lsst::meas::astrom::sip::MatchSrcToCatalogue match(cat, data, wcs, 1);
    std::vector<det::SourceMatch> sm = match.getMatches();

    //Debugging code: Print out the values of the matched variables
    std::cout << sm.size() << std::endl;
    BOOST_CHECK_MESSAGE(sm.size() == 3, "Did not find 3 matches as expected");
        
    for(unsigned int i=0; i<sm.size(); ++i) {
        det::SourceMatch s = sm[i];
        BOOST_CHECK_CLOSE((s.first)->getRa(), (s.second)->getRa(), .01);
    }

}


//Test that every data point matches at most one point in the catalogue
BOOST_AUTO_TEST_CASE( OneToManyTest)
{
    
    int size=14;
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
    img::Wcs wcs = trivialWcs();

    //Change a data point so it lies close another data point.
    //data[0] is at x=1, now data[1] is at x=1.5. Both of which are close
    //to cat[0].
    img::PointD val = wcs.raDecToXY(cat[0]->getRa(), cat[0]->getDec());
    data[1]->setXAstrom(val[0]+.5);
    data[1]->setYAstrom(val[1]+.5);

    lsst::meas::astrom::sip::MatchSrcToCatalogue match(cat, data, wcs, 1);
    std::vector<det::SourceMatch> sm = match.getMatches();
    
    std::cout << sm.size() << std::endl;
    BOOST_CHECK(sm.size() == 7);
}



//Test that every catalogue entry matches at most one entry in the data set
BOOST_AUTO_TEST_CASE( ManyToOneTest)
{
    
    int size=14;
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
    img::Wcs wcs = trivialWcs();


    cat[12]->setRa(  (cat[1])->getRa() +1e-5);
    cat[12]->setDec( (cat[1])->getDec()+1e-5);

    lsst::meas::astrom::sip::MatchSrcToCatalogue match(cat, data, wcs, 1);
    std::vector<det::SourceMatch> sm = match.getMatches();

    std::cout << sm.size() << std::endl;
    BOOST_CHECK(sm.size() == 7);
}
