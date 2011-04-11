/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testlsf1d

//The boost unit test header
#include "boost/test/unit_test.hpp"

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>
#include "boost/shared_ptr.hpp"
#include "Eigen/Core.h"

#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/image/TanWcs.h"

#include "lsst/meas/astrom/sip/CreateWcsWithSip.h"


using namespace std;
namespace except = lsst::pex::exceptions;
namespace pexLog = lsst::pex::logging;
namespace afwCoord = lsst::afw::coord;
namespace afwGeom = lsst::afw::geom;
namespace afwImg = lsst::afw::image;
namespace afwDet = lsst::afw::detection;
namespace det = lsst::afw::detection;
namespace math = lsst::afw::math;


afwImg::TanWcs::Ptr createWcsPtr(afwGeom::Point2D crval, afwGeom::Point2D crpix, Eigen::Matrix2d CD)
{
    afwImg::TanWcs::Ptr wcs = afwImg::TanWcs::Ptr(new  afwImg::TanWcs(crval, crpix, CD));
    return wcs;
}


afwImg::TanWcs::Ptr createDefaultWcsPtr()
{

    afwGeom::Point2D crval = afwGeom::Point2D(44., 45.);
    afwGeom::Point2D crpix = afwGeom::Point2D(0,0);   //(0,0) in lsst coords
    
    double arcsecPerPixel = 1/3600.;
    Eigen::Matrix2d CD;
    CD(0,0) = arcsecPerPixel;
    CD(0,1) = 0;
    CD(1,0) = 0;
    CD(1,1) = arcsecPerPixel;

    return createWcsPtr(crval, crpix, CD);
}



vector<afwDet::SourceMatch> generateSourceSet(afwImg::TanWcs::Ptr wcsPtr)
{
    int num = 4000;
    int step= 1000;
    vector<afwDet::SourceMatch> sourceMatchSet;

    for(int i=0; i<=num; i+=step)
    {
        for(int j=0; j<=num; j+=step)
        {
            afwDet::Source::Ptr src = afwDet::Source::Ptr(new afwDet::Source());
            afwDet::Source::Ptr cat = afwDet::Source::Ptr(new afwDet::Source());

            cat->setXAstrom(i);
            cat->setYAstrom(j);

            src->setXAstrom(i);
            src->setYAstrom(j);

            afwCoord::Fk5Coord c = wcsPtr->pixelToSky(i, j)->toFk5();
            
            cat->setRa(c.getRa(afwCoord::DEGREES));
            cat->setDec(c.getDec(afwCoord::DEGREES));
            
            double dist = hypot(src->getXAstrom()-cat->getXAstrom(), src->getYAstrom() - cat->getYAstrom());
            sourceMatchSet.push_back(afwDet::SourceMatch(cat, src, dist));
        }

    }

    return sourceMatchSet;


}


void checkResults(afwImg::Wcs::Ptr wcsPtr, afwImg::TanWcs::Ptr sipWcsPtr, vector<afwDet::SourceMatch> sourceMatchSet)
{
    for(unsigned int i=0; i< sourceMatchSet.size(); ++i)
    {
        double catX = sourceMatchSet[i].first->getXAstrom();
        double catY = sourceMatchSet[i].first->getYAstrom();    
        
        double srcX = sourceMatchSet[i].second->getXAstrom();
        double srcY = sourceMatchSet[i].second->getYAstrom();    
    
        afwCoord::Fk5Coord catCoo = wcsPtr->pixelToSky(catX, catY)->toFk5();
        afwCoord::Fk5Coord srcCoo = sipWcsPtr->pixelToSky(srcX, srcY)->toFk5();
        
        double catA = sourceMatchSet[i].first->getRa(); //catCoo.getRa(afwCoord::DEGREES);
        double catD = sourceMatchSet[i].first->getDec(); //catCoo.getDec(afwCoord::DEGREES);
        
        double srcA = srcCoo.getRa(afwCoord::DEGREES);
        double srcD = srcCoo.getDec(afwCoord::DEGREES);

        //Forward transform (units are degrees)    
        BOOST_CHECK_CLOSE(catA, srcA, 1e-4);         
        BOOST_CHECK_CLOSE(catD, srcD, 1e-4);         
        
        //Reverse tranform (units are pixels)
        afwGeom::Point2D sipxy = sipWcsPtr->skyToPixel(catA, catD);
        BOOST_CHECK_CLOSE(srcX+1, sipxy[0]+1, 1e-8);         
        BOOST_CHECK_CLOSE(srcY+1, sipxy[1]+1, 1e-8);         

    }

}


BOOST_AUTO_TEST_CASE(trivial)
{

    afwImg::TanWcs::Ptr wcsPtr = createDefaultWcsPtr();
    vector<afwDet::SourceMatch> sourceMatchSet = generateSourceSet(wcsPtr);
    
    //Add no distortion

    lsst::meas::astrom::sip::CreateWcsWithSip sipObject(sourceMatchSet, wcsPtr, 2);
    std::cout << "trivial" << std::endl;
    std::cout << sipObject.getNewWcs()->getFitsMetadata()->toString() << std::endl;
    checkResults(wcsPtr, sipObject.getNewWcs(), sourceMatchSet);
}


BOOST_AUTO_TEST_CASE(offset)
{

    afwImg::TanWcs::Ptr wcsPtr = createDefaultWcsPtr();
    vector<afwDet::SourceMatch> sourceMatchSet = generateSourceSet(wcsPtr);
    
    //Add an offset to each point
    for(unsigned int i = 0; i<sourceMatchSet.size(); ++i)
    {
        double x = sourceMatchSet[i].second->getXAstrom();
        double y = sourceMatchSet[i].second->getYAstrom();
        
        sourceMatchSet[i].second->setXAstrom(x+5);
        sourceMatchSet[i].second->setYAstrom(y+7);
    }

    lsst::meas::astrom::sip::CreateWcsWithSip sipObject(sourceMatchSet, wcsPtr, 2);
    std::cout << "offset" << std::endl;
    std::cout << sipObject.getNewWcs()->getFitsMetadata()->toString() << std::endl;
    checkResults(wcsPtr, sipObject.getNewWcs(), sourceMatchSet);
}



BOOST_AUTO_TEST_CASE(linearX)
{

    afwImg::TanWcs::Ptr wcsPtr = createDefaultWcsPtr();
    vector<afwDet::SourceMatch> sourceMatchSet = generateSourceSet(wcsPtr);
    
    //Add an offset to each point
    for(unsigned int i = 0; i<sourceMatchSet.size(); ++i)
    {
        double x = sourceMatchSet[i].second->getXAstrom();
        double y = sourceMatchSet[i].second->getYAstrom();
        
        sourceMatchSet[i].second->setXAstrom(2*x);
        sourceMatchSet[i].second->setYAstrom(y+7);
    }

    lsst::meas::astrom::sip::CreateWcsWithSip sipObject(sourceMatchSet, wcsPtr, 2);
    std::cout << "linearX" << std::endl;
    std::cout << sipObject.getNewWcs()->getFitsMetadata()->toString() << std::endl;
    checkResults(wcsPtr, sipObject.getNewWcs(), sourceMatchSet);
}


BOOST_AUTO_TEST_CASE(linearXY)
{

    afwImg::TanWcs::Ptr wcsPtr = createDefaultWcsPtr();
    vector<afwDet::SourceMatch> sourceMatchSet = generateSourceSet(wcsPtr);
    
    //Add an offset to each point
    for(unsigned int i = 0; i<sourceMatchSet.size(); ++i)
    {
        double x = sourceMatchSet[i].second->getXAstrom();
        double y = sourceMatchSet[i].second->getYAstrom();
        
        sourceMatchSet[i].second->setXAstrom(2*x);
        sourceMatchSet[i].second->setYAstrom(3*y);
    }

    lsst::meas::astrom::sip::CreateWcsWithSip sipObject(sourceMatchSet, wcsPtr, 2);
    std::cout << "linearXY" << std::endl;
    std::cout << sipObject.getNewWcs()->getFitsMetadata()->toString() << std::endl;
    checkResults(wcsPtr, sipObject.getNewWcs(), sourceMatchSet);
}


BOOST_AUTO_TEST_CASE(linearYX)
{

    afwImg::TanWcs::Ptr wcsPtr = createDefaultWcsPtr();
    vector<afwDet::SourceMatch> sourceMatchSet = generateSourceSet(wcsPtr);
    
    //Add an offset to each point
    for(unsigned int i = 0; i<sourceMatchSet.size(); ++i)
    {
        double x = sourceMatchSet[i].second->getXAstrom();
        double y = sourceMatchSet[i].second->getYAstrom();
        
        sourceMatchSet[i].second->setXAstrom(x + .2*y);
        sourceMatchSet[i].second->setYAstrom(y + .3*x);
    }

    lsst::meas::astrom::sip::CreateWcsWithSip sipObject(sourceMatchSet, wcsPtr, 2);
    std::cout << "linearYX" << std::endl;
    std::cout << sipObject.getNewWcs()->getFitsMetadata()->toString() << std::endl;
    checkResults(wcsPtr, sipObject.getNewWcs(), sourceMatchSet);
}



BOOST_AUTO_TEST_CASE(quadraticX)
{

    afwImg::TanWcs::Ptr wcsPtr = createDefaultWcsPtr();
    vector<afwDet::SourceMatch> sourceMatchSet = generateSourceSet(wcsPtr);
    
    //Add an offset to each point
    for(unsigned int i = 0; i<sourceMatchSet.size(); ++i)
    {
        double x = sourceMatchSet[i].second->getXAstrom();
        double y = sourceMatchSet[i].second->getYAstrom();
        
        sourceMatchSet[i].second->setXAstrom(x + 1e-5*x*x);
        sourceMatchSet[i].second->setYAstrom(y);
    }

    lsst::meas::astrom::sip::CreateWcsWithSip sipObject(sourceMatchSet, wcsPtr, 3);
    std::cout << "quadraticX" << std::endl;
    std::cout << sipObject.getNewWcs()->getFitsMetadata()->toString() << std::endl;
    checkResults(wcsPtr, sipObject.getNewWcs(), sourceMatchSet);
}



