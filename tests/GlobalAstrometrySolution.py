#!/usr/bin/env python
import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.afw.image as afwImage
import lsst.meas.astrom.net as net
import lsst.utils.tests as utilsTests
import lsst.afw.image.imageLib as img
import lsst.afw.detection.detectionLib as detect
try:
    type(verbose)
except NameError:
    verbose = 0

dataDir = eups.productDir("astrometry_net_data")
if not dataDir:
    raise RuntimeError("Must set up astrometry_net_data to run these tests")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def loadXYFromFile(filename):
    """Load a list of positions from a file"""
    f= open(filename)
    
    s1=detect.SourceSet()
    i=0
    for line in f:
        #Split the row into an array
        line = re.sub("^\s+", "", line)
        elts = re.split("\s+", line)
        
        #Swig requires floats, but elts is an array of strings
        x=float(elts[0])
        y=float(elts[1])
        flux=float(elts[2])

        source = detect.Source()

        source.setSourceId(i)
        source.setXAstrom(x); source.setXAstromErr(0.1)
        source.setYAstrom(y); source.setYAstromErr(0.1)
        source.setPsfFlux(flux)

        s1.append(source)
        
        i=i + 1
    f.close()
    
    return s1


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class WCSTestCaseNet(unittest.TestCase):
    """A test case for WCS from astrometry.net"""


    def setUp(self):
        pass

    def tearDown(self):
        gas.reset()



    def testCfhtField(self):
        plateScale = .185
        gas.setImageScaleArcsecPerPixel(plateScale)

        #Set starlist  
        starlist=  os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")
        starlist = loadXYFromFile(starlist)
        gas.setStarlist(starlist)

        crval = afwImage.PointD(334.3036567, -17.3301660)
        self.assertTrue(gas.solve(crval))
        
        xy = afwImage.PointD(1,1)
        rd = afwImage.PointD(334.2762566, -17.2097954)
        self.checkMatch(rd, xy)

        xy = afwImage.PointD(1000,1000)
        rd = afwImage.PointD(334.3299915, -17.2602685)
        self.checkMatch(rd, xy)

        xy = afwImage.PointD(25, 912)
        rd = afwImage.PointD(334.2779355, -17.2561995)
        self.checkMatch(rd, xy)

    def checkMatch(self, rdIn, xyIn):
        raDec = gas.xyToRaDec(xyIn.getX(), xyIn.getY())
        #print "Ra's are (%.7f v %.7f)" %(raDec.getX(), rdIn.getX())
        #print "Dec's are (%.7f v %.7f)" %(raDec.getY(), rdIn.getY())

        self.assertAlmostEqual(raDec.getX(), rdIn.getX(), 6, "Ra doesn't match (%.7f v %.7f)" 
            %(raDec.getX(), rdIn.getX()))
        self.assertAlmostEqual(raDec.getY(), rdIn.getY(), 6, "Dec doesn't match (%.7f v %.7f)" 
            %(raDec.getX(), rdIn.getX()))

        #Test the reverse operation
        xy = gas.raDecToXY(rdIn.getX(), rdIn.getY())
        self.assertAlmostEqual(xy.getX(), xyIn.getX(), 2, "x doesn't match (%.7f v %.7f)" 
            %(xy.getX(), xyIn.getX()))
        self.assertAlmostEqual(xy.getY(), xyIn.getY(), 2, "y doesn't match (%.7f v %.7f)" 
            %(xy.getY(), xyIn.getY()))
        
    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(WCSTestCaseNet)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)


#Create a globally accessible instance of a GAS
policyFile=eups.productDir("astrometry_net_data")
policyFile=os.path.join(policyFile, "metadata.paf")
gas = net.GlobalAstrometrySolution(policyFile)
 
if __name__ == "__main__":
    #print "Warning, tests turned off"
    #return 0
    run(True)
