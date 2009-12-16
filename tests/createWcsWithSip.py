#!/usr/bin/env python
import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.meas.astrom.net as net
import lsst.afw.detection as det
import lsst.afw.math as afwMath
import lsst.utils.tests as utilsTests

import lsst.meas.astrom.sip as sip
import lsst.meas.astrom.sip.genDistortedImage as distort
import lsst.meas.astrom.sip.cleanBadPoints as cleanBadPoints

import sourceSetIO




dataDir = eups.productDir("astrometry_net_data")
if not dataDir:
    raise RuntimeError("Must set up astrometry_net_data to run these tests")


#Create a globally accessible instance of a GAS. This takes a few seconds
#to load, so we don't want to do it everytime we setup a test case
policyFile=eups.productDir("astrometry_net_data")
policyFile=os.path.join(policyFile, "metadata.paf")
GLOBALGAS = net.GlobalAstrometrySolution(policyFile)


class DistortedImageTestCase(unittest.TestCase):
    def setUp(self):
        path=eups.productDir("meas_astrom")
        self.filename=os.path.join(path, "tests", "cat.xy.list")
        self.tolArcsec = .3 
        self.tolPixel = .01

    def tearDown(self):
        pass
        
    def testLinearXDistort(self):
        self.singleTestInstance(self.filename, distort.linearXDistort, 
            GLOBALGAS)

    def testLinearYDistort(self):
        self.singleTestInstance(self.filename, distort.linearYDistort, 
            GLOBALGAS)

    def testQuadraticDistort(self):
        self.singleTestInstance(self.filename, distort.linearYDistort, 
            GLOBALGAS)


    #
    # Implementation functions
    #
    
    def singleTestInstance(self, filename, distortFunc, gas):
        cat = self.loadCatalogue(self.filename, GLOBALGAS)
        img = distort.distortList(cat, distortFunc)

        #Get a wcs
        gas.setStarlist(img)
        flag = gas.solve()
        self.assertTrue(flag, "Failed to solve distorted image. Too much distortion?")
        imgWcs = gas.getWcs()
        gas.reset()
    
        #Create a wcs with sip
        matchList = self.matchSrcAndCatalogue(cat, img, imgWcs)
        sipObject = sip.CreateWcsWithSip(matchList, imgWcs, 3)
        imgWcs = sipObject.getNewWcs()

        scatter = sipObject.getScatterInArcsec()
        self.assertTrue(scatter < self.tolArcsec, "Scatter exceeds tolerance in arcsec")

        scatter = sipObject.getScatterInPixels()
        self.assertTrue(scatter < self.tolPixel, "Scatter exceeds tolerance in pixels")
        

    def loadCatalogue(self, filename, gas):
        """Load a list of xy points from a file, solve for position, and
        return a SourceSet of points"""


        cat = sourceSetIO.read(filename)

        gas.setStarlist(cat)
        flag = gas.solve()
        if flag == True:
            catWcs = gas.getWcs()
            gas.reset()        
        else:
            gas.reset()
            raise "Failed to solve catalogue"

        #Set catalogue ra and decs
        for src in cat:
            raDec = catWcs.xyToRaDec(src.getXAstrom(), src.getYAstrom())
            src.setRa(raDec[0])
            src.setDec(raDec[1])

        return cat

    def matchSrcAndCatalogue(self, cat, img, imgWcs, distInArcsec=1.0, cleanParam=3):
        """Given an input catalogue, match a list of objects in an image, given
        their x,y position and a wcs solution.

        Return: A list of x, y, dx and dy. Each element of the list is itself a list
        """

        matcher = sip.MatchSrcToCatalogue(cat, img, imgWcs, distInArcsec)    
        matchList = matcher.getMatches()

        mList = cleanBadPoints.clean(matchList, imgWcs, order=cleanParam)
        return mList
        



        

        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(DistortedImageTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)



 
if __name__ == "__main__":
    run(True)
