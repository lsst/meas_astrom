#!/usr/bin/env python
import re
import os
import sys
import glob
import math
import unittest

import eups
import lsst.meas.astrom.net as net
import lsst.afw.detection as det
import lsst.utils.tests as utilsTests

import lsst.meas.astrom.sip as sip
import lsst.meas.astrom.sip.genDistortedImage as distort
import lsst.meas.astrom.sip.cleanBadPoints as cleanBadPoints

import sourceSetIO

#meas_astrom sets up astrometry_net_data if it's available. If it isn't
#available, we don't run the test, because the test would fail, and packages
#shouldn't depend on data packages to build successfully.
dataDir = eups.productDir("astrometry_net_data")
if not dataDir:
    print >> sys.stderr, "Warning: distortedImageTest not run because astrometry_net_data is not set up"
    sys.exit(1)

try:
    if eups.getSetupVersion("astrometry_net_data") != "cfhttemplate":
        setupCmds = eups.setup("astrometry_net_data", "cfhttemplate")
except AttributeError, e:              # old eups
    if eups.Eups().findSetupVersion("astrometry_net_data")[0] == "cfhttemplate":
        setupCmds = ['true']
    else:
        setupCmds = ['false']

if len(setupCmds) == 1 and setupCmds[0] == 'false':
    print >> sys.stderr, "Warning: createWcsWithSip not run because " \
          "astrometry_net_data cfttemplate is unavailable"
    sys.exit(1)

#Create a globally accessible instance of a GAS. This takes a few seconds
#to load, so we don't want to do it everytime we setup a test case
policyFile=eups.productDir("astrometry_net_data")
policyFile=os.path.join(policyFile, "metadata.paf")
GLOBALGAS = net.GlobalAstrometrySolution(policyFile)

class DistortedImageTestCase(unittest.TestCase):
    def setUp(self):
        path=eups.productDir("meas_astrom")
        self.filename=os.path.join(path, "tests", "cat.xy.list")
        self.orderlist = range(2,10) #Fit distortion with every order in this range

    def tearDown(self):
        pass
        
    def testLinearXDistort(self):
        self.singleTestInstance(self.filename, distort.linearXDistort, 
            GLOBALGAS, self.orderlist)

    def testLinearYDistort(self):
        self.singleTestInstance(self.filename, distort.linearYDistort, 
            GLOBALGAS, self.orderlist)

    def testQuadraticDistort(self):
        self.singleTestInstance(self.filename, distort.linearYDistort, 
            GLOBALGAS, self.orderlist)

    #
    # Implementation functions
    #
    
    def singleTestInstance(self, filename, distortFunc, gas, order):
        cat = self.loadCatalogue(self.filename, GLOBALGAS)
        img = distort.distortList(cat, distortFunc)

        #Get a wcs
        gas.setStarlist(img)
        if gas.solve() == True:
            imgWcs = gas.getWcs()
        else:
            print "Failed to solve distorted image. Too much distortion?"
            return
        gas.reset()
    
        #Returns a clean list of matches (i.e with outliers removed)
        matchList = self.matchSrcAndCatalogue(cat, img, imgWcs)    
        x, y, dx, dy = self.get1dVectorsFromSourceMatchVector(matchList)

        s = list(.1 for i in range(len(x)))
        for i in order:
            lsfX = sip.LeastSqFitter2dPoly(x, y, dx, s, i);
            lsfY = sip.LeastSqFitter2dPoly(x, y, dy, s, i);
            
            for j in range(len(x)):
                self.assertAlmostEqual(dx[j], lsfX.valueAt(x[j], y[j]), 3)
                self.assertAlmostEqual(dy[j], lsfY.valueAt(x[j], y[j]), 3)
        

    def loadCatalogue(self, filename, gas):
        """Load a list of xy points from a file, solve for position, and
        return a SourceSet of points"""


        cat = sourceSetIO.read(filename)

        gas.setStarlist(cat)
        if gas.solve() is True:
            catWcs = gas.getWcs()
            gas.reset()        
        else:
            gas.reset()
            raise "Failed to solve catalogue"

        #Set catalogue ra and decs
        for src in cat:
            raDec = catWcs.pixelToSky(src.getXAstrom(), src.getYAstrom())
            src.setRa(math.degrees(raDec[0]))
            src.setDec(math.degrees(raDec[1]))

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
        


    def get1dVectorsFromSourceMatchVector(self, mList):
        """Does what it says on the tin"""

        xFunc = lambda x: (x[0]).getXAstrom()
        x = map(xFunc, mList)

        yFunc = lambda y: (y[0]).getYAstrom()
        y = map(yFunc, mList)

        dxFunc = lambda x: (x[0]).getXAstrom() - (x[1]).getXAstrom()
        dx = map(dxFunc, mList)

        dyFunc = lambda y: (y[0]).getYAstrom() - (y[1]).getYAstrom()
        dy = map(dyFunc, mList)

        return x, y, dx, dy

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
    return utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
