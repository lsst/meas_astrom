#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import re
import os
import sys
import glob
import math
import unittest

import eups
import lsst.meas.astrom.net as net
import lsst.afw.detection as det
import lsst.afw.math as afwMath
import lsst.utils.tests as utilsTests
import lsst.afw.geom as afwGeom

import lsst.meas.astrom.sip as sip
import lsst.meas.astrom.sip.genDistortedImage as distort
import lsst.meas.astrom.sip.cleanBadPoints as cleanBadPoints

import sourceSetIO

############################
# Set up local astrometry_net_data
meas_astrom_dir = eups.productDir("meas_astrom")
datapath = os.path.join(meas_astrom_dir, 'tests', 'astrometry_net_data', 'cfhttemplate')

# Work around lame scons bug (doesn't pass HOME)
os.environ['HOME'] = 'iswheretheheartis' 
eupsObj = eups.Eups(root=datapath)

ok, version, reason = eupsObj.setup('astrometry_net_data')

if not ok:
    raise ValueError("Can't find astrometry_net_data version cfhttemplate (from path: %s): %s" %
                     (datapath, reason))

#Create a globally accessible instance of a GAS. This takes a few seconds
#to load, so we don't want to do it everytime we setup a test case
policyFile=eups.productDir("astrometry_net_data")
policyFile=os.path.join(policyFile, "metadata.paf")
print "GLOBALGAS"
GLOBALGAS = net.GlobalAstrometrySolution(policyFile)
print "...done"


class CreateWcsWithSipCase(unittest.TestCase):
    def setUp(self):
        path=eups.productDir("meas_astrom")
        self.filename=os.path.join(path, "tests", "cat.xy.list")
        self.tolArcsec = .4 
        self.tolPixel = .1

    def tearDown(self):
        pass
        
    def testLinearXDistort(self):
	print "linearXDistort"
        self.singleTestInstance(self.filename, distort.linearXDistort, 
            GLOBALGAS)

    def testLinearYDistort(self):
	print "linearYDistort"
        self.singleTestInstance(self.filename, distort.linearYDistort, 
            GLOBALGAS)

    def testQuadraticDistort(self):
	print "linearQuadraticDistort"
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

        print 'matchList:'
        for m in matchList:
            print m
        
        sipObject = sip.CreateWcsWithSip(matchList, imgWcs, 3)
        imgWcs = sipObject.getNewWcs()

        print 'number of matches:', len(matchList), sipObject.getNPoints()
        scatter = sipObject.getScatterOnSky().asArcseconds()
        print "Scatter in arcsec is %g" % (scatter)

        if scatter >= self.tolArcsec:
            print 'matches:'
            for m in matchList:
                print '  ', m
                
        self.assertTrue(scatter < self.tolArcsec, "Scatter exceeds tolerance in arcsec")

        if False:
            scatter = sipObject.getScatterInPixels()
            self.assertTrue(scatter < self.tolPixel, "Scatter exceeds tolerance in pixels: %g" %(scatter))
        

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
            src.setRaDecFromXy(catWcs)
        return cat

    def matchSrcAndCatalogue(self, cat, img, imgWcs, distInArcsec=1.0, cleanParam=3):
        """Given an input catalogue, match a list of objects in an image, given
        their x,y position and a wcs solution.

        Return: A list of x, y, dx and dy. Each element of the list is itself a list
        """

        matcher = sip.MatchSrcToCatalogue(cat, img, imgWcs, distInArcsec * afwGeom.arcseconds)
        matchList = matcher.getMatches()

        mList = cleanBadPoints.clean(matchList, imgWcs, order=cleanParam)
        return mList
        



        

        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CreateWcsWithSipCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)



 
if __name__ == "__main__":
    run(True)
