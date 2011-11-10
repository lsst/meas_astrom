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

import numpy as np

import eups
import lsst.meas.astrom as measAstrom
import lsst.meas.astrom.net as net
import lsst.afw.detection as det
import lsst.afw.math as afwMath
import lsst.afw.image as afwImg
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.utils.tests as utilsTests
import lsst.pex.policy as pexPolicy
import lsst.meas.photocal as photocal

from lsst.pex.exceptions import LsstCppException

import sourceSetIO as ssi


class TagAlongTest(unittest.TestCase):

    def setUp(self):
        self.defaultPolicy = pexPolicy.Policy.createPolicy(pexPolicy.PolicyString(
        """#<?cfg paf policy?>     
        inputExposureKey: visitExposure
        inputSourceSetKey: sourceSet
        allowDistortion: true
        matchThreshold: 22
        blindSolve: false
        outputWcsKey: measuredWcs
        outputMatchListKey: matchList
        distanceForCatalogueMatchinArcsec: 1.0
        cleaningParameter: 3
        calculateSip: false
        numBrightStars: 75
        defaultFilterName: mag
        sipOrder: 4
        """
        ))
        
        #Load sample input from disk
        mypath = eups.productDir("meas_astrom")
        path = os.path.join(mypath, "examples")
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci"))
        self.srcSet = ssi.read(os.path.join(path, "v695833-e0-c000.xy.txt"))

        # The .xy.txt file has sources in the range ~ [0,2000],[0,4500], but
        # the exposure is only one amp -- 1024x1153.  Work around.
        print 'Exposure image size: %i x %i' % (self.exposure.getWidth(), self.exposure.getHeight())
        self.forceImageSize = (2048, 4612) # approximately; 2x4 x (1024 x 1153)
        print 'Forcing image size to %i x %i to match source list.' % (self.forceImageSize[0], self.forceImageSize[1])

        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'testTagAlong')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Couldn't set up local photocal version of astrometry_net_data (from path: %s): %s" % (datapath, reason))

    def tearDown(self):
        del self.defaultPolicy
        del self.exposure
        del self.srcSet

                        
    def test1(self):
        path = os.path.join(os.environ['ASTROMETRY_NET_DATA_DIR'], "metadata.paf")
        gas = measAstrom.net.GlobalAstrometrySolution(path)
        matchThreshold = self.defaultPolicy.get('matchThreshold')
        gas.setMatchThreshold(matchThreshold)

        #matches, wcs = measAstrom.determineWcs(self.defaultPolicy, self.exposure,
        #                                       self.srcSet, forceImageSize=self.forceImageSize,
        #                                       solver=gas)
        #c = wcs.pixelToSky(self.forceImageSize[0]/2, self.forceImageSize[1]/2)
        #ra,dec = c.getLongitude(afwCoord.DEGREES), c.getLatitude(afwCoord.DEGREES)
        #radius = wcs.pixArea(afwGeom.Point2D(self.forceImageSize[0]/2, self.forceImageSize[1]/2))
        #(xyz,radec,inds,tag) = gas.getIndexStars(ra, dec, radius)

        (ra,dec,radius) = (-145, 53, 0.02)
        print 'Searching RA,Dec %g,%g, radius %g deg' % (ra,dec,radius)
        (xyz,radec,inds,tag) = gas.getIndexStars(ra,dec,radius)
        print 'Found %i index stars' % len(xyz)
        print 'Found tag-along columns:', tag.keys()

        self.assertEqual(9, len(xyz))
        self.assertEqual(len(xyz), len(radec))
        self.assertEqual(len(xyz), len(inds))
        self.assertTrue('mag' in tag)
        self.assertEqual(len(tag['mag']), len(xyz))
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TagAlongTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
