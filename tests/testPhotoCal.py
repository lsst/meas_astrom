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
import lsst.utils.tests as utilsTests
import lsst.pex.policy as pexPolicy
import lsst.meas.photocal as photocal

from lsst.pex.exceptions import LsstCppException

import sourceSetIO as ssi


class PhotoCalTest(unittest.TestCase):

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
        path = os.path.join(eups.productDir("meas_astrom"), "examples")
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci"))
        self.srcSet = ssi.read(os.path.join(path, "v695833-e0-c000.xy.txt"))

        # The .xy.txt file has sources in the range ~ [0,2000],[0,4500], but
        # the exposure is onlyl one amp -- 1024x1153.  Work around.
        print 'Exposure image size: %i x %i' % (self.exposure.getWidth(), self.exposure.getHeight())
        self.forceImageSize = (2048, 4612) # approximately; 2x4 x (1024 x 1153)
        print 'Forcing image size to %i x %i to match source list.' % (self.forceImageSize[0], self.forceImageSize[1])

        #Setup up astrometry_net_data
        #print "Setting up meas_astrom cfhtlsDeep"
        eupsObj = eups.Eups()

        ok, version, reason = eupsObj.setup("astrometry_net_data", versionName="cfhtlsDeep")
        if not ok:
            raise ValueError("Couldn't set up cfhtlsDeep version of astrometry_net_data: %s" %(reason))

    def tearDown(self):
        del self.defaultPolicy
        del self.exposure
        del self.srcSet

                        
    def test1(self):
        matches, wcs = measAstrom.determineWcs(self.defaultPolicy, self.exposure, \
                self.srcSet, forceImageSize=self.forceImageSize)
        
           
        pCal = photocal.calcPhotoCal(matches)
        print pCal


        diff=[]
        for m in matches:
            catFlux = m[0].getPsfFlux()     #Catalogue flux
            catMag = -2.5*np.log10(catFlux) #Cat mag
            instFlux = m[1].getPsfFlux()    #Instrumental Flux
            mag = pCal.getMag(instFlux)     #Instrumental mag
            
            diff.append(mag-catMag)


        #A very loose test, but the input data has a lot of scatter

        diff = np.array(diff)
        self.assertAlmostEqual(np.mean(diff), 0, 0)
        
        
    def test2(self):
        """Check that negative fluxes dealt with properly"""
        
        matches, wcs = measAstrom.determineWcs(self.defaultPolicy, self.exposure, \
                self.srcSet, forceImageSize=self.forceImageSize)

        matches[0].first.setPsfFlux(0)
        matches[1].first.setPsfFlux(-1)
        matches[2].second.setPsfFlux(0)
        
        pCal = photocal.calcPhotoCal(matches)
        print pCal


        diff=[]
        for m in matches:
            catFlux = m[0].getPsfFlux()     #Catalogue flux
            instFlux = m[1].getPsfFlux()    #Instrumental Flux
            
            if catFlux > 0 and instFlux > 0:
                catMag = -2.5*np.log10(catFlux) #Cat mag
                mag = pCal.getMag(instFlux)     #Instrumental mag
                diff.append(mag-catMag)


        #A very loose test, but the input data has a lot of scatter

        diff = np.array(diff)
        self.assertAlmostEqual(np.mean(diff), 0, 0)
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PhotoCalTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)



 
if __name__ == "__main__":
    run(True)
