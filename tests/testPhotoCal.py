#!/usr/bin/env python
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
                self.srcSet)
        
           
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
