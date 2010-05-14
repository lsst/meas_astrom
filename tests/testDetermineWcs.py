#!/usr/bin/env python
import re
import os
import sys
import glob
import math
import unittest

import eups
import lsst.meas.astrom as measAstrom
import lsst.meas.astrom.net as net
import lsst.afw.detection as det
import lsst.afw.math as afwMath
import lsst.afw.image as afwImg
import lsst.utils.tests as utilsTests
import lsst.pex.policy as pexPolicy

from lsst.pex.exceptions import LsstCppException



class chooseFilterNameTest(unittest.TestCase):
    """The logic for determining which tag-along data to extract from
    the astrometry net catalogue is a little complicated, so 
    this suite tests that it does it right"""

    def setUp(self):
        #Load sample input from disk
        path = os.path.join(eups.productDir("meas_astrom"), "examples")
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci"))
        
        #Set the filter properties.
        afwImg.Filter.reset()
        afwImg.FilterProperty.reset() #Both of these are required to reset
        
        fp = afwImg.FilterProperty("mag")
        afwImg.Filter.define(fp)
        filt = afwImg.Filter("mag")
        self.exposure.setFilter(filt)
        
        #Setup up astrometry_net_data
        print "Setting up meas_astrom cfhtlsDeep"
        eupsObj = eups.Eups()

        ok, version, reason = eupsObj.setup("astrometry_net_data", versionName="cfhtlsDeep")
        if not ok:
            raise ValueError("Couldn't set up cfhtlsDeep version of astrometry_net_data: %s" %(reason))

        #
        ##Load a solver
        metaFile = os.path.join(eups.productDir("astrometry_net_data"), "metadata.paf")
        self.solver = net.GlobalAstrometrySolution(metaFile)
        #
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
        calculateSip: true
        numBrightStars: 75
        defaultFilterName: mag
        wcsToleranceInArcsec: .3
        maxSipOrder: 9
        """
        ))
        
                    
    def tearDown(self):
        del self.exposure
        del self.solver
        del self.defaultPolicy
        pass
                        
    def test1(self):
        """The exposures filtername is one of the filters stored in the catalogue"""
        log=None
        filt = measAstrom.chooseFilterName(self.exposure, self.defaultPolicy, self.solver, log)
        self.assertEqual(filt, "mag")
        

    def test2(self):
        """The exposures filtername is not one of the filters stored in the catalogue, so a default
        is loaded instead
        """

        #Set the filter properties.
        fp = afwImg.FilterProperty("strangelyNamedFilter")
        afwImg.Filter.define(fp)
        filt = afwImg.Filter("strangelyNamedFilter")
        self.exposure.setFilter(filt)

        
        log=None
        filt = measAstrom.chooseFilterName(self.exposure, self.defaultPolicy, self.solver, log)
        self.assertEqual(filt, "mag")
        
        
    def test3(self):
        """Test that we can override the default with another policy"""
        #Set the filter properties.
        fp = afwImg.FilterProperty("strangelyNamedFilter")
        afwImg.Filter.define(fp)
        filt = afwImg.Filter("strangelyNamedFilter")
        self.exposure.setFilter(filt)

        newPolicy = pexPolicy.Policy()
        newPolicy.set("defaultFilterName", "V")
        
        #self.defaultPolicy.mergeDefaults(newPolicy.getDictionary())
        newPolicy.mergeDefaults(self.defaultPolicy)
        
        self.assertEqual( newPolicy.get("defaultFilterName"), "V")

        log=None
        self.assertRaises(ValueError, measAstrom.chooseFilterName, self.exposure, 
                newPolicy, self.solver, log)
                
                
    def test4(self):
        """Test what happens when defaultFilterName isn't defined"""

        #Set the filter properties.
        fp = afwImg.FilterProperty("strangelyNamedFilter")
        afwImg.Filter.define(fp)
        filt = afwImg.Filter("strangelyNamedFilter")
        self.exposure.setFilter(filt)
        
        self.defaultPolicy.remove("defaultFilterName")
        log=None
        filt = measAstrom.chooseFilterName(self.exposure, self.defaultPolicy, self.solver, log)
        self.assertEqual(filt, '')
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(chooseFilterNameTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)



 
if __name__ == "__main__":
    run(True)
