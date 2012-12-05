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
import lsst.meas.astrom as measAstrom
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
        mypath = eups.productDir("meas_astrom")
        path = os.path.join(mypath, "examples")
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci.fits"))
        
        # Set the filter properties.
        afwImg.Filter.reset()
        afwImg.FilterProperty.reset() #Both of these are required to reset
        
        #fp = afwImg.FilterProperty("mag")
        #afwImg.Filter.define(fp)
        #filt = afwImg.Filter("mag")
        #self.exposure.setFilter(filt)

        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))
        # This "photocal" index has ugriz columns.

        self.conf = measAstrom.MeasAstromConfig()
        self.astrom = measAstrom.Astrometry(self.conf)
                    
    def tearDown(self):
        del self.exposure
        del self.astrom
        import lsst.meas.astrom.astrometry_net as an
        an.stop_an_logging()
                        
    def test1(self):
        """The exposures filtername is one of the filters stored in the catalogue"""

        print 'filter map:', self.astrom.andConfig.magColumnMap

        filterName = self.exposure.getFilter().getName()
        print 'filter name:', filterName
        print type(filterName)
        filt = self.astrom.getCatalogFilterName(filterName)
        self.assertEqual(filt, 'i')

    def test2(self):
        """The exposures filtername is not one of the filters stored in the catalogue, so a default
        is loaded instead
        """

        #Set the filter properties.
        fp = afwImg.FilterProperty("strangelyNamedFilter")
        afwImg.Filter.define(fp)
        filt = afwImg.Filter("strangelyNamedFilter")
        self.exposure.setFilter(filt)

        filterName = self.exposure.getFilter().getName()
        filt = self.astrom.getCatalogFilterName(filterName)
        self.assertEqual(filt, 'r')
        
    def test3(self):
        """Test that we can override the default with another policy"""
        #Set the filter properties.
        fp = afwImg.FilterProperty("strangelyNamedFilter")
        afwImg.Filter.define(fp)
        filt = afwImg.Filter("strangelyNamedFilter")
        self.exposure.setFilter(filt)

        # Update config.
        self.astrom.andConfig.magColumnMap['strangelyNamedFilter'] = 'z'

        filterName = self.exposure.getFilter().getName()
        filt = self.astrom.getCatalogFilterName(filterName)
        self.assertEqual(filt, 'z')
                
        
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
