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

import lsst.utils.tests as utilsTests
from lsst.meas.photocal.PhotometricMagnitude import PhotometricMagnitude

class PhotometricMagnitudeTest(unittest.TestCase):

    def setUp(self):
        self.obj = PhotometricMagnitude(zeroFlux=15000, zeroMag=12.0)
    
    def testStrRep(self):
        testStr = "PhotometricMagnitude: magnitude 12 at 15000"
        self.assertEqual( str(self.obj), testStr)
        
    def testToMagScalar(self):
        mag = self.obj.getMag(2e4)
        self.assertAlmostEqual(mag, 11.687, 2)

    def testToMagList(self):
        mag = self.obj.getMag( [2e3, 2e4, 2e5])
        expected = [14.188, 11.687, 9.188]
        
        for i in range(len(mag)):
            self.assertAlmostEqual(mag[i], expected[i], 2)


    def testToMagNumpyArray(self):
        flux = np.array([2e3, 2e4, 2e5])
        mag = self.obj.getMag(flux)
        expected = [14.188, 11.687, 9.188]
        
        self.assertEqual(type(mag), type(flux))
        for i in range(len(mag)):
            self.assertAlmostEqual(mag[i], expected[i], 2)



    def testToFluxList(self):
        mag = [14, 16, 18]
        expected = [2377.3, 376.8, 59.7]
        
        flux = self.obj.getFlux(mag)
        for i in range(len(mag)):
            self.assertAlmostEqual(flux[i], expected[i], 1)
   
    def testFunctor(self):
        """Test the object-as-a-function interface"""
        mag = self.obj(15000)
        self.assertAlmostEqual(mag, 12.0, 6)     
        
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PhotometricMagnitudeTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)



 
if __name__ == "__main__":
    run(True)
