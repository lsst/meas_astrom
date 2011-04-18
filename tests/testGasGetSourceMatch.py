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


"""
Tests of the GlobalAstrometrySolution interface, without actually trying to solve a field.
Solving a field can be slow, but these tests run quickly, so are easier to run to test
small bits of the class. For an end-to-end test, use testGAS.py
"""


import re
import os
import math
import sys
import unittest

import eups
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.meas.astrom.net as net
import lsst.utils.tests as utilsTests
import lsst.afw.image as afwImg
import lsst.afw.detection.detectionLib as detect
import lsst.pex.exceptions as pexExcept
try:
    type(verbose)
except NameError:
    verbose = 0

verbose=True


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


class GlobalAstrometrySolutionTest(unittest.TestCase):
    """These tests work the interface to the class, but don't actually try to solve anything"""


    def setUp(self):
        eupsObj = eups.Eups()

        ok, version, reason = eupsObj.setup("astrometry_net_data", versionName="cfhttemplate")
        if not ok:
            raise ValueError("Couldn't set up cfht version of astrometry_net_data: %s" %(reason))
        
        metaFile = os.path.join(eups.productDir("astrometry_net_data"), "metadata.paf")
        self.gas = net.GlobalAstrometrySolution(metaFile)

        starlistFile = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")
        self.starlist = loadXYFromFile(starlistFile)
        self.gas.setStarlist(self.starlist)

    def tearDown(self):
        del self.gas


    def testOutput(self):
        self.gas.setLogLevel(3)
        flag = self.gas.solve(334.303012, -17.233988)
        #flag = self.gas.solve(wcs)
        self.assertTrue(flag, "Failed to solve")
        sm = self.gas.getMatchedSources()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(GlobalAstrometrySolutionTest)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    verbose = 3
    run(True)
