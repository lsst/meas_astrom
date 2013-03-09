#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2013 Dustin Lang.
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

import resource

import eups
import lsst.meas.astrom as measAstrom
import lsst.utils.tests as utilsTests
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImg
import lsst.pex.logging as pexLog


class OpenFilesTest(unittest.TestCase):

    def setUp(self):
        print 'NOFILE rlimit:', resource.getrlimit(resource.RLIMIT_NOFILE)
        resource.setrlimit(resource.RLIMIT_NOFILE, (25, -1))
        print 'NOFILE rlimit:', resource.getrlimit(resource.RLIMIT_NOFILE)

        #R = resource.getrusage(resource.RUSAGE_SELF)
        #print 'NOFILE usage:', 

        conf = measAstrom.MeasAstromConfig()
        mypath = os.path.dirname(os.path.dirname(__file__))
        path = os.path.join(mypath, "examples")
        self.srcCat = afwTable.SourceCatalog.readFits(os.path.join(path, "v695833-e0-c000.xy.fits"))
        self.srcCat.table.defineApFlux("flux.psf")
        # The .xy.fits file has sources in the range ~ [0,2000],[0,4500]
        self.imageSize = (2048, 4612) # approximate
        #self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci"))

        andpath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        os.environ['ASTROMETRY_NET_DATA_DIR'] = andpath
        andcfn = os.path.join(andpath, 'andConfigOpenFiles.py')

        andconfig = measAstrom.AstrometryNetDataConfig()
        andconfig.load(andcfn)

        self.astrom = measAstrom.Astrometry(config=conf, andConfig=andconfig,
                                            logLevel=pexLog.Log.DEBUG)

    def tearDown(self):
        del self.astrom
        del self.imageSize
        del self.srcCat
        
        
    def test1(self):
        res = self.astrom.determineWcs2(self.srcCat, imageSize=self.imageSize,
                                        filterName='i')
        print 'Got result', res
    

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(OpenFilesTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
    


