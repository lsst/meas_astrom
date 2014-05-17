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
import lsst.meas.astrom            as measAstrom
import lsst.afw.detection          as afwDet
import lsst.afw.table              as afwTable
import lsst.afw.math               as afwMath
import lsst.afw.image              as afwImg
import lsst.utils.tests            as utilsTests
import lsst.pex.config             as pexConfig
from lsst.pex.logging import Log
from lsst.pex.exceptions import LsstCppException

from lsst.meas.astrom import AstrometryNetDataConfig

class MultiIndexTest(unittest.TestCase):

    def setUp(self):
        self.conf = measAstrom.MeasAstromConfig()

        # Load sample input from disk
        mypath = eups.productDir("meas_astrom")
        path = os.path.join(mypath, "examples")
        self.srcCat = afwTable.SourceCatalog.readFits(os.path.join(path, "v695833-e0-c000.xy.fits"))
        self.srcCat.table.defineApFlux("flux.psf")
        
        # The .xy.fits file has sources in the range ~ [0,2000],[0,4500]
        self.imageSize = (2048, 4612) # approximate
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci.fits"))
        
        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))
        self.an_data_dir = datapath
        
    def tearDown(self):
        del self.srcCat
        del self.conf
        del self.exposure

    def getAstrometrySolution(self, loglvl = Log.INFO, andConfig=None):
        astrom = measAstrom.Astrometry(self.conf, logLevel=loglvl,
                                       andConfig=andConfig)
        res = astrom.determineWcs(self.srcCat, self.exposure,
                                  imageSize=self.imageSize)
        del astrom
        return res

    def _testGetSolution(self, **kwargs):
        res = self.getAstrometrySolution(loglvl=Log.DEBUG, **kwargs)
        self.assertTrue(res is not None)
        self.assertTrue(len(res.getMatches()) > 50)

    # This is the "vanilla" no-multiIndex setup
    def testMultiIndexA(self):
        andConfig = AstrometryNetDataConfig()
        fn = os.path.join(self.an_data_dir, 'andConfig2.py')
        andConfig.load(fn)
        self._testGetSolution(andConfig=andConfig)

    # This is a multiIndex setup with two indices
    def testMultiIndexB(self):
        andConfig = AstrometryNetDataConfig()
        fn = os.path.join(self.an_data_dir, 'andConfig3.py')
        andConfig.load(fn)
        self._testGetSolution(andConfig=andConfig)

    # This is a multiIndex setup with two indices, one that has no star kdtree
    def testMultiIndexC(self):
        andConfig = AstrometryNetDataConfig()
        fn = os.path.join(self.an_data_dir, 'andConfig4.py')
        andConfig.load(fn)
        self._testGetSolution(andConfig=andConfig)

    # This one has a multiIndex and a normal index.
    def testMultiIndexD(self):
        andConfig = AstrometryNetDataConfig()
        fn = os.path.join(self.an_data_dir, 'andConfig5.py')
        andConfig.load(fn)
        self._testGetSolution(andConfig=andConfig)

    # This one has a multiIndex with stars-only and index-only parts
    def testMultiIndexE(self):
        andConfig = AstrometryNetDataConfig()
        fn = os.path.join(self.an_data_dir, 'andConfig6.py')
        andConfig.load(fn)
        self._testGetSolution(andConfig=andConfig)


    # Test that creating an Astrometry object with many index files
    # does not use up a lot of memory or file descriptors.
    # FIXME -- there are no tests on memory usage -- not clear exactly
    # what to enforce.
    def testResources(self):
        from astrometry.util import ttime
        mem0 = ttime.get_memusage()
        fd0 = ttime.count_file_descriptors()
        print
        print 'Mem0:'
        print ttime.memusage()
        print 'FD0:', fd0

        andConfig = AstrometryNetDataConfig()
        fn = os.path.join(self.an_data_dir, 'andConfig6.py')
        andConfig.load(fn)
        andConfig.multiIndexFiles = andConfig.multiIndexFiles * 100
        print len(andConfig.multiIndexFiles), 'multi-index files'

        astrom = measAstrom.Astrometry(self.conf, logLevel=Log.INFO,
                                       andConfig=andConfig)

        mem1 = ttime.get_memusage()
        fd1 = ttime.count_file_descriptors()
        print
        print 'Mem1:'
        print ttime.memusage()
        print 'FD1:', fd1

        # Number of used file descriptors should not grow.  Magic 10
        # is just margin from other things going on in the python process
        self.assertTrue(fd1 < (fd0 + 10))

        res = astrom.determineWcs(self.srcCat, self.exposure,
                                  imageSize=self.imageSize)

        mem2 = ttime.get_memusage()
        fd2 = ttime.count_file_descriptors()
        print
        print 'Mem2:'
        print ttime.memusage()
        print 'FD2:', fd2

        # Number of used file descriptors should not grow.  Magic 10
        # is just margin from other things going on in the python process
        self.assertTrue(fd2 < (fd0 + 10))

        del astrom

        mem3 = ttime.get_memusage()
        fd3 = ttime.count_file_descriptors()
        print
        print 'Mem3:'
        print ttime.memusage()
        print 'FD3:', fd3

        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(MultiIndexTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
