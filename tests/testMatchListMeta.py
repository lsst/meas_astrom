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
import tempfile
import numpy as np
from math import ceil

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
from lsst.pex.logging import Log

from astrometry.util.pyfits_utils import fits_table
import sourceSetIO as ssi



class MatchListMetaTest(unittest.TestCase):

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
        defaultFilterName: rmag
        sipOrder: 4
        """
        ))

        mypath = eups.productDir('meas_astrom')
        tests = os.path.join(mypath, 'tests')
        self.testdir = tests

        # Create a fake blank image of the appropriate size
        (W,H) = (2048,1489)
        self.exposure = afwImg.ExposureF(W, H)

        # Grab source list
        T = fits_table(os.path.join(tests, 'fpobjc-0745-3-0564-cut.fits'))
        sourceSet = det.SourceSet()
        for i,(x,y,f) in enumerate(zip(T.x, T.y, T.flux)):
            s = det.Source()
            s.setId(i)
            s.setFlagForDetection(8320)
            s.setRa(0.)
            s.setDec(0.)
            s.setXAstrom(x)
            s.setYAstrom(y)
            s.setPsfFlux(f)
            sourceSet.append(s)
        self.srcSet = sourceSet

        print 'Exposure image size: %i x %i' % (self.exposure.getWidth(), self.exposure.getHeight())
        print '%i sources' % len(sourceSet)

        # Set up local astrometry_net_data
        datapath = os.path.join(tests, 'astrometry_net_data', 'tagalong')
        # Work around lame scons bug (doesn't pass HOME)
        os.environ['HOME'] = 'iswheretheheartis'
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Couldn't set up local 'tagalong' version of astrometry_net_data (from path: %s): %s" % (datapath, reason))

    def tearDown(self):
        del self.defaultPolicy
        del self.exposure
        del self.srcSet

                        
    def test1(self):
        log = Log.getDefaultLog()
        gas = measAstrom.createSolver(self.defaultPolicy, log)

        astrom = measAstrom.determineWcs(self.defaultPolicy, self.exposure,
                                         self.srcSet, solver=gas)
        matchMeta = astrom.getMatchMetadata()
        print 'matchMeta:', matchMeta.toString()

        self.assertEqual(matchMeta.getAsString('REFCAT'), 'sdss-tsobj-0745-3-40-0564')
        self.assertEqual(matchMeta.getAsString('REFCAMD5'), 'none')
        self.assertEqual(matchMeta.getAsInt('ANINDID'), 9999)
        self.assertEqual(matchMeta.getAsInt('ANINDHP'), -1)
        indexname = matchMeta.getAsString('ANINDNM')
        self.assertTrue('tests/astrometry_net_data/tagalong/tsobj-0745-3-40-0564.index' in indexname)
        self.assertAlmostEqual(matchMeta.getAsDouble('RA'), 243.1766, 2)
        self.assertAlmostEqual(matchMeta.getAsDouble('DEC'), -0.1095, 2)
        self.assertAlmostEqual(matchMeta.getAsDouble('RADIUS'), 0.1392, 2)





#-=-=-=-=-=-=-=-=-=-=-=-=-=-=  silly boilerplate -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(MatchListMetaTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
