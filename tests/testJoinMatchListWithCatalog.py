#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import os
import unittest

import lsst.utils.tests as utilsTests
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImg
import lsst.afw.table as afwTable
import lsst.meas.astrom as measAstrom
from lsst.pex.logging import Log

import testFindAstrometryNetDataDir as helper

class joinMatchListTestCase(unittest.TestCase):
    def setUp(self):
        # Load sample input from disk
        testDir=os.path.dirname(__file__)
        # Set up local astrometry_net_data
        helper.setupAstrometryNetDataDir('photocal', verbose=True)

        self.srcSet = afwTable.SourceCatalog.readFits(os.path.join(testDir, "v695833-e0-c000.xy.fits"))
        self.bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(2048, 4612)) # approximate
        self.exposure = afwImg.ExposureF(os.path.join(testDir, "v695833-e0-c000-a00.sci.fits"))

        config = measAstrom.ANetBasicAstrometryConfig()
        logLevel = Log.INFO
        self.astrom = measAstrom.ANetBasicAstrometryTask(config)
        self.astrom.log.setThreshold(logLevel)

    def tearDown(self):
        del self.srcSet
        del self.bbox
        del self.exposure
        del self.astrom

    def getAstrometrySolution(self):
        return self.astrom.determineWcs(self.srcSet, self.exposure, bbox=self.bbox)

    def testJoin(self):
        res = self.getAstrometrySolution()

        matches = res.matches
        matchmeta = res.matchMeta

        normalized = afwTable.packMatches(matches)
        normalized.table.setMetadata(matchmeta)

        # Test function with AstrometryTask passed as caller initilized and None
        for astrom in (self.astrom, None):
            matches2 = measAstrom.joinMatchListWithCatalog(normalized, self.srcSet, astrom)

            self.assertEqual(len(matches2), len(matches))
            for i in xrange(len(matches)):
                self.assertEqual(matches2[i].second.table, matches[i].second.table)
                self.assertEqual(matches2[i].second.getId(), matches[i].second.getId())
                self.assertEqual(matches2[i].second, matches[i].second) # no deep copy, so we can compare ptrs
                self.assertEqual(matches2[i].first.getId(), matches[i].first.getId())
                self.assertEqual(matches2[i].first.getRa().asDegrees(), matches[i].first.getRa().asDegrees())
                self.assertEqual(matches2[i].first.getDec().asDegrees(), matches[i].first.getDec().asDegrees())
                self.assertEqual(matches2[i].first.get("i_flux"), matches[i].first.get("i_flux"))

    def testJoinAllFluxes(self):
        """Test that we can read all the fluxes back from an a.n.d catalogue"""
        res = self.getAstrometrySolution()

        matches = res.matches
        matchmeta = res.matchMeta

        normalized = afwTable.packMatches(matches)
        normalized.table.setMetadata(matchmeta)

        # Test function with AstrometryTask passed as caller initilized and None
        for astrom in (self.astrom, None):
            matches2 = measAstrom.joinMatchListWithCatalog(normalized, self.srcSet, astrom)
            self.assertGreater(len(matches2), 0)
            ref = matches2[0][0]

            names = ref.getSchema().getNames()
            for b in ("u", "g", "r", "i", "z"):
                self.assertTrue("%s_flux" % (b,) in names)
                self.assertTrue("%s_fluxSigma" % (b,) in names)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(joinMatchListTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the utilsTests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
