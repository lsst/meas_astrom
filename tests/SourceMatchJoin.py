#!/usr/bin/env python

import os
import unittest
import lsst.utils.tests as utilsTests
import lsst.afw.table as afwTable
import lsst.meas.astrom as measAstrom
from lsst.pex.logging import Log
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImg

import testFindAstrometryNetDataDir as helper

try:
    type(verbose)
except NameError:
    display = False
    verbose = 0


class matchlistTestCase(unittest.TestCase):
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
        #logLevel = Log.DEBUG
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

        matches2 = self.astrom.joinMatchListWithCatalog(normalized, self.srcSet)

        self.assertEqual(len(matches2), len(matches))
        for i in xrange(len(matches)):
            self.assertEqual(matches2[i].second.table, matches[i].second.table)
            self.assertEqual(matches2[i].second.getId(), matches[i].second.getId())
            self.assertEqual(matches2[i].second, matches[i].second) # no deep copying, so we can compare ptrs
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

        matches2 = self.astrom.joinMatchListWithCatalog(normalized, self.srcSet)
        self.assertGreater(len(matches2), 0)
        ref = matches2[0][0]

        names = ref.getSchema().getNames()
        for b in ("u", "g", "r", "i", "z"):
            self.assertTrue("%s_flux" % (b,) in names)
            self.assertTrue("%s_fluxSigma" % (b,) in names)
            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= silly boilerplate -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(matchlistTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
