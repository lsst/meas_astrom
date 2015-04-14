#!/usr/bin/env python

import os
import unittest
import eups
import lsst.utils.tests as utilsTests
import lsst.afw.table as afwTable
import lsst.meas.astrom as measAstrom
from lsst.pex.logging import Log
import lsst.afw.image as afwImg

try:
    type(verbose)
except NameError:
    display = False
    verbose = 0


class matchlistTestCase(unittest.TestCase):
    def setUp(self):
        # Load sample input from disk
        mypath = eups.productDir("meas_astrom")
        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        print 'Setting up astrometry_net_data:', datapath
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Couldn't set up local photocal version of astrometry_net_data (from path: %s): %s" % (datapath, reason))

        path = os.path.join(mypath, "examples")
        self.srcSet = afwTable.SourceCatalog.readFits(os.path.join(path, "v695833-e0-c000.xy.fits"))
        self.imageSize = (2048, 4612) # approximate
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci.fits"))

        conf = measAstrom.ANetBasicAstrometryConfig()
        loglvl = Log.INFO
        #loglvl = Log.DEBUG
        self.astrom = measAstrom.ANetBasicAstrometryTask(conf, logLevel=loglvl)

    def tearDown(self):
        del self.srcSet
        del self.imageSize
        del self.exposure
        del self.astrom
        
    def getAstrometrySolution(self):
        return self.astrom.determineWcs(self.srcSet, self.exposure, imageSize=self.imageSize)

    def testJoin(self):
        res = self.getAstrometrySolution()

        matches = res.matches
        matchmeta = res.matchMetadata

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
        matchmeta = res.matchMetadata

        normalized = afwTable.packMatches(matches)
        normalized.table.setMetadata(matchmeta)

        matches2 = self.astrom.joinMatchListWithCatalog(normalized, self.srcSet)
        self.assertTrue(len(matches2) > 0)
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
