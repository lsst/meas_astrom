from builtins import range
#!/usr/bin/env python
#
# LSST Data Management System
#
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

import lsst.afw.geom as afwGeom
from lsst.afw.image import ExposureF
from lsst.afw.table import packMatches, SourceCatalog
import lsst.utils.tests
from lsst.daf.persistence import Butler
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
from lsst.meas.astrom import AstrometryTask
from lsst.log import Log


class joinMatchListWithCatalogTestCase(unittest.TestCase):

    def setUp(self):
        # Load sample input from disk
        testDir = os.path.dirname(__file__)

        self.srcSet = SourceCatalog.readFits(os.path.join(testDir, "v695833-e0-c000.xy.fits"))

        self.bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(2048, 4612))  # approximate
        # create an exposure with the right metadata; the closest thing we have is
        # apparently v695833-e0-c000-a00.sci.fits, which is much too small
        smallExposure = ExposureF(os.path.join(testDir, "v695833-e0-c000-a00.sci.fits"))
        self.exposure = ExposureF(self.bbox)
        self.exposure.setWcs(smallExposure.getWcs())
        self.exposure.setFilter(smallExposure.getFilter())
        # copy the pixels we can, in case the user wants a debug display
        mi = self.exposure.getMaskedImage()
        mi.assign(smallExposure.getMaskedImage(), smallExposure.getBBox())

        logLevel = Log.INFO
        #logLevel = Log.DEBUG
        refCatDir = os.path.join(testDir, "data", "sdssrefcat")
        butler = Butler(refCatDir)
        refObjLoader = LoadIndexedReferenceObjectsTask(butler=butler)
        astrometryConfig = AstrometryTask.ConfigClass()
        self.astrom = AstrometryTask(config=astrometryConfig, refObjLoader=refObjLoader)
        self.astrom.log.setLevel(logLevel)
        # Since our sourceSelector is a registry object we have to wait for it to be created
        # before setting default values.
        self.astrom.matcher.sourceSelector.config.minSnr = 0

    def tearDown(self):
        del self.srcSet
        del self.bbox
        del self.exposure
        del self.astrom

    def getAstrometrySolution(self):
        return self.astrom.solve(exposure=self.exposure, sourceCat=self.srcSet)

    def testJoin(self):
        res = self.getAstrometrySolution()

        matches = res.matches
        matchmeta = res.matchMeta

        normalized = packMatches(matches)
        normalized.table.setMetadata(matchmeta)

        matches2 = self.astrom.refObjLoader.joinMatchListWithCatalog(normalized, self.srcSet)

        self.assertEqual(len(matches2), len(matches))
        for i in range(len(matches)):
            self.assertEqual(matches2[i].second.table, matches[i].second.table)
            self.assertEqual(matches2[i].second.getId(), matches[i].second.getId())
            self.assertEqual(matches2[i].second, matches[i].second)  # no deep copying, so we can compare ptrs
            self.assertEqual(matches2[i].first.getId(), matches[i].first.getId())
            self.assertEqual(matches2[i].first.getRa().asDegrees(), matches[i].first.getRa().asDegrees())
            self.assertEqual(matches2[i].first.getDec().asDegrees(), matches[i].first.getDec().asDegrees())
            self.assertEqual(matches2[i].first.get("i_flux"), matches[i].first.get("i_flux"))

    def testJoinAllFluxes(self):
        """Test that we can read all the fluxes back from an a.n.d catalogue"""
        res = self.getAstrometrySolution()

        matches = res.matches
        matchmeta = res.matchMeta

        normalized = packMatches(matches)
        normalized.table.setMetadata(matchmeta)

        matches2 = self.astrom.refObjLoader.joinMatchListWithCatalog(normalized, self.srcSet)
        self.assertGreater(len(matches2), 0)
        ref = matches2[0][0]

        names = ref.getSchema().getNames()
        for b in ("u", "g", "r", "i", "z"):
            self.assertIn("%s_flux" % (b,), names)
            self.assertIn("%s_fluxSigma" % (b,), names)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
