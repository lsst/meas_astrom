# This file is part of meas_astrom.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import unittest

import numpy as np

import lsst.pipe.base
import lsst.utils.tests
import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from lsst.meas.algorithms import LoadReferenceObjectsTask
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.astrom import FitAffineWcsTask, TransformedSkyWcsMaker


class BaseTestCase:

    """A test case for FitAffineWcsTask.

    Use involves setting one class attribute:
    * MatchClass: match class, e.g. ReferenceMatch or SourceMatch
    """
    MatchClass = None

    def setUp(self):
        crval = lsst.geom.SpherePoint(44, 45, lsst.geom.degrees)
        crpix = lsst.geom.Point2D(15000, 4000)

        scale = 1 * lsst.geom.arcseconds
        cdMatrix = afwGeom.makeCdMatrix(scale=scale, flipX=True)
        self.tanWcs = afwGeom.makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cdMatrix)
        self.loadData()

    def loadData(self, rangePix=3000, numPoints=25):
        """Load catalogs and make the match list

        This is a separate function so data can be reloaded if fitting more than once
        (each time a WCS is fit it may update the source catalog, reference catalog and match list)
        """
        refSchema = LoadReferenceObjectsTask.makeMinimalSchema(
            filterNameList=["r"], addIsPhotometric=True, addCentroid=True)
        self.refCat = afwTable.SimpleCatalog(refSchema)
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        SingleFrameMeasurementTask(schema=srcSchema)
        self.srcCoordKey = afwTable.CoordKey(srcSchema["coord"])
        self.srcCentroidKey = afwTable.Point2DKey(srcSchema["slot_Centroid"])
        self.srcCentroidKey_xErr = srcSchema["slot_Centroid_xErr"].asKey()
        self.srcCentroidKey_yErr = srcSchema["slot_Centroid_yErr"].asKey()
        self.sourceCat = afwTable.SourceCatalog(srcSchema)

        self.matches = []

        for i in np.linspace(0., rangePix, numPoints):
            for j in np.linspace(0., rangePix, numPoints):
                src = self.sourceCat.addNew()
                refObj = self.refCat.addNew()

                src.set(self.srcCentroidKey, lsst.geom.Point2D(i, j))
                src.set(self.srcCentroidKey_xErr, 0.1)
                src.set(self.srcCentroidKey_yErr, 0.1)

                c = self.tanWcs.pixelToSky(i, j)
                refObj.setCoord(c)

                self.matches.append(self.MatchClass(refObj, src, 0.0))

    def tearDown(self):
        del self.refCat
        del self.sourceCat
        del self.matches
        del self.tanWcs

    def checkResults(self, fitRes, catsUpdated):
        """Check results.
        """
        self.assertLess(fitRes.scatterOnSky.asArcseconds(), 0.001)
        affineWcs = fitRes.wcs
        maxAngSep = 0*lsst.geom.radians
        maxPixSep = 0
        refCoordKey = afwTable.CoordKey(self.refCat.schema["coord"])
        if catsUpdated:
            refCentroidKey = afwTable.Point2DKey(self.refCat.schema["centroid"])
        maxDistErr = 0*lsst.geom.radians
        for refObj, src, distRad in self.matches:
            srcPixPos = src.get(self.srcCentroidKey)
            refCoord = refObj.get(refCoordKey)
            if catsUpdated:
                refPixPos = refObj.get(refCentroidKey)
                srcCoord = src.get(self.srcCoordKey)
            else:
                refPixPos = affineWcs.skyToPixel(refCoord)
                srcCoord = affineWcs.pixelToSky(srcPixPos)

            angSep = refCoord.separation(srcCoord)
            dist = distRad*lsst.geom.radians
            distErr = abs(dist - angSep)
            maxDistErr = max(maxDistErr, distErr)
            maxAngSep = max(maxAngSep, angSep)

            pixSep = math.hypot(*(srcPixPos - refPixPos))
            maxPixSep = max(maxPixSep, pixSep)

        print("max angular separation = %0.4f arcsec" % (maxAngSep.asArcseconds(),))
        print("max pixel separation = %0.3f" % (maxPixSep,))
        self.assertLess(maxAngSep.asArcseconds(), 0.001)
        self.assertLess(maxPixSep, 0.005)
        if catsUpdated:
            allowedDistErr = 1e-7
        else:
            allowedDistErr = 0.001
        self.assertLess(maxDistErr.asArcseconds(), allowedDistErr,
                        "Computed distance in match list is off by %s arcsec" % (maxDistErr.asArcseconds(),))

    def doTest(self, name, func):
        """Apply func(x, y) to each source in self.sourceCat, then fit and
        check the resulting WCS.
        """
        bbox = lsst.geom.Box2I()
        for refObj, src, d in self.matches:
            origPos = src.get(self.srcCentroidKey)
            x, y = func(*origPos)
            distortedPos = lsst.geom.Point2D(x, y)
            src.set(self.srcCentroidKey, distortedPos)
            bbox.include(lsst.geom.Point2I(lsst.geom.Point2I(distortedPos)))

        fitter = FitAffineWcsTask()
        fitRes = fitter.fitWcs(
            matches=self.matches,
            initWcs=self.tanWcs,
            bbox=bbox,
            refCat=self.refCat,
            sourceCat=self.sourceCat,
        )

        self.checkResults(fitRes, catsUpdated=True)

    def doTestAffine(self, name, offset, matrix):
        """Apply func(x, y) to each source in self.sourceCat, then fit and
        check the resulting WCS.
        """
        wcsMaker = TransformedSkyWcsMaker(self.tanWcs)

        newWcs = wcsMaker.makeWcs(offset, matrix)
        bbox = lsst.geom.Box2I()
        for refObj, src, d in self.matches:
            origPos = src.get(self.srcCentroidKey)
            newCoord = newWcs.pixelToSky(origPos)
            src.setCoord(newCoord)
            bbox.include(lsst.geom.Point2I(lsst.geom.Point2I(origPos)))

        fitter = FitAffineWcsTask()
        fitRes = fitter.fitWcs(
            matches=self.matches,
            initWcs=newWcs,
            bbox=bbox,
            refCat=self.refCat,
            sourceCat=self.sourceCat,
        )

        self.checkResults(fitRes, catsUpdated=True)

    def doTestAffineReverse(self, name, offset, matrix):
        """Apply func(x, y) to each source in self.sourceCat, then fit and
        check the resulting WCS.
        """
        wcsMaker = TransformedSkyWcsMaker(self.tanWcs)

        newWcs = wcsMaker.makeWcs(offset, matrix)
        bbox = lsst.geom.Box2I()
        for refObj, src, d in self.matches:
            origCoord = src.get(self.srcCoordKey)
            newCentroid = newWcs.skyToPixel(origCoord)
            src.set(self.srcCentroidKey, newCentroid)
            bbox.include(lsst.geom.Point2I(lsst.geom.Point2I(newCentroid)))

        fitter = FitAffineWcsTask()
        fitRes = fitter.fitWcs(
            matches=self.matches,
            initWcs=self.tanWcs,
            bbox=bbox,
            refCat=self.refCat,
            sourceCat=self.sourceCat,
        )

        self.checkResults(fitRes, catsUpdated=True)


class SideLoadTestCases:

    """Base class implementations of testing methods.

    Explicitly does not inherit from unittest.TestCase"""

    def testTrivial(self):
        """Add no distortion"""
        self.doTest("testTrivial", lambda x, y: (x, y))

    def testOffset(self):
        """Add an offset"""
        self.doTest("testOffset", lambda x, y: (x + 5, y + 7))

    def testSkyOffset(self):
        """Add an on sky offset"""
        self.doTestAffine("testSkyOffset",
                          [77, -200],
                          np.array([[1.0, 0.0], [0.0, 1.0]]))
        self.doTestAffineReverse("testSkyOffsetRev",
                                 [77, -200],
                                 np.array([[1.0, 0.0], [0.0, 1.0]]))

    def testAffine(self):
        """Add an Affine (shear + scale + rot) distortion"""
        self.doTestAffine("testAffine",
                          [0, 0],
                          np.array([[0.4, 0.1], [-0.21, 2.0]]))
        self.doTestAffineReverse("testAffineRev",
                                 [0, 0],
                                 np.array([[0.4, 0.1], [-0.21, 2.0]]))

    def testAffineAndOffset(self):
        """Add a transform and offset"""
        self.doTestAffine("testAffineAndOffset",
                          [30, 100],
                          np.array([[0.5, 0.01], [-0.2, 0.3]]))
        self.doTestAffineReverse("testAffineAndOffsetRev",
                                 [30, 100],
                                 np.array([[0.5, 0.01], [-0.2, 0.3]]))


# The test classes inherit from two base classes and differ in the match
# class being used.

class FitAffineWcsTaskTestCaseReferenceMatch(BaseTestCase,
                                             SideLoadTestCases,
                                             lsst.utils.tests.TestCase):
    MatchClass = afwTable.ReferenceMatch


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
