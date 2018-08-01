#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
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

import os.path
import math
import unittest

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
from lsst.daf.persistence import Butler
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
from lsst.meas.astrom import AstrometryTask


class TestAstrometricSolver(lsst.utils.tests.TestCase):

    def setUp(self):
        refCatDir = os.path.join(os.path.dirname(__file__), "data", "sdssrefcat")

        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(3001, 3001))
        crpix = lsst.geom.Box2D(self.bbox).getCenter()
        self.tanWcs = afwGeom.makeSkyWcs(crpix=crpix,
                                         crval=lsst.geom.SpherePoint(215.5, 53.0, lsst.geom.degrees),
                                         cdMatrix=afwGeom.makeCdMatrix(scale=5.1e-5*lsst.geom.degrees))
        self.exposure = afwImage.ExposureF(self.bbox)
        self.exposure.setWcs(self.tanWcs)
        self.exposure.setFilter(afwImage.Filter("r", True))
        butler = Butler(refCatDir)
        self.refObjLoader = LoadIndexedReferenceObjectsTask(butler=butler)

    def tearDown(self):
        del self.tanWcs
        del self.exposure
        del self.refObjLoader

    def testTrivial(self):
        """Test fit with no distortion
        """
        self.doTest(afwGeom.makeIdentityTransform())

    def testRadial(self):
        """Test fit with radial distortion

        The offset comes from the fact that the CCD is not centered
        """
        self.doTest(afwGeom.makeRadialTransform([0, 1.01, 1e-7]))

    def testUsedFlag(self):
        """Test that the solver will record number of sources used to table
           if it is passed a schema on initialization.
        """
        self.exposure.setWcs(self.tanWcs)
        loadRes = self.refObjLoader.loadPixelBox(bbox=self.bbox, wcs=self.tanWcs, filterName="r")
        refCat = loadRes.refCat
        refCentroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
        refFluxRKey = refCat.schema["r_flux"].asKey()

        sourceSchema = afwTable.SourceTable.makeMinimalSchema()
        measBase.SingleFrameMeasurementTask(schema=sourceSchema)  # expand the schema
        config = AstrometryTask.ConfigClass()
        config.wcsFitter.order = 2
        config.wcsFitter.numRejIter = 0
        # schema must be passed to the solver task constructor
        solver = AstrometryTask(config=config, refObjLoader=self.refObjLoader, schema=sourceSchema)
        sourceCat = afwTable.SourceCatalog(sourceSchema)
        sourceCat.reserve(len(refCat))
        sourceCentroidKey = afwTable.Point2DKey(sourceSchema["slot_Centroid"])
        sourceFluxKey = sourceSchema["slot_ApFlux_flux"].asKey()
        sourceFluxErrKey = sourceSchema["slot_ApFlux_fluxErr"].asKey()

        for refObj in refCat:
            src = sourceCat.addNew()
            src.set(sourceCentroidKey, refObj.get(refCentroidKey))
            src.set(sourceFluxKey, refObj.get(refFluxRKey))
            src.set(sourceFluxErrKey, refObj.get(refFluxRKey)/100)

        results = solver.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        # check that the used flag is set the right number of times
        count = 0
        for source in sourceCat:
            if source.get('calib_astrometryUsed'):
                count += 1
        self.assertEqual(count, len(results.matches))

    def doTest(self, pixelsToTanPixels, order=3):
        """Test using pixelsToTanPixels to distort the source positions
        """
        distortedWcs = afwGeom.makeModifiedWcs(pixelTransform=pixelsToTanPixels, wcs=self.tanWcs,
                                               modifyActualPixels=False)
        self.exposure.setWcs(distortedWcs)
        sourceCat = self.makeSourceCat(distortedWcs)
        config = AstrometryTask.ConfigClass()
        config.wcsFitter.order = order
        config.wcsFitter.numRejIter = 0
        solver = AstrometryTask(config=config, refObjLoader=self.refObjLoader)
        results = solver.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        fitWcs = self.exposure.getWcs()
        self.assertRaises(Exception, self.assertWcsAlmostEqualOverBBox, fitWcs, distortedWcs)
        self.assertWcsAlmostEqualOverBBox(distortedWcs, fitWcs, self.bbox,
                                          maxDiffSky=0.01*lsst.geom.arcseconds, maxDiffPix=0.02)

        srcCoordKey = afwTable.CoordKey(sourceCat.schema["coord"])
        refCoordKey = afwTable.CoordKey(results.refCat.schema["coord"])
        refCentroidKey = afwTable.Point2DKey(results.refCat.schema["centroid"])
        maxAngSep = 0*lsst.geom.radians
        maxPixSep = 0
        for refObj, src, d in results.matches:
            refCoord = refObj.get(refCoordKey)
            refPixPos = refObj.get(refCentroidKey)
            srcCoord = src.get(srcCoordKey)
            srcPixPos = src.getCentroid()

            angSep = refCoord.separation(srcCoord)
            maxAngSep = max(maxAngSep, angSep)

            pixSep = math.hypot(*(srcPixPos-refPixPos))
            maxPixSep = max(maxPixSep, pixSep)
        print("max angular separation = %0.4f arcsec" % (maxAngSep.asArcseconds(),))
        print("max pixel separation = %0.3f" % (maxPixSep,))
        self.assertLess(maxAngSep.asArcseconds(), 0.0026)
        self.assertLess(maxPixSep, 0.015)

        # try again, but without fitting the WCS
        config.forceKnownWcs = True
        solverNoFit = AstrometryTask(config=config, refObjLoader=self.refObjLoader)
        self.exposure.setWcs(distortedWcs)
        resultsNoFit = solverNoFit.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        self.assertIsNone(resultsNoFit.scatterOnSky)

        # fitting should result in matches that are at least as good
        # (strictly speaking fitting might result in a larger match list with
        # some outliers, but in practice this test passes)
        meanFitDist = np.mean([match.distance for match in results.matches])
        meanNoFitDist = np.mean([match.distance for match in resultsNoFit.matches])
        self.assertLessEqual(meanFitDist, meanNoFitDist)

    def makeSourceCat(self, distortedWcs):
        """Make a source catalog by reading the position reference stars and distorting the positions
        """
        loadRes = self.refObjLoader.loadPixelBox(bbox=self.bbox, wcs=distortedWcs, filterName="r")
        refCat = loadRes.refCat
        refCentroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
        refFluxRKey = refCat.schema["r_flux"].asKey()

        sourceSchema = afwTable.SourceTable.makeMinimalSchema()
        measBase.SingleFrameMeasurementTask(schema=sourceSchema)  # expand the schema
        sourceCat = afwTable.SourceCatalog(sourceSchema)
        sourceCentroidKey = afwTable.Point2DKey(sourceSchema["slot_Centroid"])
        sourceFluxKey = sourceSchema["slot_ApFlux_flux"].asKey()
        sourceFluxErrKey = sourceSchema["slot_ApFlux_fluxErr"].asKey()

        sourceCat.reserve(len(refCat))
        for refObj in refCat:
            src = sourceCat.addNew()
            src.set(sourceCentroidKey, refObj.get(refCentroidKey))
            src.set(sourceFluxKey, refObj.get(refFluxRKey))
            src.set(sourceFluxErrKey, refObj.get(refFluxRKey)/100)
        return sourceCat


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
