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
import glob

from astropy import units
import scipy.stats
import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms.testUtils import MockLoadReferenceObjects
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
        self.exposure.setFilterLabel(afwImage.FilterLabel(band="r", physical="rTest"))
        filenames = sorted(glob.glob(os.path.join(refCatDir, 'ref_cats', 'cal_ref_cat', '??????.fits')))
        self.refObjLoader = MockLoadReferenceObjects(filenames, convert=True)

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
        config = AstrometryTask.ConfigClass()
        config.wcsFitter.order = 2
        config.wcsFitter.numRejIter = 0

        sourceSchema = afwTable.SourceTable.makeMinimalSchema()
        measBase.SingleFrameMeasurementTask(schema=sourceSchema)  # expand the schema
        # schema must be passed to the solver task constructor
        solver = AstrometryTask(config=config, refObjLoader=self.refObjLoader, schema=sourceSchema)
        sourceCat = self.makeSourceCat(self.tanWcs, sourceSchema=sourceSchema)

        results = solver.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        # check that the used flag is set the right number of times
        count = 0
        for source in sourceCat:
            if source.get('calib_astrometry_used'):
                count += 1
        self.assertEqual(count, len(results.matches))

    def testMaxMeanDistance(self):
        """If the astrometric fit does not satisfy the maxMeanDistanceArcsec
        threshold, ensure task raises an lsst.pipe.base.TaskError.
        """
        self.exposure.setWcs(self.tanWcs)
        config = AstrometryTask.ConfigClass()
        config.maxMeanDistanceArcsec = 0.0  # To ensure a "deemed" WCS failure
        solver = AstrometryTask(config=config, refObjLoader=self.refObjLoader)
        sourceCat = self.makeSourceCat(self.tanWcs, doScatterCentroids=True)

        with self.assertRaisesRegex(pipeBase.TaskError, "Fatal astrometry failure detected"):
            solver.run(sourceCat=sourceCat, exposure=self.exposure)

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
        self.assertLess(maxAngSep.asArcseconds(), 0.0038)
        self.assertLess(maxPixSep, 0.021)

        # try again, invoking the reference selector
        config.referenceSelector.doUnresolved = True
        config.referenceSelector.unresolved.name = 'resolved'
        solverRefSelect = AstrometryTask(config=config, refObjLoader=self.refObjLoader)
        self.exposure.setWcs(distortedWcs)
        resultsRefSelect = solverRefSelect.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        self.assertLess(len(resultsRefSelect.matches), len(results.matches))

        # try again, allowing magnitude outlier rejection.
        config.doMagnitudeOutlierRejection = True
        solverMagOutlierRejection = AstrometryTask(config=config, refObjLoader=self.refObjLoader)
        self.exposure.setWcs(distortedWcs)
        resultsMagOutlierRejection = solverMagOutlierRejection.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        self.assertLess(len(resultsMagOutlierRejection.matches), len(resultsRefSelect.matches))
        config.doMagnitudeOutlierRejection = False

        # try again, but without fitting the WCS, no reference selector
        config.referenceSelector.doUnresolved = False
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

        # try once again, without fitting the WCS, with the reference selector
        # (this goes through a different code path)
        config.referenceSelector.doUnresolved = True
        solverNoFitRefSelect = AstrometryTask(config=config, refObjLoader=self.refObjLoader)
        resultsNoFitRefSelect = solverNoFitRefSelect.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        self.assertLess(len(resultsNoFitRefSelect.matches), len(resultsNoFit.matches))

    def makeSourceCat(self, wcs, sourceSchema=None, doScatterCentroids=False):
        """Make a source catalog by reading the position reference stars using
        the proviced WCS.

        Optionally provide a schema for the source catalog (to allow
        AstrometryTask in the test methods to update it with the
        "calib_astrometry_used" flag).  Otherwise, a minimal SourceTable
        schema will be created.

        Optionally, via doScatterCentroids, add some scatter to the centroids
        assiged to the source catalog (otherwise they will be identical to
        those of the reference catalog).
        """
        loadRes = self.refObjLoader.loadPixelBox(bbox=self.bbox, wcs=wcs, filterName="r")
        refCat = loadRes.refCat

        if sourceSchema is None:
            sourceSchema = afwTable.SourceTable.makeMinimalSchema()
            measBase.SingleFrameMeasurementTask(schema=sourceSchema)  # expand the schema
        sourceCat = afwTable.SourceCatalog(sourceSchema)

        sourceCat.resize(len(refCat))
        scatterFactor = 1.0
        if doScatterCentroids:
            np.random.seed(12345)
            scatterFactor = np.random.uniform(0.999, 1.001, len(sourceCat))
        sourceCat["slot_Centroid_x"] = scatterFactor*refCat["centroid_x"]
        sourceCat["slot_Centroid_y"] = scatterFactor*refCat["centroid_y"]
        sourceCat["slot_ApFlux_instFlux"] = refCat["r_flux"]
        sourceCat["slot_ApFlux_instFluxErr"] = refCat["r_flux"]/100

        # Deliberately add some outliers to check that the magnitude
        # outlier rejection code is being run.
        sourceCat["slot_ApFlux_instFlux"][0: 4] *= 1000.0

        return sourceCat


class TestMagnitudeOutliers(lsst.utils.tests.TestCase):
    def testMagnitudeOutlierRejection(self):
        """Test rejection of magnitude outliers.

        This test only tests the outlier rejection, and not any other
        part of the matching or astrometry fitter.
        """
        config = AstrometryTask.ConfigClass()
        config.doMagnitudeOutlierRejection = True
        config.magnitudeOutlierRejectionNSigma = 4.0
        solver = AstrometryTask(config=config, refObjLoader=None)

        nTest = 100

        refSchema = lsst.afw.table.SimpleTable.makeMinimalSchema()
        refSchema.addField('refFlux', 'F')
        refCat = lsst.afw.table.SimpleCatalog(refSchema)
        refCat.resize(nTest)

        srcSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        srcSchema.addField('srcFlux', 'F')
        srcCat = lsst.afw.table.SourceCatalog(srcSchema)
        srcCat.resize(nTest)

        np.random.seed(12345)
        refMag = np.full(nTest, 20.0)
        srcMag = np.random.normal(size=nTest, loc=0.0, scale=1.0)

        # Determine the sigma of the random sample
        zp = np.median(refMag[: -4] - srcMag[: -4])
        sigma = scipy.stats.median_abs_deviation(srcMag[: -4], scale='normal')

        # Deliberately alter some magnitudes to be outliers.
        srcMag[-3] = (config.magnitudeOutlierRejectionNSigma + 0.1)*sigma + (20.0 - zp)
        srcMag[-4] = -(config.magnitudeOutlierRejectionNSigma + 0.1)*sigma + (20.0 - zp)

        refCat['refFlux'] = (refMag*units.ABmag).to_value(units.nJy)
        srcCat['srcFlux'] = 10.0**(srcMag/(-2.5))

        # Deliberately poison some reference fluxes.
        refCat['refFlux'][-1] = np.inf
        refCat['refFlux'][-2] = np.nan

        matchesIn = []
        for ref, src in zip(refCat, srcCat):
            matchesIn.append(lsst.afw.table.ReferenceMatch(first=ref, second=src, distance=0.0))

        matchesOut = solver._removeMagnitudeOutliers('srcFlux', 'refFlux', matchesIn)

        # We should lose the 4 outliers we created.
        self.assertEqual(len(matchesOut), len(matchesIn) - 4)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
