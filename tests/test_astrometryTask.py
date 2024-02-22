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

import astropy.units as u
import scipy.stats
import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
from lsst.meas.algorithms.testUtils import MockReferenceObjectLoaderFromFiles
from lsst.meas.astrom import AstrometryTask, exceptions
import lsst.pipe.base as pipeBase


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
        self.exposure.info.setVisitInfo(afwImage.VisitInfo(date=lsst.daf.base.DateTime(60000)))
        self.exposure.setFilter(afwImage.FilterLabel(band="r", physical="rTest"))
        filenames = sorted(glob.glob(os.path.join(refCatDir, 'ref_cats', 'cal_ref_cat', '??????.fits')))
        self.refObjLoader = MockReferenceObjectLoaderFromFiles(filenames, htmLevel=8)

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
        schema = self._makeSourceCatalogSchema()
        # schema must be passed to the solver task constructor
        solver = AstrometryTask(config=config, refObjLoader=self.refObjLoader, schema=schema)
        sourceCat = self.makeSourceCat(self.tanWcs, schema)

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


    def doTest(self, pixelsToTanPixels, order=3):
        """Test using pixelsToTanPixels to distort the source positions
        """
        distortedWcs = afwGeom.makeModifiedWcs(pixelTransform=pixelsToTanPixels, wcs=self.tanWcs,
                                               modifyActualPixels=False)
        self.exposure.setWcs(distortedWcs)
        sourceCat = self.makeSourceCat(distortedWcs, self._makeSourceCatalogSchema())
        config = AstrometryTask.ConfigClass()
        config.wcsFitter.order = order
        config.wcsFitter.numRejIter = 0
        # This test is from before rough magnitude rejection was implemented.
        config.doMagnitudeOutlierRejection = False
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

    @staticmethod
    def _makeSourceCatalogSchema():
        """Return a catalog schema with all necessary fields added.
        """
        schema = afwTable.SourceTable.makeMinimalSchema()
        measBase.SingleFrameMeasurementTask(schema=schema)  # expand the schema
        afwTable.CoordKey.addErrorFields(schema)
        schema.addField("deblend_nChild", type=np.int32,
                        doc="Number of children this object has (defaults to 0)")
        schema.addField("detect_isPrimary", type=np.int32,
                        doc="true if source has no children and is not a sky source")
        return schema

    def makeSourceCat(self, wcs, schema, doScatterCentroids=False):
        """Make a source catalog by reading the position reference stars using
        the proviced WCS.

        Optionally, via doScatterCentroids, add some scatter to the centroids
        assiged to the source catalog (otherwise they will be identical to
        those of the reference catalog).
        """
        loadRes = self.refObjLoader.loadPixelBox(bbox=self.bbox, wcs=wcs, filterName="r")
        refCat = loadRes.refCat

        sourceCat = afwTable.SourceCatalog(schema)

        sourceCat.resize(len(refCat))
        scatterFactor = 1.0
        if doScatterCentroids:
            np.random.seed(12345)
            scatterFactor = np.random.uniform(0.999, 1.001, len(sourceCat))
        sourceCat["slot_Centroid_x"] = scatterFactor*refCat["centroid_x"]
        sourceCat["slot_Centroid_y"] = scatterFactor*refCat["centroid_y"]
        sourceCat["slot_PsfFlux_instFlux"] = refCat["r_flux"]
        sourceCat["slot_PsfFlux_instFluxErr"] = refCat["r_flux"]/100
        # All of these sources are primary.
        sourceCat['detect_isPrimary'] = 1

        # Deliberately add some outliers to check that the magnitude
        # outlier rejection code is being run.
        sourceCat["slot_PsfFlux_instFlux"][0: 4] *= 1000.0

        return sourceCat

    def testBadAstrometry(self):
        """Test that an appropriately informative exception is raised for a
        bad quality fit.
        """
        catalog = self.makeSourceCat(self.tanWcs, self._makeSourceCatalogSchema())
        # Fake match list with 10" match distance for every source to force
        # a bad quality fit result.
        matches = [afwTable.ReferenceMatch(r, r, (10*u.arcsecond).to_value(u.radian)) for r in catalog]
        result = pipeBase.Struct(
            matches=matches,
            wcs=None,
            scatterOnSky=20.0*lsst.geom.arcseconds,
            matchTolerance=None
        )
        with unittest.mock.patch("lsst.meas.astrom.AstrometryTask._matchAndFitWcs",
                                 return_value=result, autospec=True):
            with self.assertRaises(exceptions.BadAstrometryFit):
                task = AstrometryTask(refObjLoader=self.refObjLoader)
                task.run(catalog, self.exposure)
            self.assertIsNone(self.exposure.wcs)
            self.assertTrue(np.all(np.isnan(catalog["coord_ra"])))
            self.assertTrue(np.all(np.isnan(catalog["coord_dec"])))

    def testMatcherFails(self):
        """Test that a matcher exception has additional metadata attached.
        """
        catalog = self.makeSourceCat(self.tanWcs, self._makeSourceCatalogSchema())
        with unittest.mock.patch("lsst.meas.astrom.AstrometryTask._matchAndFitWcs", autospec=True,
                                 side_effect=exceptions.MatcherFailure("some matcher problem")):
            with self.assertRaises(exceptions.MatcherFailure) as cm:
                task = AstrometryTask(refObjLoader=self.refObjLoader)
                task.run(catalog, self.exposure)
            self.assertEqual(cm.exception.metadata["iterations"], 1)
            self.assertIsNone(self.exposure.wcs)
            self.assertTrue(np.all(np.isnan(catalog["coord_ra"])))
            self.assertTrue(np.all(np.isnan(catalog["coord_dec"])))

    def testExceptions(self):
        """Test that the custom astrometry exceptions are well behaved.
        """
        error = exceptions.AstrometryError("something", blah=10)
        self.assertEqual(error.metadata["blah"], 10)
        self.assertIn("something", str(error))
        self.assertIn("'blah': 10", str(error))

        # Metadata cannot contain an astropy unit.
        error = exceptions.AstrometryError("something", blah=10*u.arcsecond)
        with self.assertRaisesRegex(TypeError, "blah is of type <class 'astropy.units.quantity.Quantity'>"):
            error.metadata


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

        refCat['refFlux'] = (refMag*u.ABmag).to_value(u.nJy)
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
