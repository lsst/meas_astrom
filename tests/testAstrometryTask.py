#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

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
import math
import os
import unittest

import numpy as np

import eups
import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
import lsst.meas.astrom as measAstrom

class TestAstrometricSolver(unittest.TestCase):

    def setUp(self):
        mypath = eups.productDir("meas_astrom")
        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))
        self.datapath = datapath

        self.bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(3001, 3001))
        self.ctrPix = afwGeom.Point2I(1500, 1500)
        metadata = dafBase.PropertySet()
        metadata.set("RADECSYS", "FK5")
        metadata.set("EQUINOX", 2000.0)
        metadata.set("CTYPE1", "RA---TAN")
        metadata.set("CTYPE2", "DEC--TAN")
        metadata.set("CUNIT1", "deg")
        metadata.set("CUNIT2", "deg")
        metadata.set("CRVAL1", 215.5)
        metadata.set("CRVAL2",  53.0)
        metadata.set("CRPIX1", self.ctrPix[0] + 1)
        metadata.set("CRPIX2", self.ctrPix[1] + 1)
        metadata.set("CD1_1",  5.1e-05)
        metadata.set("CD1_2",  0.0)
        metadata.set("CD2_2", -5.1e-05)
        metadata.set("CD2_1",  0.0)
        self.tanWcs = afwImage.cast_TanWcs(afwImage.makeWcs(metadata))

    def tearDown(self):
        del self.ctrPix
        del self.tanWcs

    def testTrivial(self):
        """Test fit with no distortion
        """
        self.doTest(afwGeom.IdentityXYTransform())

    def testRadial(self):
        """Test fit with radial distortion

        The offset comes from the fact that the CCD is not centered
        """
        self.doTest(afwGeom.RadialXYTransform([0, 1.01, 1e-7]))

    def doTest(self, pixelsToTanPixels, order=3):
        """Test using pixelsToTanPixels to distort the source positions
        """
        distortedWcs = afwImage.DistortedTanWcs(self.tanWcs, pixelsToTanPixels)
        sourceCat = self.makeSourceCat(distortedWcs)
        config = measAstrom.AstrometryTask.ConfigClass()
        config.wcsFitter.order = order
        solver = measAstrom.AstrometryTask(config=config)
        results = solver.solve(
            sourceCat = sourceCat,
            bbox = self.bbox,
            initWcs = distortedWcs,
            filterName = 'r'
        )
        self.assertTrue(results.initWcs is distortedWcs)
        self.assertFalse(results.wcs is distortedWcs)
        self.assertWcssAlmostEqual(distortedWcs, results.wcs, self.bbox)

        srcCoordKey = sourceCat.schema["coord"].asKey()
        refCoordKey = results.refCat.schema["coord"].asKey()
        refCentroidKey = afwTable.Point2DKey(results.refCat.schema["centroid"])
        maxAngSep = afwGeom.Angle(0)
        maxPixSep = 0
        for refObj, src, d in results.matches:
            refCoord = refObj.get(refCoordKey)
            refPixPos = refObj.get(refCentroidKey)
            srcCoord = src.get(srcCoordKey)
            srcPixPos = src.getCentroid()

            angSep = refCoord.angularSeparation(srcCoord)
            maxAngSep = max(maxAngSep, angSep)
            self.assertLess(refCoord.angularSeparation(srcCoord).asArcseconds(), 0.0025)

            pixSep = math.hypot(*(srcPixPos-refPixPos))
            maxPixSep = max(maxPixSep, pixSep)
            self.assertLess(pixSep, 0.015)
        print("max angular separation = %0.4f arcsec" % (maxAngSep.asArcseconds(),))
        print("max pixel separation = %0.3f" % (maxPixSep,))

        # try again, but without fitting the WCS
        config.forceKnownWcs = True
        solverNoFit = measAstrom.AstrometryTask(config=config)
        resultsNoFit = solverNoFit.solve(
            sourceCat = sourceCat,
            bbox = self.bbox,
            initWcs = distortedWcs,
            filterName = 'r'
        )
        self.assertTrue(resultsNoFit.wcs is distortedWcs)
        self.assertTrue(resultsNoFit.initWcs is distortedWcs)
        self.assertTrue(resultsNoFit.scatterOnSky is None)

        # fitting should improve the quality of the matches
        meanFitDist = np.mean([match.distance for match in results.matches])
        meanNoFitDist = np.mean([match.distance for match in results.matches])
        self.assertLessEqual(meanFitDist, meanNoFitDist)

    def makeSourceCat(self, distortedWcs):
        """Make a source catalog by reading the position reference stars and distorting the positions
        """
        loaderConfig = measAstrom.LoadAstrometryNetObjectsTask.ConfigClass()
        loader = measAstrom.LoadAstrometryNetObjectsTask(config=loaderConfig)
        loadRes = loader.loadPixelBox(bbox=self.bbox, wcs=distortedWcs, filterName="r")
        refCat = loadRes.refCat
        refCentroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
        refFluxRKey = refCat.schema["r_flux"].asKey()

        sourceSchema = afwTable.SourceTable.makeMinimalSchema()
        measBase.SingleFrameMeasurementTask(schema=sourceSchema) # expand the schema
        sourceCat = afwTable.SourceCatalog(sourceSchema)
        sourceCentroidKey = afwTable.Point2DKey(sourceSchema["slot_Centroid"])
        sourceFluxKey = sourceSchema["slot_ApFlux_flux"].asKey()

        for refObj in refCat:
            src = sourceCat.addNew()
            src.set(sourceCentroidKey, refObj.get(refCentroidKey))
            src.set(sourceFluxKey, refObj.get(refFluxRKey))
        return sourceCat

    def assertWcssAlmostEqual(self, wcs0, wcs1, bbox, maxSkyErr=0.01 * afwGeom.arcseconds, maxPixErr = 0.02,
        nx = 5, ny = 5, msg="WCSs differ"):
        """Assert that two WCSs give nearly equal results for pixelToSky and skyToPixel

        @param[in] wcs0  WCS 0 (an lsst.afw.image.Wcs)
        @param[in] wcs1  WCS 1 (an lsst.afw.image.Wcs)
        @param[in] bbox  boundaries of pixel grid over which to compare the WCSs (an lsst.afw.geom.Box2I or Box2D)
        @param[in] maxSkyErr  maximum separation between sky positions computed using Wcs.pixelToSky
            (an lsst.afw.geom.Angle)
        @param[in] maxPixErr  maximum separation between pixel positions computed using Wcs.skyToPixel
        @param[in] nx  number of points in x for the grid of pixel positions
        @param[in] ny  number of points in y for the grid of pixel positions
        @param[in] msg  prefix for error messages; details of the error are appended after ": "

        @throw AssertionError if the two WCSs do not match sufficiently closely
        """
        bboxd = afwGeom.Box2D(bbox)
        xList = np.linspace(bboxd.getMinX(), bboxd.getMaxX(), nx)
        yList = np.linspace(bboxd.getMinY(), bboxd.getMaxY(), ny)
        maxSkyErrPixPos = [afwGeom.Angle(0), None]
        maxPixErrSkyPos = [0, None]
        for x in xList:
            for y in yList:
                fromPixPos = afwGeom.Point2D(x, y)
                sky0 = wcs0.pixelToSky(fromPixPos)
                sky1 = wcs1.pixelToSky(fromPixPos)
                skyErr = sky0.angularSeparation(sky1)
                if skyErr > maxSkyErrPixPos[0]:
                    maxSkyErrPixPos = (skyErr, fromPixPos)

                toPixPos0 = wcs0.skyToPixel(sky0)
                toPixPos1 = wcs1.skyToPixel(sky0)
                pixErr = math.hypot(*(toPixPos0 - toPixPos1))
                if pixErr > maxPixErrSkyPos[0]:
                    maxPixErrSkyPos = (pixErr, sky0)

        errStrList = []
        if maxSkyErrPixPos[0] > maxSkyErr:
            skyErr, pixPos = maxSkyErrPixPos
            errStrList.append("%f arcsec sky error > %f arcsec max sky error at pixPos=%s" %
                (skyErr.asArcseconds(), maxSkyErr.asArcseconds(), pixPos))
        if maxPixErrSkyPos[0] > maxPixErr:
            pixErr, skyPos = maxPixErrSkyPos
            errStrList.append("%f pix error > %f max pix error at skyPos=%s" %
                    (pixErr, maxPixErr, skyPos))
        if errStrList:
            raise AssertionError("%s: %s" % (msg, "; ".join(errStrList)))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(TestAstrometricSolver)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
