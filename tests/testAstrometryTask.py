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

import sys
import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
import lsst.meas.astrom as measAstrom

import testFindAstrometryNetDataDir as helper


class TestAstrometricSolver(utilsTests.TestCase):

    def setUp(self):
        # Set up local astrometry_net_data
        self.datapath = helper.setupAstrometryNetDataDir('photocal')

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
        self.exposure = afwImage.ExposureF(self.bbox)
        self.exposure.setWcs(self.tanWcs)
        self.exposure.setFilter(afwImage.Filter("r", True))

    def tearDown(self):
        del self.ctrPix
        del self.tanWcs
        del self.exposure

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
        self.exposure.setWcs(distortedWcs)
        sourceCat = self.makeSourceCat(distortedWcs)
        config = measAstrom.AstrometryTask.ConfigClass()
        config.wcsFitter.order = order
        solver = measAstrom.AstrometryTask(config=config)
        results = solver.run(
            sourceCat = sourceCat,
            exposure = self.exposure,
        )
        fitWcs = self.exposure.getWcs()
        self.assertRaises(Exception, self.assertWcsNearlyEqualOverBBox, fitWcs, distortedWcs)
        self.assertWcsNearlyEqualOverBBox(distortedWcs, fitWcs, self.bbox,
            maxDiffSky=0.01*afwGeom.arcseconds, maxDiffPix=0.02)

        srcCoordKey = afwTable.CoordKey(sourceCat.schema["coord"])
        refCoordKey = afwTable.CoordKey(results.refCat.schema["coord"])
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
        self.exposure.setWcs(distortedWcs)
        resultsNoFit = solverNoFit.run(
            sourceCat = sourceCat,
            exposure = self.exposure,
        )
        self.assertTrue(resultsNoFit.scatterOnSky is None)

        # fitting should result in matches that are at least as good
        # (strictly speaking fitting might result in a larger match list with
        # some outliers, but in practice this test passes)
        meanFitDist = np.mean([match.distance for match in results.matches])
        meanNoFitDist = np.mean([match.distance for match in resultsNoFit.matches])
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
