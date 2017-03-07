#!/usr/bin/env python
from __future__ import division
import unittest

import numpy as np

import lsst.utils.tests
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.meas.astrom import approximateWcs


class ApproximateWcsTestCase(lsst.utils.tests.TestCase):

    """A test case for CreateWcsWithSip

    Use involves setting one class attribute:
    * MatchClass: match class, e.g. ReferenceMatch or SourceMatch
    """

    def setUp(self):
        metadata = dafBase.PropertySet()

        self.crPix = afwGeom.Point2D(15000, 4000)
        dimd = afwGeom.Extent2D(4000, 4000)
        bboxd = afwGeom.Box2D(self.crPix - dimd/2, dimd)
        self.bbox = afwGeom.Box2I(bboxd)
        metadata.set("RADECSYS", 'ICRS')
        metadata.set("EQUINOX", 2000.0)
        metadata.setDouble("CRVAL1", 215.60)
        metadata.setDouble("CRVAL2", 53.16)
        metadata.setDouble("CRPIX1", self.crPix[0])
        metadata.setDouble("CRPIX2", self.crPix[1])
        metadata.set("CTYPE1", "RA---TAN")
        metadata.set("CTYPE2", "DEC--TAN")
        metadata.setDouble("CD1_1", 5.10808596133527E-05)
        metadata.setDouble("CD1_2", 1.85579539217196E-07)
        metadata.setDouble("CD2_2", -5.10281493481982E-05)
        metadata.setDouble("CD2_1", -8.27440751733828E-07)
        self.tanWcs = afwImage.makeWcs(metadata)

    def tearDown(self):
        del self.tanWcs

    def testTrivial(self):
        """Add no distortion"""
        for order in (3, 4, 5, 6):
            self.doTest("testTrivial", afwGeom.IdentityXYTransform(), order=order, doPlot=False)

    def testRadial(self):
        """Add a radial transform"""
        for order in (4, 5, 6):
            self.doTest("testRadial", afwGeom.RadialXYTransform([0, 1.001, 0.000003]), order=order,
                        doPlot=False)

    def testWarnings(self):
        """Test that approximateWcs raises a UserWarning when it cannot achieve desired tolerance"""
        radialTransform = afwGeom.RadialXYTransform([0, 2.0, 3.0])
        wcs = afwImage.DistortedTanWcs(self.tanWcs, radialTransform)
        with self.assertRaises(UserWarning):
            approximateWcs(wcs=wcs, bbox=self.bbox, order=2)

    def doTest(self, name, xyTransform, order=3, doPlot=False):
        """Create a DistortedTanWcs from the specified transform and fit it
        """
        wcs = afwImage.DistortedTanWcs(self.tanWcs, xyTransform)

        fitWcs = approximateWcs(
            wcs=wcs,
            bbox=self.bbox,
            order=order,
        )

        if doPlot:
            self.plotWcs(wcs, fitWcs, self.bbox, xyTransform)

        msg = "ERROR: %s failed with order %s" % (name, order)
        self.assertWcsNearlyEqualOverBBox(wcs, fitWcs, self.bbox,
                                          maxDiffSky=0.001*afwGeom.arcseconds, maxDiffPix=0.02, msg=msg)

    def plotWcs(self, wcs0, wcs1, bbox, xyTransform):
        import matplotlib.pyplot as plt
        bboxd = afwGeom.Box2D(bbox)
        x0Arr = []
        y0Arr = []
        x1Arr = []
        y1Arr = []
        x2Arr = []
        y2Arr = []
        for x in np.linspace(bboxd.getMinX(), bboxd.getMaxX(), 10):
            for y in np.linspace(bboxd.getMinY(), bboxd.getMaxY(), 10):
                pixelPos0 = afwGeom.Point2D(x, y)
                skyCoord = wcs0.pixelToSky(pixelPos0)
                pixelPos1 = wcs1.skyToPixel(skyCoord)
                distortedPos = xyTransform.forwardTransform(pixelPos0)
                x0Arr.append(pixelPos0[0])
                y0Arr.append(pixelPos0[1])
                x1Arr.append(pixelPos1[0])
                y1Arr.append(pixelPos1[1])
                x2Arr.append(distortedPos[0])
                y2Arr.append(distortedPos[1])
        plt.plot(x0Arr, y0Arr, 'b+', x1Arr, y1Arr, 'rx', x2Arr, y2Arr, 'g.')

        plt.show()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
