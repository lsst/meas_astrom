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
import unittest

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.geom as afwGeom
from lsst.meas.astrom import approximateWcs


class ApproximateWcsTestCase(lsst.utils.tests.TestCase):

    """A test case for CreateWcsWithSip

    Use involves setting one class attribute:
    * MatchClass: match class, e.g. ReferenceMatch or SourceMatch
    """

    def setUp(self):
        self.crPix = lsst.geom.Point2D(15000, 4000)
        dimd = lsst.geom.Extent2D(4000, 4000)
        bboxd = lsst.geom.Box2D(self.crPix - dimd/2, dimd, invert=False)
        self.bbox = lsst.geom.Box2I(bboxd)
        self.tanWcs = afwGeom.makeSkyWcs(crpix=self.crPix,
                                         crval=lsst.geom.SpherePoint(215.5, 53.0, lsst.geom.degrees),
                                         cdMatrix=np.array([[5.10808596133527E-05, 1.85579539217196E-07],
                                                            [-8.27440751733828E-07, -5.10281493481982E-05]]))

    def tearDown(self):
        del self.tanWcs

    def testTrivial(self):
        """Add no distortion"""
        for order in (3, 4, 5, 6):
            self.doTest("testTrivial", afwGeom.makeIdentityTransform(), order=order, doPlot=False)

    def testRadial(self):
        """Add a radial transform"""
        for order in (4, 5, 6):
            self.doTest("testRadial", afwGeom.makeRadialTransform([0, 1.001, 0.000003]), order=order,
                        doPlot=False)

    def testWarnings(self):
        """Test that approximateWcs raises a UserWarning when it cannot achieve desired tolerance"""
        radialTransform = afwGeom.makeRadialTransform([0, 2.0, 3.0])
        wcs = afwGeom.makeModifiedWcs(pixelTransform=radialTransform, wcs=self.tanWcs,
                                      modifyActualPixels=False)
        with self.assertRaises(UserWarning):
            approximateWcs(wcs=wcs, bbox=self.bbox, order=2)

    def doTest(self, name, transform, order=3, doPlot=False):
        """Add the specified distorting transform to a TAN WCS and fit it

        The resulting WCS pixelToSky method acts as follows:
            pixelToSky(transform.applyForward(pixels))
        """
        wcs = afwGeom.makeModifiedWcs(pixelTransform=transform,
                                      wcs=self.tanWcs,
                                      modifyActualPixels=False)

        fitWcs = approximateWcs(
            wcs=wcs,
            bbox=self.bbox,
            order=order,
        )

        if doPlot:
            self.plotWcs(wcs, fitWcs, self.bbox, transform)

        msg = "ERROR: %s failed with order %s" % (name, order)
        self.assertWcsAlmostEqualOverBBox(wcs, fitWcs, self.bbox,
                                          maxDiffSky=0.001*lsst.geom.arcseconds, maxDiffPix=0.02, msg=msg)

    def plotWcs(self, wcs0, wcs1, bbox, transform):
        import matplotlib.pyplot as plt
        bboxd = lsst.geom.Box2D(bbox)
        x0Arr = []
        y0Arr = []
        x1Arr = []
        y1Arr = []
        x2Arr = []
        y2Arr = []
        for x in np.linspace(bboxd.getMinX(), bboxd.getMaxX(), 10):
            for y in np.linspace(bboxd.getMinY(), bboxd.getMaxY(), 10):
                pixelPos0 = lsst.geom.Point2D(x, y)
                skyCoord = wcs0.pixelToSky(pixelPos0)
                pixelPos1 = wcs1.skyToPixel(skyCoord)
                distortedPos = transform.applyForward(pixelPos0)
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
