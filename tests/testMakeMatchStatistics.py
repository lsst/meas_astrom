#
# LSST Data Management System
# Copyright 2015 LSST Corporation.
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
from __future__ import absolute_import, division, print_function
from builtins import zip
import unittest
import math

import numpy as np

import lsst.utils.tests
from lsst.daf.base import PropertySet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
from lsst.meas.base import SingleFrameMeasurementTask
import lsst.meas.astrom as measAstrom
from functools import reduce


class TestAstrometricSolver(lsst.utils.tests.TestCase):

    def setUp(self):
        # make a nominal match list where the distances are 0; test can then modify
        # source centroid, reference coord or distance field for each match, as desired
        ctrPix = afwGeom.Point2I(1500, 1500)
        metadata = PropertySet()
        metadata.set("RADECSYS", "FK5")
        metadata.set("EQUINOX", 2000.0)
        metadata.set("CTYPE1", "RA---TAN")
        metadata.set("CTYPE2", "DEC--TAN")
        metadata.set("CUNIT1", "deg")
        metadata.set("CUNIT2", "deg")
        metadata.set("CRVAL1", 215.5)
        metadata.set("CRVAL2", 53.0)
        metadata.set("CRPIX1", ctrPix[0] + 1)
        metadata.set("CRPIX2", ctrPix[1] + 1)
        metadata.set("CD1_1", 5.1e-05)
        metadata.set("CD1_2", 0.0)
        metadata.set("CD2_2", -5.1e-05)
        metadata.set("CD2_1", 0.0)
        self.wcs = afwImage.makeWcs(metadata)
        self.bboxD = afwGeom.Box2D(afwGeom.Point2D(10, 100), afwGeom.Extent2D(1000, 1500))
        self.numMatches = 25

        sourceSchema = afwTable.SourceTable.makeMinimalSchema()
        # add centroid (and many other unwanted fields) to sourceSchema
        SingleFrameMeasurementTask(schema=sourceSchema)
        self.sourceCentroidKey = afwTable.Point2DKey(sourceSchema["slot_Centroid"])
        self.sourceCat = afwTable.SourceCatalog(sourceSchema)

        refSchema = afwTable.SourceTable.makeMinimalSchema()
        self.refCoordKey = afwTable.CoordKey(refSchema["coord"])
        self.refCat = afwTable.SourceCatalog(refSchema)

        self.matchList = []

        np.random.seed(5)
        pixPointList = [afwGeom.Point2D(pos) for pos in
                        np.random.random_sample([self.numMatches, 2])*self.bboxD.getDimensions() +
                        self.bboxD.getMin()]
        for pixPoint in pixPointList:
            src = self.sourceCat.addNew()
            src.set(self.sourceCentroidKey, pixPoint)
            ref = self.refCat.addNew()
            ref.set(self.refCoordKey, self.wcs.pixelToSky(pixPoint))

            match = afwTable.ReferenceMatch(ref, src, 0)
            self.matchList.append(match)

    def tearDown(self):
        del self.wcs
        del self.bboxD
        del self.sourceCentroidKey
        del self.sourceCat
        del self.refCoordKey
        del self.refCat
        del self.matchList

    def testMakeMatchStatistics(self):
        """Test makeMatchStatistics
        """
        np.random.seed(47)
        distList = list((np.random.random_sample([self.numMatches]) - 0.5) * 10)
        for dist, match in zip(distList, self.matchList):
            match.distance = dist
        itemList = (afwMath.MEDIAN, afwMath.MEANCLIP, afwMath.STDEV)
        itemMask = reduce(lambda a, b: a | b, itemList)
        distStats = measAstrom.makeMatchStatistics(self.matchList, itemMask)
        directStats = afwMath.makeStatistics(distList, itemMask)
        for item in itemList:
            self.assertAlmostEqual(distStats.getValue(item), directStats.getValue(item))

    def testMakeMatchStatisticsInRadians(self):
        """Test makeMatchStatisticsInRadians
        """
        np.random.seed(164)
        offLenList = [val*afwGeom.radians for val in np.random.random_sample([self.numMatches])]
        offDirList = [val*afwGeom.radians for val in np.random.random_sample([self.numMatches])*math.pi*2]
        for offLen, offDir, match in zip(offLenList, offDirList, self.matchList):
            coord = match.first.get(self.refCoordKey)
            coord.offset(offDir, offLen)  # an in-place operation
            match.first.set(self.refCoordKey, coord)
        itemList = (afwMath.MEDIAN, afwMath.MEANCLIP, afwMath.IQRANGE)
        itemMask = reduce(lambda a, b: a | b, itemList)
        distStats = measAstrom.makeMatchStatisticsInRadians(self.wcs, self.matchList, itemMask)
        distRadiansList = [val.asRadians() for val in offLenList]
        directStats = afwMath.makeStatistics(distRadiansList, itemMask)
        for item in itemList:
            self.assertAlmostEqual(distStats.getValue(item), directStats.getValue(item))

    def testMakeMatchStatisticsInPixels(self):
        """Test testMakeMatchStatisticsInPixels
        """
        np.random.seed(164)
        offList = [afwGeom.Extent2D(val) for val in (np.random.random_sample([self.numMatches, 2])-0.5)*10]
        for off, match in zip(offList, self.matchList):
            centroid = match.second.get(self.sourceCentroidKey)
            offCentroid = centroid + off
            match.second.set(self.sourceCentroidKey, offCentroid)
        itemList = (afwMath.MEDIAN, afwMath.MEAN, afwMath.STDEVCLIP)
        itemMask = reduce(lambda a, b: a | b, itemList)
        distStats = measAstrom.makeMatchStatisticsInPixels(self.wcs, self.matchList, itemMask)
        distList = [math.hypot(*val) for val in offList]
        directStats = afwMath.makeStatistics(distList, itemMask)
        for item in itemList:
            self.assertAlmostEqual(distStats.getValue(item), directStats.getValue(item))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
