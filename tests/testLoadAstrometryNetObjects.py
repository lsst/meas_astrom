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

import os
import unittest

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.astrom as measAstrom

import testFindAstrometryNetDataDir as helper


class TestLoadAstrometryNetObjects(unittest.TestCase):

    def setUp(self):
        # Set up local astrometry_net_data
        self.datapath = helper.setupAstrometryNetDataDir('photocal')
        self.config = measAstrom.LoadAstrometryNetObjectsTask.ConfigClass()

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
        self.wcs = afwImage.makeWcs(metadata)
        self.desNumStarsInPixelBox = 270
        self.desNumStarsInSkyCircle = 410

    def tearDown(self):
        del self.ctrPix
        del self.wcs
        del self.config

    def testLoadPixelBox(self):
        """Test loadPixelBox
        """
        loadANetObj = measAstrom.LoadAstrometryNetObjectsTask(config=self.config)

        loadRes = loadANetObj.loadPixelBox(bbox=self.bbox, wcs=self.wcs, filterName="r")
        refCat = loadRes.refCat
#        self.plotStars(refCat, bbox=self.bbox)
        self.assertEqual(loadRes.fluxField, "r_flux")
        self.assertEqual(len(refCat), self.desNumStarsInPixelBox)
        self.assertObjInBBox(refCat=refCat, bbox=self.bbox, wcs=self.wcs)
        schema = refCat.getSchema()
        filterNameList = ['u', 'g', 'r', 'i', 'z']
        for filterName in filterNameList:
            schema.find(filterName + "_flux").key
            schema.find(filterName + "_fluxSigma").key
        for fieldName in ("coord_ra", "coord_dec", "centroid_x", "centroid_y", "hasCentroid",
                          "photometric", "resolved"):
            schema.find(fieldName)

    def testLoadSkyCircle(self):
        loadANetObj = measAstrom.LoadAstrometryNetObjectsTask(config=self.config)

        ctrCoord = self.wcs.pixelToSky(afwGeom.Point2D(self.ctrPix))
        radius = ctrCoord.angularSeparation(self.wcs.pixelToSky(afwGeom.Box2D(self.bbox).getMin()))

        loadRes = loadANetObj.loadSkyCircle(ctrCoord=ctrCoord, radius=radius, filterName="r")
        self.assertEqual(len(loadRes.refCat), self.desNumStarsInSkyCircle)

    def testNoMagErrs(self):
        """Exclude magnitude errors from the found catalog
        """
        andConfig = measAstrom.AstrometryNetDataConfig()
        andConfig.load(os.path.join(self.datapath, 'andConfig2.py'))
        andConfig.magErrorColumnMap = {}
        loadANetObj = measAstrom.LoadAstrometryNetObjectsTask(config = self.config, andConfig = andConfig)

        loadRes = loadANetObj.loadPixelBox(bbox=self.bbox, wcs=self.wcs, filterName="r")
        refCat = loadRes.refCat
        self.assertEqual(loadRes.fluxField, "r_flux")
        self.assertEqual(len(refCat), self.desNumStarsInPixelBox)
        self.assertObjInBBox(refCat=refCat, bbox=self.bbox, wcs=self.wcs)
        schema = refCat.getSchema()
        for filterName in ['u', 'g', 'r', 'i', 'z']:
            schema.find(filterName + "_flux")
            self.assertRaises(KeyError, schema.find, filterName + "_fluxSigma")

    def testRequestForeignFilter(self):
        """The user requests a filter not in the astrometry.net catalog.

        In that case, we must specify a mapping in the AstrometryConfig to point
        to an alternative filter (e.g., g instead of B).
        We should expect the returned catalog to contain references
        to the filterNameList that are in the catalog.
        """
        filterNameList = ['u', 'g', 'r', 'i', 'z']
        andConfig = measAstrom.AstrometryNetDataConfig()
        andConfig.load(os.path.join(self.datapath, 'andConfig2.py'))
        self.config.filterMap = dict(('my_'+b, b) for b in filterNameList)
        loadANetObj = measAstrom.LoadAstrometryNetObjectsTask(config = self.config, andConfig = andConfig)

        loadRes = loadANetObj.loadPixelBox(bbox=self.bbox, wcs=self.wcs, filterName="my_r")
        refCat = loadRes.refCat
        self.assertEqual(loadRes.fluxField, "my_r_camFlux")
        self.assertEqual(len(refCat), self.desNumStarsInPixelBox)
        self.assertObjInBBox(refCat=refCat, bbox=self.bbox, wcs=self.wcs)
        schema = refCat.getSchema()
        for filterName in filterNameList:
            schema.find(filterName + "_flux")
            schema.find(filterName + '_fluxSigma')

    def testDifferentMagNames(self):
        """The astrometry.net catalog's magnitude columns are not named after filters.

        In that case, the AstrometryNetDataConfig has a mapping to point to the correct columns.
        We should expect that the returned catalog refers to the filter
        requested (not the implementation-dependent column names).
        """
        andConfig = measAstrom.AstrometryNetDataConfig()
        andConfig.load(os.path.join(self.datapath, 'andConfig2.py'))
        baseNameList = ('u', 'g', 'r', 'i', 'z')
        filterNameList = ["my_" + b for b in baseNameList]
        andConfig.magColumnMap = dict(("my_" + b, b) for b in baseNameList)
        andConfig.magErrorColumnMap = dict([('my_' + b, b + "_err") for b in baseNameList])
        loadANetObj = measAstrom.LoadAstrometryNetObjectsTask(config = self.config, andConfig = andConfig)

        loadRes = loadANetObj.loadPixelBox(bbox=self.bbox, wcs=self.wcs, filterName="my_r")
        refCat = loadRes.refCat
        self.assertEqual(loadRes.fluxField, "my_r_flux")
        self.assertEqual(len(refCat), self.desNumStarsInPixelBox)
        self.assertObjInBBox(refCat=refCat, bbox=self.bbox, wcs=self.wcs)
        schema = refCat.getSchema()
        for nm in filterNameList:
            schema.find(nm + "_flux")
            schema.find(nm + '_fluxSigma')

    def assertObjInBBox(self, refCat, bbox, wcs):
        """Assert that all reference objects are inside the specified pixel bounding box plus a margin

        @param[in] refCat  reference object catalog, an lsst.afw.table.SimpleCatalog or compatible;
            the only fields read are "centroid_x/y" and "coord_ra/dec"
        @param[in] bbox  pixel bounding box coordinates, an lsst.afw.geom.Box2I or Box2D;
            the supplied box is grown by self.config.pixelMargin before testing the stars
        @param[in] wcs  WCS, an lsst.afw.image.Wcs
        """
        bbox = afwGeom.Box2D(bbox)
        bbox.grow(self.config.pixelMargin)
        centroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
        coordKey = afwTable.CoordKey(refCat.schema["coord"])
        for refObj in refCat:
            point = refObj.get(centroidKey)
            if not bbox.contains(point):
                coord = refObj.get(coordKey)
                self.fail("refObj at RA, Dec %0.3f, %0.3f point %s is not in bbox %s" % \
                    (coord[0].asDegrees(), coord[1].asDegrees(), point, bbox))

    def plotStars(self, refCat, bbox=None):
        """Plot the centroids of reference objects, and the bounding box (if specified)
        """
        import matplotlib.pyplot as pyplot
        if bbox is not None:
            cornerList = list(afwGeom.Box2D(bbox).getCorners())
            cornerList.append(cornerList[0]) # show 4 sides of the box by going back to the beginning
            xc, yc = zip(*cornerList)
            pyplot.plot(xc, yc, '-')

        centroidKey = refCat.schema["centroid"].asKey()
        centroidList = [rec.get(centroidKey) for rec in refCat]
        xp, yp = zip(*centroidList)
        pyplot.plot(xp, yp, '.')
        pyplot.show()


            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(TestLoadAstrometryNetObjects)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
