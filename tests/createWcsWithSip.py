#!/usr/bin/env python

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

from lsst.afw.table import SimpleCatalog, SourceCatalog
import lsst.utils.tests
import lsst.afw.geom as afwGeom
from lsst.pex.logging import Log
from lsst.meas.astrom import ANetBasicAstrometryTask
import lsst.meas.astrom.sip as sip
import lsst.meas.astrom.sip.genDistortedImage as distort
import lsst.meas.astrom.sip.cleanBadPoints as cleanBadPoints
import testFindAstrometryNetDataDir as helper


class CreateWcsWithSipCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up local astrometry_net_data
        helper.setupAstrometryNetDataDir('cfhttemplate')

    def setUp(self):

        self.config = ANetBasicAstrometryTask.ConfigClass()
        self.config.defaultFilter = "r"
        self.astrom = ANetBasicAstrometryTask(config=self.config)

        testDir = os.path.dirname(__file__)
        self.filename = os.path.join(testDir, "cat.xy.fits")
        self.tolArcsec = .4
        self.tolPixel = .1

    def tearDown(self):
        del self.config
        del self.astrom

    def testBigXy0(self):
        # test for ticket #2710
        log = Log.getDefaultLog()
        log.setThreshold(Log.DEBUG)
        self.astrom.log = log
        x0, y0 = 200000, 500000
        cx = 500
        a2 = 1e-5
        cat = SourceCatalog.readFits(self.filename)
        print 'Catalog size', len(cat)
        # Source x,y positions are ~ (500,1500) x (500,1500)
        xKey = cat.table.getCentroidKey().getX()
        yKey = cat.table.getCentroidKey().getY()
        for src in cat:
            x = src.get(xKey) - 500
            dx = x - cx
            x += a2 * (dx**2)
            src.set(xKey, x + x0)
            src.set(yKey, src.get(yKey) - 500 + y0)
        bbox = afwGeom.Box2I(afwGeom.Point2I(x0, y0), afwGeom.Extent2I(1000, 1000))
        res = self.astrom.determineWcs2(cat, bbox=bbox)
        self.assertIsNotNone(res.sipWcs, "Failed to fit SIP terms")
        print 'Got result', res
        print 'SIP:', res.sipWcs.getFitsMetadata().toString()

        wcs = res.wcs
        for src in cat:
            rd = wcs.pixelToSky(src.getCentroid())
            xy = wcs.skyToPixel(rd)
            #print 'src', src.getX(), src.getY()
            #print 'rd', rd
            #print 'xy', xy
            #print 'dx,dy', xy[0] - src.getX(), xy[1] - src.getY()
            self.assertLess(abs(xy[0] - src.getX()), 0.1)
            self.assertLess(abs(xy[1] - src.getY()), 0.1)

    def testLinearXDistort(self):
        print "linearXDistort"
        self.singleTestInstance(self.filename, distort.linearXDistort)

    def testLinearYDistort(self):
        print "linearYDistort"
        self.singleTestInstance(self.filename, distort.linearYDistort)

    def testQuadraticDistort(self):
        print "linearQuadraticDistort"
        self.singleTestInstance(self.filename, distort.linearYDistort)

    def singleTestInstance(self, filename, distortFunc):
        cat = self.loadCatalogue(self.filename)
        img = distort.distortList(cat, distortFunc)
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(1000, 1000))
        res = self.astrom.determineWcs2(img, bbox=bbox)
        imgWcs = res.getWcs()

        def printWcs(wcs):
            print "WCS metadata:"
            md = wcs.getFitsMetadata()
            for name in md.names():
                print "%s: %r" % (name, md.get(name))

        printWcs(imgWcs)

        #Create a wcs with sip
        cat = cat.cast(SimpleCatalog, False)
        matchList = self.matchSrcAndCatalogue(cat, img, imgWcs)
        print "*** num matches =", len(matchList)
        return
        sipObject = sip.makeCreateWcsWithSip(matchList, imgWcs, 3)

        #print 'Put in TAN Wcs:'
        #print imgWcs.getFitsMetadata().toString()
        imgWcs = sipObject.getNewWcs()
        #print 'Got SIP Wcs:'
        #print imgWcs.getFitsMetadata().toString()

        # Write out the SIP header
        #afwImage.fits_write_imageF('createWcsWithSip.sip', afwImage.ImageF(0,0),
        #imgWcs.getFitsMetadata())

        print 'number of matches:', len(matchList), sipObject.getNPoints()
        scatter = sipObject.getScatterOnSky().asArcseconds()
        print "Scatter in arcsec is %g" % (scatter)

        self.assertLess(scatter, self.tolArcsec, "Scatter exceeds tolerance in arcsec")

        if False:
            scatter = sipObject.getScatterInPixels()
            self.assertLess(scatter, self.tolPixel, "Scatter exceeds tolerance in pixels: %g" % (scatter,))

    def loadCatalogue(self, filename):
        """Load a list of xy points from a file, solve for position, and
        return a SourceSet of points"""

        cat = SourceCatalog.readFits(filename)

        # Source x,y positions are ~ (500,1500) x (500,1500)
        xKey = cat.table.getCentroidKey().getX()
        yKey = cat.table.getCentroidKey().getY()
        for src in cat:
            src.set(xKey, src.get(xKey) - 500)
            src.set(yKey, src.get(yKey) - 500)

        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(1000, 1000))
        res = self.astrom.determineWcs2(cat, bbox=bbox)
        catWcs = res.getWcs()

        #Set catalogue ra and decs
        for src in cat:
            src.updateCoord(catWcs)
        return cat

    def matchSrcAndCatalogue(self, cat, img, imgWcs, dist=1.*afwGeom.arcseconds, cleanParam=3):
        """Given an input catalogue, match a list of objects in an image, given
        their x,y position and a wcs solution.

        Return: A list of x, y, dx and dy. Each element of the list is itself a list
        """

        matcher = sip.MatchSrcToCatalogue(cat, img, imgWcs, dist)
        matchList = matcher.getMatches()

        mList = cleanBadPoints.clean(matchList, imgWcs, order=cleanParam)
        return mList


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
