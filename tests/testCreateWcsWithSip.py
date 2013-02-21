#!/usr/bin/env python

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
 
"""
Tests for TanSip

Run with:
   CreateWcsWithSip.py
or
   python
   >>> import CreateWcsWithSip; CreateWcsWithSip.run()
"""

import math, os, sys
import unittest
import lsst.utils.tests as tests

import lsst.pex.logging as pexLog
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.astrom.sip as astromSip

class BaseTestCase(unittest.TestCase):
    """A test case for CreateWcsWithSip

    Use involves setting a couple class attributes:
    * CatTableClass: Reference catalog Table class, e.g., SimpleTable
    * MatchClass: Match class (which should be compatible with the CatTableClass), e.g., ReferenceMatch
    """
    CatTableClass = None
    MatchClass = None

    def setUp(self):
        crval = afwCoord.Coord(afwGeom.PointD(44., 45.))
        crpix = afwGeom.PointD(0, 0)
        
        arcsecPerPixel = 1/3600.0
        CD11 = arcsecPerPixel
        CD12 = 0
        CD21 = 0
        CD22 = arcsecPerPixel
        
        self.wcs = afwImage.makeWcs(crval, crpix, CD11, CD12, CD21, CD22)
        self.sipObject = None

        num = 4000
        step = 1000

        catTable = afwTable.SimpleTable.make(afwTable.SimpleTable.makeMinimalSchema())
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        key = srcSchema.addField("centroid", type="PointD")
        srcSchema.addField("centroid.flags", type="Flag")
        srcSchema.addField("centroid.err", type="CovPointF")
        srcTable = afwTable.SourceTable.make(srcSchema)
        srcTable.defineCentroid("centroid")
        self.sourceMatchSet = []
        for i in range(0, num, step):
            for j in range(0, num, step):
                src = srcTable.makeRecord()
                cat = catTable.makeRecord()

                src.set(key.getX(), i)
                src.set(key.getY(), j)

                c = self.wcs.pixelToSky(i, j);
                cat.setCoord(c);

                if False:
                    print "RA,Dec = (%.3f, %.3f) deg" %  (c.toFk5().getRa().asDegrees(),
                                                          c.toFk5().getDec().asDegrees())

                self.sourceMatchSet.append(afwTable.ReferenceMatch(cat, src, 0.0))

    def tearDown(self):
        del self.sourceMatchSet
        del self.wcs
        del self.sipObject
        
    def assertAlmostEqualAngle(self, a1, a2):
        self.assertAlmostEqual(a1.asArcseconds(), a2.asArcseconds(), 3) # 1 mas tolerance

    def checkResults(self, sipWcs):
        for cat, src, d in self.sourceMatchSet:
            srcX = src.getX()
            srcY = src.getY()
            catRa  = cat.getRa();
            catDec = cat.getDec();
            srcCoord = sipWcs.pixelToSky(srcX, srcY).toFk5();

            if False:
                print "cat RA,Dec = (%.5f, %.5f) deg" % (catRa.asDegrees(), catDec.asDegrees())
                print "src RA,Dec = (%.5f, %.5f) deg" % (srcCoord.getRa().asDegrees(), srcCoord.getDec().asDegrees())
            self.assertAlmostEqualAngle(catRa, srcCoord.getRa())
            self.assertAlmostEqualAngle(catDec, srcCoord.getDec())
            # these are in pixels.
            catxy = sipWcs.skyToPixel(cat.getCoord());
            if False:
                print "cat X,Y = (%.3f, %.3f)" % (catxy[0], catxy[1])
                print "src X,Y = (%.3f, %.3f)" % (srcX, srcY);
            self.assertAlmostEqual(srcX, catxy[0], 3) # within a milli-pixel
            self.assertAlmostEqual(srcY, catxy[1], 3)

    def doTest(self, name, func, order=2):
        """Apply func(x, y) to each point"""

        for cat, src, d in self.sourceMatchSet:
            x, y = func(src.getX(), src.getY())
            src.set("centroid.x", x); src.set("centroid.y", y)

        self.sipObject = astromSip.makeCreateWcsWithSip(self.sourceMatchSet, self.wcs, order)

        if False:
            print name
            print self.sipObject.getNewWcs().getFitsMetadata().toString()

        self.checkResults(self.sipObject.getNewWcs())

    def testTrivial(self):
        """Add no distortion"""
        self.doTest(__doc__, lambda x, y: (x, y))

    def testOffset(self):
        """Add an offset"""
        self.doTest(__doc__, lambda x, y: (x + 5, y + 7))

    def testLinearX(self):
        """Scale x, offset y"""
        self.doTest(__doc__, lambda x, y: (2*x, y + 7))

    def testLinearXY(self):
        """Scale x and y"""
        self.doTest(__doc__, lambda x, y: (2*x, 3*y))

    def testLinearYX(self):
        """Add an offset to each point; scale in y and x"""
        self.doTest(__doc__, lambda x, y: (x + 0.2*y, y + 0.3*x))

    def testQuadraticX(self):
        """Add quadratic distortion in x"""
        order = 3
        self.doTest(__doc__, lambda x, y: (x + 1e-5*x**2, y), order=order)
        if False:
            afwImage.ImageF().writeFits("bad%d.fits" % order, self.sipObject.getNewWcs().getFitsMetadata())


def makeTestCase(_CatTableClass, _MatchClass):
    class CreateWcsWithSipTestCase(BaseTestCase):
        CatTableClass = _CatTableClass
        MatchClass = _MatchClass
    return CreateWcsWithSipTestCase

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
#    suites += unittest.makeSuite(CreateWcsWithSipTestCase)
    suites += unittest.makeSuite(makeTestCase(afwTable.SimpleTable, afwTable.ReferenceMatch))
    suites += unittest.makeSuite(makeTestCase(afwTable.SourceTable, afwTable.SourceMatch))
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
