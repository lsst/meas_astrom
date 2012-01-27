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
import lsst.afw.detection as afwDet
import lsst.meas.astrom.sip as astromSip

class CreateWcsWithSipTestCase(unittest.TestCase):
    """A test case for Calib"""

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

        self.sourceMatchSet = []
        for i in range(0, num, step):
            for j in range(0, num, step):
                src = afwDet.Source()
                cat = afwDet.Source()

                cat.setXAstrom(i)
                cat.setYAstrom(j)

                src.setXAstrom(i)
                src.setYAstrom(j)

                c = self.wcs.pixelToSky(i, j);
                cat.setRaDec(c);
                cat.setRaDecAstrom(c);
                if False:
                    print "RA,Dec = (%.3f, %.3f) deg" %  (c.toFk5().getRa().asDegrees(),
                                                          c.toFk5().getDec().asDegrees())
                dist = math.hypot(src.getXAstrom() - cat.getXAstrom(), src.getYAstrom() - cat.getYAstrom())
                self.sourceMatchSet.append(afwDet.SourceMatch(cat, src, dist))

    def tearDown(self):
        del self.sourceMatchSet
        del self.wcs
        del self.sipObject
        
    def assertAlmostEqualAngle(self, a1, a2):
        self.assertAlmostEqual(a1.asDegrees(), a2.asDegrees())

    def checkResults(self, sipWcs):
        for cat, src, d in self.sourceMatchSet:
            srcX = src.getXAstrom()
            srcY = src.getYAstrom()
            catRa  = cat.getRa();
            catDec = cat.getDec();
            srcRaDec = sipWcs.pixelToSky(srcX, srcY).toFk5();

            if False:
                print "cat RA,Dec = (%.5f, %.5f) deg" % (catRa.asDegrees(), catDec.asDegrees())
                print "src RA,Dec = (%.5f, %.5f) deg" % (srcRaDec.getRa().asDegrees(), srcRaDec.getDec().asDegrees())
            self.assertAlmostEqualAngle(catRa, srcRaDec.getRa())
            self.assertAlmostEqualAngle(catDec, srcRaDec.getDec())
            # these are in pixels.
            catxy = sipWcs.skyToPixel(cat.getRaDec());
            if False:
                print "cat X,Y = (%.3f, %.3f)" % (catxy[0], catxy[1])
                print "src X,Y = (%.3f, %.3f)" % (srcX, srcY);
            self.assertAlmostEqual(srcX, catxy[0])
            self.assertAlmostEqual(srcY, catxy[1])

    def doTest(self, name, func, order=2):
        """Apply func(x, y) to each point"""

        for cat, src, d in self.sourceMatchSet:
            x, y = func(src.getXAstrom(), src.getYAstrom())
            src.setXAstrom(x); src.setYAstrom(y)

        self.sipObject = astromSip.CreateWcsWithSip(self.sourceMatchSet, self.wcs, order)

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

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(CreateWcsWithSipTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
        
