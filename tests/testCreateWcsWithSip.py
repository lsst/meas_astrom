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

import matplotlib
matplotlib.use('Agg')
import pylab as plt
pnum = 1

import math, os, sys
import unittest
import numpy as np
import lsst.utils.tests as tests

import lsst.pex.logging as pexLog
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.astrom.sip as astromSip

from lsst.pex.logging import Log
log = Log.getDefaultLog()
log.setThreshold(Log.DEBUG)


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

        S = 4000
        N = 10

        catTable = afwTable.SimpleTable.make(afwTable.SimpleTable.makeMinimalSchema())
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        key = srcSchema.addField("centroid", type="PointD")
        srcSchema.addField("centroid.flags", type="Flag")
        srcSchema.addField("centroid.err", type="CovPointF")
        srcTable = afwTable.SourceTable.make(srcSchema)
        srcTable.defineCentroid("centroid")
        self.sourceMatchSet = []
        for i in np.linspace(0., S, N):
            for j in np.linspace(0., S, N):
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



        self.bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(S,S))

    def tearDown(self):
        del self.sourceMatchSet
        del self.wcs
        del self.sipObject
        
    def assertAlmostEqualAngle(self, a1, a2):
        self.assertAlmostEqual(a1.asArcseconds(), a2.asArcseconds(), 1) # 100 mas tolerance

    def checkResults(self, sipWcs):
        for cat, src, d in self.sourceMatchSet:
            srcX = src.getX()
            srcY = src.getY()
            catRa  = cat.getRa();
            catDec = cat.getDec();
            srcCoord = sipWcs.pixelToSky(srcX, srcY).toFk5();

            if True:
                print "cat RA,Dec = (%.8f, %.8f) deg" % (catRa.asDegrees(), catDec.asDegrees())
                print "src RA,Dec = (%.8f, %.8f) deg" % (srcCoord.getRa().asDegrees(), srcCoord.getDec().asDegrees())
            self.assertAlmostEqualAngle(catRa, srcCoord.getRa())
            self.assertAlmostEqualAngle(catDec, srcCoord.getDec())
            # these are in pixels.
            catxy = sipWcs.skyToPixel(cat.getCoord());
            if True:
                print "cat X,Y = (%.3f, %.3f)" % (catxy[0], catxy[1])
                print "src X,Y = (%.3f, %.3f)" % (srcX, srcY);
            self.assertAlmostEqual(srcX, catxy[0], 3) # within a milli-pixel
            self.assertAlmostEqual(srcY, catxy[1], 3)

    def doTest(self, name, func, order=2):
        """Apply func(x, y) to each point"""

        print
        print name
        print

        keepmatches = []
        for cat, src, d in self.sourceMatchSet:
            x, y = func(src.getX(), src.getY())
            src.set("centroid.x", x); src.set("centroid.y", y)
            if self.bbox.contains(afwGeom.Point2I(int(np.round(x)),int(np.round(y)))):
                keepmatches.append(afwTable.ReferenceMatch(cat, src, d))
        self.sourceMatchSet = keepmatches

        self.sipObject = astromSip.makeCreateWcsWithSip(self.sourceMatchSet, self.wcs, order, self.bbox)

        if False:
            print name
            print self.sipObject.getNewWcs().getFitsMetadata().toString()

        wcs = self.sipObject.getNewWcs()
        for cat, src, d in self.sourceMatchSet:
            src.updateCoord(wcs)

        print 'Got WCS', self.sipObject.getNewWcs().getFitsMetadata().toString()

        xs,ys, xc,yc = [],[],[],[]
        rs,ds, rc,dc = [],[],[],[]
        for cat, src, d in self.sourceMatchSet:
            xs.append(src.getX())
            ys.append(src.getY())
            catxy = wcs.skyToPixel(cat.getCoord())
            xc.append(catxy[0])
            yc.append(catxy[1])
            rc.append(cat.getRa())
            dc.append(cat.getDec())
            srd = wcs.pixelToSky(src.getX(), src.getY()).toFk5()
            rs.append(srd.getRa())
            ds.append(srd.getDec())
        xs = np.array(xs)
        ys = np.array(ys)
        xc = np.array(xc)
        yc = np.array(yc)
            
        global pnum
        plt.clf()
        plt.plot(xs, ys, 'r.')
        plt.plot(xc, yc, 'bx')
        fn = 'check-%i.png' % pnum
        plt.savefig(fn)
        print 'Wrote', fn
        pnum += 1

        plt.clf()
        plt.plot(xs, xc-xs, 'b.')
        fn = 'check-%i.png' % pnum
        plt.xlabel('x(source)')
        plt.ylabel('x(cat - src)')
        plt.savefig(fn)
        print 'Wrote', fn
        pnum += 1

        plt.clf()
        plt.plot(rs, ds, 'r.')
        plt.plot(rc, dc, 'bx')
        fn = 'check-%i.png' % pnum
        plt.savefig(fn)
        print 'Wrote', fn
        pnum += 1

        plt.clf()
        for y in [0., 1000., 2000., 3000., 4000.]:
            x0,y0 = [],[]
            x1,y1 = [],[]
            for x in np.linspace(0., 4000., 400):
                x0.append(x)
                y0.append(y)
                rd = wcs.pixelToSky(x, y)
                xy = wcs.skyToPixel(rd)
                x1.append(xy[0])
                y1.append(xy[1])
            x0 = np.array(x0)
            x1 = np.array(x1)
            plt.plot(x0, x1-x0, 'b-')
        fn = 'check-%i.png' % pnum
        plt.savefig(fn)
        print 'Wrote', fn
        pnum += 1

        plt.clf()
        for y in [0., 1000., 2000., 3000., 4000.]:
            x0,y0 = [],[]
            x1,y1 = [],[]
            x2,y2 = [],[]
            for x in np.linspace(0., 4000., 400):
                x0.append(x)
                y0.append(y)
                xy = wcs.undistortPixel(afwGeom.Point2D(x+1,y+1))
                x1.append(xy[0]-1)
                y1.append(xy[1]-1)
                xy = wcs.distortPixel(xy)
                x2.append(xy[0]-1)
                y2.append(xy[1]-1)
            x0 = np.array(x0)
            x1 = np.array(x1)
            x2 = np.array(x2)
            plt.plot(x0, x1-x0, 'b-')
            plt.plot(x0, x1-x2, 'r-')
            plt.plot(x0, (x2-x0)*1e5, 'g-')
        plt.xlabel('x (orig)')
        plt.ylabel('dx (undistorted)')
        fn = 'check-%i.png' % pnum
        plt.savefig(fn)
        print 'Wrote', fn
        pnum += 1
        
        self.checkResults(self.sipObject.getNewWcs())

    def tstTrivial(self):
        """Add no distortion"""
        self.doTest(__doc__, lambda x, y: (x, y))

    def tstOffset(self):
        """Add an offset"""
        self.doTest(__doc__, lambda x, y: (x + 5, y + 7))

    def tstLinearX(self):
        """Scale x, offset y"""
        self.doTest(__doc__, lambda x, y: (2*x, y + 7))

    def tstLinearXY(self):
        """Scale x and y"""
        self.doTest(__doc__, lambda x, y: (2*x, 3*y))

    def tstLinearYX(self):
        """Add an offset to each point; scale in y and x"""
        self.doTest(__doc__, lambda x, y: (x + 0.2*y, y + 0.3*x))

    def testQuadraticX(self):
        """Add quadratic distortion in x"""
        order = 3
        self.doTest('quadraticX', lambda x, y: (x + 1e-5*x**2, y), order=order)
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
