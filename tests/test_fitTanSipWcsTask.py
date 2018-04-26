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
# The classes in this test are a little non-standard to reduce code
# duplication and support automated unittest discovery.
# A base class includes all the code that implements the testing and
# itself inherits from unittest.TestCase. unittest automated discovery
# will scan all classes that inherit from unittest.TestCase and invoke
# any test methods found. To prevent this base class from being executed
# the test methods are placed in a different class that does not inherit
# from unittest.TestCase. The actual test classes then inherit from
# both the testing class and the implementation class allowing test
# discovery to only run tests found in the subclasses.

import math
import unittest

import numpy as np

import lsst.pipe.base
import lsst.utils.tests
import lsst.afw.geom as afwGeom
from lsst.afw.geom.wcsUtils import makeTanSipMetadata
import lsst.afw.table as afwTable
from lsst.meas.algorithms import LoadReferenceObjectsTask
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.astrom import FitTanSipWcsTask, setMatchDistance
from lsst.meas.astrom.sip import makeCreateWcsWithSip


class BaseTestCase(object):

    """A test case for CreateWcsWithSip

    Use involves setting one class attribute:
    * MatchClass: match class, e.g. ReferenceMatch or SourceMatch

    This test is a bit messy because it exercises two templatings of makeCreateWcsWithSip,
    the underlying TAN-SIP WCS fitter, but only one of those is supported by FitTanSipWcsTask
    """
    MatchClass = None

    def setUp(self):
        crval = afwGeom.SpherePoint(44, 45, afwGeom.degrees)
        crpix = afwGeom.Point2D(15000, 4000)

        scale = 1 * afwGeom.arcseconds
        cdMatrix = afwGeom.makeCdMatrix(scale=scale, flipX=True)
        self.tanWcs = afwGeom.makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cdMatrix)
        self.loadData()

    def loadData(self, rangePix=3000, numPoints=25):
        """Load catalogs and make the match list

        This is a separate function so data can be reloaded if fitting more than once
        (each time a WCS is fit it may update the source catalog, reference catalog and match list)
        """
        if self.MatchClass == afwTable.ReferenceMatch:
            refSchema = LoadReferenceObjectsTask.makeMinimalSchema(
                filterNameList=["r"], addFluxSigma=True, addIsPhotometric=True)
            self.refCat = afwTable.SimpleCatalog(refSchema)
        elif self.MatchClass == afwTable.SourceMatch:
            refSchema = afwTable.SourceTable.makeMinimalSchema()
            self.refCat = afwTable.SourceCatalog(refSchema)
        else:
            raise RuntimeError("Unsupported MatchClass=%r" % (self.MatchClass,))
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        SingleFrameMeasurementTask(schema=srcSchema)
        self.srcCoordKey = afwTable.CoordKey(srcSchema["coord"])
        self.srcCentroidKey = afwTable.Point2DKey(srcSchema["slot_Centroid"])
        self.srcCentroidKey_xSigma = srcSchema["slot_Centroid_xSigma"].asKey()
        self.srcCentroidKey_ySigma = srcSchema["slot_Centroid_ySigma"].asKey()
        self.sourceCat = afwTable.SourceCatalog(srcSchema)

        self.matches = []

        for i in np.linspace(0., rangePix, numPoints):
            for j in np.linspace(0., rangePix, numPoints):
                src = self.sourceCat.addNew()
                refObj = self.refCat.addNew()

                src.set(self.srcCentroidKey, afwGeom.Point2D(i, j))
                src.set(self.srcCentroidKey_xSigma, 0.1)
                src.set(self.srcCentroidKey_ySigma, 0.1)

                c = self.tanWcs.pixelToSky(i, j)
                refObj.setCoord(c)

                if False:
                    print("x,y = (%.1f, %.1f) pixels -- RA,Dec = (%.3f, %.3f) deg" %
                          (i, j, c.toFk5().getRa().asDegrees(), c.toFk5().getDec().asDegrees()))

                self.matches.append(self.MatchClass(refObj, src, 0.0))

    def tearDown(self):
        del self.refCat
        del self.sourceCat
        del self.matches
        del self.tanWcs

    def checkResults(self, fitRes, catsUpdated):
        """Check results

        @param[in] fitRes  a object with two fields:
            - wcs  fit TAN-SIP WCS, an lsst.afw.geom.SkyWcs
            - scatterOnSky  median on-sky scatter, an lsst.afw.geom.Angle
        @param[in] catsUpdated  if True then coord field of self.sourceCat and centroid fields of self.refCat
            have been updated
        """
        self.assertLess(fitRes.scatterOnSky.asArcseconds(), 0.001)
        tanSipWcs = fitRes.wcs
        maxAngSep = afwGeom.Angle(0)
        maxPixSep = 0
        refCoordKey = afwTable.CoordKey(self.refCat.schema["coord"])
        if catsUpdated:
            refCentroidKey = afwTable.Point2DKey(self.refCat.schema["centroid"])
        maxDistErr = afwGeom.Angle(0)
        for refObj, src, distRad in self.matches:
            srcPixPos = src.get(self.srcCentroidKey)
            refCoord = refObj.get(refCoordKey)
            if catsUpdated:
                refPixPos = refObj.get(refCentroidKey)
                srcCoord = src.get(self.srcCoordKey)
            else:
                refPixPos = tanSipWcs.skyToPixel(refCoord)
                srcCoord = tanSipWcs.pixelToSky(srcPixPos)

            angSep = refCoord.separation(srcCoord)
            dist = distRad*afwGeom.radians
            distErr = abs(dist - angSep)
            maxDistErr = max(maxDistErr, distErr)
            maxAngSep = max(maxAngSep, angSep)

            pixSep = math.hypot(*(srcPixPos - refPixPos))
            maxPixSep = max(maxPixSep, pixSep)

        print("max angular separation = %0.4f arcsec" % (maxAngSep.asArcseconds(),))
        print("max pixel separation = %0.3f" % (maxPixSep,))
        self.assertLess(maxAngSep.asArcseconds(), 0.001)
        self.assertLess(maxPixSep, 0.005)
        if catsUpdated:
            allowedDistErr = 1e-7
        else:
            allowedDistErr = 0.001
        self.assertLess(maxDistErr.asArcseconds(), allowedDistErr,
                        "Computed distance in match list is off by %s arcsec" % (maxDistErr.asArcseconds(),))

    def doTest(self, name, func, order=3, numIter=4, specifyBBox=False, doPlot=False, doPrint=False):
        """Apply func(x, y) to each source in self.sourceCat, then fit and check the resulting WCS
        """
        bbox = afwGeom.Box2I()
        for refObj, src, d in self.matches:
            origPos = src.get(self.srcCentroidKey)
            x, y = func(*origPos)
            distortedPos = afwGeom.Point2D(*func(*origPos))
            src.set(self.srcCentroidKey, distortedPos)
            bbox.include(afwGeom.Point2I(afwGeom.Point2I(distortedPos)))

        tanSipWcs = self.tanWcs
        for i in range(numIter):
            if specifyBBox:
                sipObject = makeCreateWcsWithSip(self.matches, tanSipWcs, order, bbox)
            else:
                sipObject = makeCreateWcsWithSip(self.matches, tanSipWcs, order)
            tanSipWcs = sipObject.getNewWcs()
        setMatchDistance(self.matches)
        fitRes = lsst.pipe.base.Struct(
            wcs=tanSipWcs,
            scatterOnSky=sipObject.getScatterOnSky(),
        )

        if doPrint:
            print("TAN-SIP metadata fit over bbox=", bbox)
            metadata = makeTanSipMetadata(
                crpix=tanSipWcs.getPixelOrigin(),
                crval=tanSipWcs.getSkyOrigin(),
                cdMatrix=tanSipWcs.getCdMatrix(),
                sipA=sipObject.getSipA(),
                sipB=sipObject.getSipB(),
                sipAp=sipObject.getSipAp(),
                sipBp=sipObject.getSipBp(),
            )
            print(metadata.toString())

        if doPlot:
            self.plotWcs(tanSipWcs, name=name)

        self.checkResults(fitRes, catsUpdated=False)

        if self.MatchClass == afwTable.ReferenceMatch:
            # reset source coord and reference centroid based on initial WCS
            afwTable.updateRefCentroids(wcs=self.tanWcs, refList=self.refCat)
            afwTable.updateSourceCoords(wcs=self.tanWcs, sourceList=self.sourceCat)

            fitterConfig = FitTanSipWcsTask.ConfigClass()
            fitterConfig.order = order
            fitterConfig.numIter = numIter
            fitter = FitTanSipWcsTask(config=fitterConfig)
            self.loadData()
            if specifyBBox:
                fitRes = fitter.fitWcs(
                    matches=self.matches,
                    initWcs=self.tanWcs,
                    bbox=bbox,
                    refCat=self.refCat,
                    sourceCat=self.sourceCat,
                )
            else:
                fitRes = fitter.fitWcs(
                    matches=self.matches,
                    initWcs=self.tanWcs,
                    bbox=bbox,
                    refCat=self.refCat,
                    sourceCat=self.sourceCat,
                )

            self.checkResults(fitRes, catsUpdated=True)

    def plotWcs(self, tanSipWcs, name=""):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fileNamePrefix = "testCreateWcsWithSip_%s_%s" % (self.MatchClass.__name__, name)
        pnum = 1

        xs, ys, xc, yc = [], [], [], []
        rs, ds, rc, dc = [], [], [], []
        for ref, src, d in self.matches:
            xs.append(src.getX())
            ys.append(src.getY())
            refPixPos = tanSipWcs.skyToPixel(ref.getCoord())
            xc.append(refPixPos[0])
            yc.append(refPixPos[1])
            rc.append(ref.getRa())
            dc.append(ref.getDec())
            srd = tanSipWcs.pixelToSky(src.get(self.srcCentroidKey))
            rs.append(srd.getRa())
            ds.append(srd.getDec())
        xs = np.array(xs)
        ys = np.array(ys)
        xc = np.array(xc)
        yc = np.array(yc)

        plt.clf()
        plt.plot(xs, ys, "r.")
        plt.plot(xc, yc, "bx")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        plt.savefig(fileName)
        print("Wrote", fileName)
        pnum += 1

        plt.clf()
        plt.plot(xs, xc-xs, "b.")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        plt.xlabel("x(source)")
        plt.ylabel("x(ref - src)")
        plt.savefig(fileName)
        print("Wrote", fileName)
        pnum += 1

        plt.clf()
        plt.plot(rs, ds, "r.")
        plt.plot(rc, dc, "bx")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        plt.savefig(fileName)
        print("Wrote", fileName)
        pnum += 1

        plt.clf()
        for y in np.linspace(0, 4000, 5):
            x0, y0 = [], []
            x1, y1 = [], []
            for x in np.linspace(0., 4000., 401):
                x0.append(x)
                y0.append(y)
                rd = tanSipWcs.pixelToSky(x, y)
                xy = tanSipWcs.skyToPixel(rd)
                x1.append(xy[0])
                y1.append(xy[1])
            x0 = np.array(x0)
            x1 = np.array(x1)
            plt.plot(x0, x1-x0, "b-")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        plt.savefig(fileName)
        print("Wrote", fileName)
        pnum += 1


class SideLoadTestCases(object):

    """Base class implementations of testing methods.

    Explicitly does not inherit from unittest.TestCase"""

    def testTrivial(self):
        """Add no distortion"""
        for order in (2, 4, 6):
            self.doTest("testTrivial", lambda x, y: (x, y), order=order)

    def testOffset(self):
        """Add an offset"""
        for order in (2, 4, 6):
            self.doTest("testOffset", lambda x, y: (x + 5, y + 7), order=order)

    def testLinearX(self):
        """Scale x, offset y"""
        for order in (2, 6):
            self.doTest("testLinearX", lambda x, y: (2*x, y + 7), order=order)

    def testLinearXY(self):
        """Scale x and y"""
        self.doTest("testLinearXY", lambda x, y: (2*x, 3*y))

    def testLinearYX(self):
        """Add an offset to each point; scale in y and x"""
        for order in (2, 6):
            self.doTest("testLinearYX", lambda x, y: (x + 0.2*y, y + 0.3*x), order=order)

    def testQuadraticX(self):
        """Add quadratic distortion in x"""
        for order in (4, 5):
            self.doTest("testQuadraticX", lambda x, y: (x + 1e-5*x**2, y), order=order)

    def testRadial(self):
        """Add radial distortion"""
        radialTransform = afwGeom.makeRadialTransform([0, 1.01, 1e-8])

        def radialDistortion(x, y):
            x, y = radialTransform.applyForward(afwGeom.Point2D(x, y))
            return (x, y)
        for order in (4, 5, 6):
            doPrint = order == 5
            self.doTest("testRadial", radialDistortion, order=order, doPrint=doPrint)

# The test classes inherit from two base classes and differ in the match
# class being used.


class CreateWcsWithSipTestCaseReferenceMatch(BaseTestCase, SideLoadTestCases, lsst.utils.tests.TestCase):
    MatchClass = afwTable.ReferenceMatch


class CreateWcsWithSipTestCaseSourceMatch(BaseTestCase, SideLoadTestCases, lsst.utils.tests.TestCase):
    MatchClass = afwTable.SourceMatch


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
