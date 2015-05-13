#!/usr/bin/env python
import math
import unittest

import numpy
import matplotlib.pylab as pylab

import lsst.utils.tests as tests
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.astrom.sip import makeCreateWcsWithSip

def approximateWcs(wcs, bbox, order=3, nx=20, ny=20, useTanWcs=False):
    """Approximate an existing WCS as a TAN-SIP WCS

    The fit is performed by evaluating the WCS at a uniform grid of points within a bounding box.

    @param[in] wcs  wcs to approximate
    @param[in] bbox  the region over which the WCS will be fit
    @param[in] order  order of SIP fit
    @param[in] nx  number of grid points along x
    @param[in] ny  number of grid points along y
    @param[in] useTanWcs  send a TAN version of wcs to the fitter? It is documented to require that,
        but I don't think the fitter actually cares
    @return the fit TAN-SIP WCS
    """
    if useTanWcs:
        crCoord = wcs.getSkyOrigin()
        crPix = wcs.getPixelOrigin()
        cdMat = wcs.getCDMatrix()
        tanWcs = afwImage.makeWcs(crCoord, crPix, cdMat[0,0], cdMat[0,1], cdMat[1,0], cdMat[1,1])
    else:
        tanWcs = wcs
    
    # create a matchList consisting of a grid of points covering the bbox
    refSchema = afwTable.SimpleTable.makeMinimalSchema()
    refCoordKey = refSchema["coord"].asKey()
    refCat = afwTable.SimpleCatalog(refSchema)

    sourceSchema = afwTable.SourceTable.makeMinimalSchema()
    SingleFrameMeasurementTask(schema=sourceSchema) # expand the schema
    sourceCentroidKey = afwTable.Point2DKey(sourceSchema["slot_Centroid"])
        
    sourceCat = afwTable.SourceCatalog(sourceSchema)

    matchList = afwTable.ReferenceMatchVector()

    bboxd = afwGeom.Box2D(bbox)
    for x in numpy.linspace(bboxd.getMinX(), bboxd.getMaxX(), nx):
        for y in numpy.linspace(bboxd.getMinY(), bboxd.getMaxY(), ny):
            pixelPos = afwGeom.Point2D(x, y)
            skyCoord = wcs.pixelToSky(pixelPos)

            refObj = refCat.addNew()
            refObj.set(refCoordKey, skyCoord)

            source = sourceCat.addNew()
            source.set(sourceCentroidKey, pixelPos)

            matchList.append(afwTable.ReferenceMatch(refObj, source, 0.0))
            
    # The TAN-SIP fitter is fitting x and y separately, so we have to iterate to make it converge 
    for indx in range(3) :
        sipObject = makeCreateWcsWithSip(matchList, tanWcs, order, bbox)
        tanWcs = sipObject.getNewWcs()
    return sipObject.getNewWcs()


class ApproximateWcsTestCase(unittest.TestCase):
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
        metadata.set("EQUINOX",       2000.0)
        metadata.setDouble("CRVAL1",   215.60)
        metadata.setDouble("CRVAL2",    53.16)
        metadata.setDouble("CRPIX1",  self.crPix[0])
        metadata.setDouble("CRPIX2",  self.crPix[1])
        metadata.set("CTYPE1", "RA---TAN")
        metadata.set("CTYPE2", "DEC--TAN")
        metadata.setDouble("CD1_1", 5.10808596133527E-05)
        metadata.setDouble("CD1_2", 1.85579539217196E-07)
        metadata.setDouble("CD2_2", -5.10281493481982E-05)
        metadata.setDouble("CD2_1", -8.27440751733828E-07)
        self.tanWcs = afwImage.cast_TanWcs(afwImage.makeWcs(metadata))

    def tearDown(self):
        del self.tanWcs

    def testTrivial(self):
        """Add no distortion"""
        for order in (3, 4, 5, 6):
            self.doTest("testTrivial", afwGeom.IdentityXYTransform(), order=order, doPlot=False)

    def testRadial(self):
        """Add a radial transform"""
        for order in (4, 5, 6):
            self.doTest("testRadial", afwGeom.RadialXYTransform([0, 1.001, 0.000003]), order=order, doPlot=False)

    def doTest(self, name, xyTransform, order=3, doPlot=False):
        """Create a DistortedTanWcs from the specified transform and fit it
        """
        wcs = afwImage.DistortedTanWcs(self.tanWcs, xyTransform)

        fitWcs = approximateWcs(
            wcs = wcs,
            bbox = self.bbox,
            order=order,
        )

        if doPlot:
            self.plotWcs(wcs, fitWcs, self.bbox, xyTransform)

        try:
            msg = "ERROR: %s failed with order %s" % (name, order)
            self.assertWcssAlmostEqual(wcs, fitWcs, self.bbox, msg=msg)
            print "OK: %s succeeded with order %s" % (name, order)
        except Exception, e:
            print e

    def plotWcs(self, wcs0, wcs1, bbox, xyTransform):
        bboxd = afwGeom.Box2D(bbox)
        x0Arr=[]
        y0Arr=[]
        x1Arr=[]
        y1Arr=[]
        x2Arr=[]
        y2Arr=[]
        for x in numpy.linspace(bboxd.getMinX(), bboxd.getMaxX(), 10):
            for y in numpy.linspace(bboxd.getMinY(), bboxd.getMaxY(), 10):
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
        pylab.plot(x0Arr, y0Arr, 'b+', x1Arr, y1Arr, 'rx', x2Arr, y2Arr, 'g.')

        pylab.show()

    def assertWcssAlmostEqual(self, wcs0, wcs1, bbox, maxSkyErr=0.01 * afwGeom.arcseconds, maxPixErr=0.02,
        nx=5, ny=5, msg="WCSs differ"):
        """Assert that two WCSs give nearly equal results for pixelToSky and skyToPixel

        @param[in] wcs0  WCS 0 (an lsst.afw.image.Wcs)
        @param[in] wcs1  WCS 1 (an lsst.afw.image.Wcs)
        @param[in] bbox  boundaries of pixel grid over which to compare the WCSs (an lsst.afw.geom.Box2I or Box2D)
        @param[in] maxSkyErr  maximum separation between sky positions computed using Wcs.pixelToSky
            (an lsst.afw.geom.Angle)
        @param[in] maxPixErr  maximum separation between pixel positions computed using Wcs.skyToPixel
        @param[in] nx  number of points in x for the grid of pixel positions
        @param[in] ny  number of points in y for the grid of pixel positions
        @param[in] msg  prefix for error messages; details of the error are appended after ": "

        @throw AssertionError if the two WCSs do not match sufficiently closely
        """
        bboxd = afwGeom.Box2D(bbox)
        xList = numpy.linspace(bboxd.getMinX(), bboxd.getMaxX(), nx)
        yList = numpy.linspace(bboxd.getMinY(), bboxd.getMaxY(), ny)
        maxSkyErrPixPos = [afwGeom.Angle(0), None]
        maxPixErrSkyPos = [0, None]
        for x in xList:
            for y in yList:
                fromPixPos = afwGeom.Point2D(x, y)
                sky0 = wcs0.pixelToSky(fromPixPos)
                sky1 = wcs1.pixelToSky(fromPixPos)
                skyErr = sky0.angularSeparation(sky1)
                if skyErr > maxSkyErrPixPos[0]:
                    maxSkyErrPixPos = (skyErr, fromPixPos)

                toPixPos0 = wcs0.skyToPixel(sky0)
                toPixPos1 = wcs1.skyToPixel(sky0)
                pixErr = math.hypot(*(toPixPos0 - toPixPos1))
                if pixErr > maxPixErrSkyPos[0]:
                    maxPixErrSkyPos = (pixErr, sky0)

        errStrList = []
        if maxSkyErrPixPos[0] > maxSkyErr:
            skyErr, pixPos = maxSkyErrPixPos
            errStrList.append("%f arcsec sky error > %f arcsec max sky error at pixPos=%s" %
                (skyErr.asArcseconds(), maxSkyErr.asArcseconds(), pixPos))
        if maxPixErrSkyPos[0] > maxPixErr:
            pixErr, skyPos = maxPixErrSkyPos
            errStrList.append("%f pix error > %f max pix error at skyPos=%s" %
                    (pixErr, maxPixErr, skyPos))
        if errStrList:
            raise AssertionError("%s: %s" % (msg, "; ".join(errStrList)))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(ApproximateWcsTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
