#
# LSST Data Management System
# Copyright 2016-2017 LSST/AURA.
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

import numpy as np

import lsst.utils.tests
import lsst.pex.exceptions
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.math
from lsst.afw.fits import readMetadata
from lsst.afw.geom.testUtils import makeSipIwcToPixel, makeSipPixelToIwc
from lsst.meas.astrom import (
    PolynomialTransform,
    ScaledPolynomialTransform,
    SipForwardTransform,
    SipReverseTransform,
    ScaledPolynomialTransformFitter,
    transformWcsPixels,
    rotateWcsPixelsBy90
)
from lsst.afw.geom.wcsUtils import getSipMatrixFromMetadata, getCdMatrixFromMetadata


def makeRandomCoefficientMatrix(n):
    matrix = np.random.randn(n, n)
    for i in range(1, n):
        matrix[i, (n-i):] = 0
    return matrix


def makeRandomAffineTransform():
    return lsst.afw.geom.AffineTransform(
        lsst.afw.geom.LinearTransform(np.random.randn(2, 2)),
        lsst.afw.geom.Extent2D(*np.random.randn(2))
    )


def makeRandomPolynomialTransform(order, sip=False):
    xc = makeRandomCoefficientMatrix(order + 1)
    yc = makeRandomCoefficientMatrix(order + 1)
    if sip:
        xc[0, 0] = 0
        yc[0, 0] = 0
        xc[0, 1] = 0
        yc[0, 1] = 0
        xc[1, 0] = 0
        yc[1, 0] = 0
    return PolynomialTransform(xc, yc)


def makeRandomScaledPolynomialTransform(order):
    return ScaledPolynomialTransform(
        makeRandomPolynomialTransform(order),
        makeRandomAffineTransform(),
        makeRandomAffineTransform()
    )


def makeRandomSipForwardTransform(order):
    return SipForwardTransform(
        lsst.afw.geom.Point2D(*np.random.randn(2)),
        lsst.afw.geom.LinearTransform(np.random.randn(2, 2)),
        makeRandomPolynomialTransform(order, sip=True)
    )


def makeRandomSipReverseTransform(order):
    origin = lsst.afw.geom.Point2D(*np.random.randn(2))
    cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
    poly = makeRandomPolynomialTransform(order, sip=False)
    return SipReverseTransform(origin, cd, poly)


class TransformTestMixin(object):

    def makeRandom(self):
        """Create an instance of the transform being tested with random testing.
        """
        raise NotImplementedError()

    def assertTransformsAlmostEqual(self, a, b, minval=0.0, maxval=1.0, atol=0, rtol=1E-8):
        rangeval = maxval - minval
        aArr = []
        bArr = []
        for i in range(10):
            xy = rangeval * np.random.rand(2) - minval
            point = lsst.afw.geom.Point2D(*xy)
            aArr.append(list(a(point)))
            bArr.append(list(b(point)))
        self.assertFloatsAlmostEqual(np.array(aArr), np.array(bArr), atol=atol, rtol=rtol)

    def testLinearize(self):
        """Test that the AffineTransform returned by linearize() is equivalent
        to the transform at the expansion point, and matches finite differences.
        """
        transform = self.makeRandom()
        point = lsst.afw.geom.Point2D(*np.random.randn(2))
        affine = transform.linearize(point)
        self.assertFloatsAlmostEqual(np.array(transform(point)), np.array(affine(point)), rtol=1E-14)
        delta = 1E-4
        deltaX = lsst.afw.geom.Extent2D(delta, 0.0)
        deltaY = lsst.afw.geom.Extent2D(0.0, delta)
        dtdx = (transform(point + deltaX) - transform(point - deltaX)) / (2*delta)
        dtdy = (transform(point + deltaY) - transform(point - deltaY)) / (2*delta)
        self.assertFloatsAlmostEqual(affine[affine.XX], dtdx.getX(), rtol=1E-6)
        self.assertFloatsAlmostEqual(affine[affine.YX], dtdx.getY(), rtol=1E-6)
        self.assertFloatsAlmostEqual(affine[affine.XY], dtdy.getX(), rtol=1E-6)
        self.assertFloatsAlmostEqual(affine[affine.YY], dtdy.getY(), rtol=1E-6)


class PolynomialTransformTestCase(lsst.utils.tests.TestCase, TransformTestMixin):

    def setUp(self):
        np.random.seed(50)

    def makeRandom(self):
        return makeRandomPolynomialTransform(4)

    def testArrayConstructor(self):
        """Test that construction with coefficient arrays yields an object with
        copies of those arrays, and that all dimensions must be the same.
        """
        order = 3
        xc = makeRandomCoefficientMatrix(order + 1)
        yc = makeRandomCoefficientMatrix(order + 1)
        p = PolynomialTransform(xc, yc)
        self.assertEqual(p.getOrder(), order)
        self.assertFloatsAlmostEqual(p.getXCoeffs(), xc, atol=0, rtol=0)
        self.assertFloatsAlmostEqual(p.getYCoeffs(), yc, atol=0, rtol=0)
        # Test that the coefficients are not a view.
        old = xc[0, 0]
        xc[0, 0] += 100.0
        self.assertEqual(p.getXCoeffs()[0, 0], old)
        # Test that rectangular coefficient arrays are not allowed.
        self.assertRaises(
            lsst.pex.exceptions.LengthError,
            PolynomialTransform,
            np.zeros((5, 4), dtype=float),
            np.zeros((5, 4), dtype=float)
        )
        # Test that x and y coefficient arrays must have the same shape.
        self.assertRaises(
            lsst.pex.exceptions.LengthError,
            PolynomialTransform,
            np.zeros((5, 5), dtype=float),
            np.zeros((4, 4), dtype=float)
        )

    def testConvertScaledPolynomial(self):
        """Test that we can convert a ScaledPolynomialTransform to a PolynomialTransform.
        """
        scaled = makeRandomScaledPolynomialTransform(4)
        converted = PolynomialTransform.convert(scaled)
        self.assertTransformsAlmostEqual(scaled, converted)

    def testConvertSipForward(self):
        """Test that we can convert a SipForwardTransform to a PolynomialTransform.
        """
        sipForward = makeRandomSipForwardTransform(4)
        converted = PolynomialTransform.convert(sipForward)
        self.assertTransformsAlmostEqual(sipForward, converted)

    def testConvertSipReverse(self):
        """Test that we can convert a SipForwardTransform to a PolynomialTransform.
        """
        sipReverse = makeRandomSipReverseTransform(4)
        converted = PolynomialTransform.convert(sipReverse)
        self.assertTransformsAlmostEqual(sipReverse, converted)

    def testCompose(self):
        """Test that AffineTransforms and PolynomialTransforms can be composed
        into an equivalent PolynomialTransform.
        """
        poly = makeRandomPolynomialTransform(4)
        affine = lsst.afw.geom.AffineTransform(
            lsst.afw.geom.LinearTransform(np.random.randn(2, 2)),
            lsst.afw.geom.Extent2D(*np.random.randn(2))
        )
        composed1 = lsst.meas.astrom.compose(poly, affine)
        composed2 = lsst.meas.astrom.compose(affine, poly)
        self.assertTransformsAlmostEqual(composed1, lambda p: poly(affine(p)))
        self.assertTransformsAlmostEqual(composed2, lambda p: affine(poly(p)))
        # Test that composition with an identity transform is a no-op
        composed3 = lsst.meas.astrom.compose(poly, lsst.afw.geom.AffineTransform())
        composed4 = lsst.meas.astrom.compose(lsst.afw.geom.AffineTransform(), poly)
        self.assertFloatsAlmostEqual(composed3.getXCoeffs(), poly.getXCoeffs())
        self.assertFloatsAlmostEqual(composed3.getYCoeffs(), poly.getYCoeffs())
        self.assertFloatsAlmostEqual(composed4.getXCoeffs(), poly.getXCoeffs())
        self.assertFloatsAlmostEqual(composed4.getYCoeffs(), poly.getYCoeffs())


class ScaledPolynomialTransformTestCase(lsst.utils.tests.TestCase, TransformTestMixin):

    def setUp(self):
        np.random.seed(50)

    def makeRandom(self):
        return makeRandomScaledPolynomialTransform(4)

    def testConstruction(self):
        poly = makeRandomPolynomialTransform(4)
        inputScaling = makeRandomAffineTransform()
        outputScalingInverse = makeRandomAffineTransform()
        scaled = ScaledPolynomialTransform(poly, inputScaling, outputScalingInverse)
        self.assertTransformsAlmostEqual(
            scaled,
            lambda p: outputScalingInverse(poly(inputScaling(p)))
        )

    def testConvertPolynomial(self):
        """Test that we can convert a PolynomialTransform to a ScaledPolynomialTransform.
        """
        poly = makeRandomPolynomialTransform(4)
        converted = ScaledPolynomialTransform.convert(poly)
        self.assertTransformsAlmostEqual(poly, converted)

    def testConvertSipForward(self):
        """Test that we can convert a SipForwardTransform to a ScaledPolynomialTransform.
        """
        sipForward = makeRandomSipForwardTransform(4)
        converted = ScaledPolynomialTransform.convert(sipForward)
        self.assertTransformsAlmostEqual(sipForward, converted)

    def testConvertSipReverse(self):
        """Test that we can convert a SipReverseTransform to a ScaledPolynomialTransform.
        """
        sipReverse = makeRandomSipReverseTransform(4)
        converted = ScaledPolynomialTransform.convert(sipReverse)
        self.assertTransformsAlmostEqual(sipReverse, converted)


class SipForwardTransformTestCase(lsst.utils.tests.TestCase, TransformTestMixin):

    def setUp(self):
        np.random.seed(50)

    def makeRandom(self):
        return makeRandomSipForwardTransform(4)

    def testConstruction(self):
        poly = makeRandomPolynomialTransform(4, sip=True)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = SipForwardTransform(crpix, cd, poly)
        self.assertTransformsAlmostEqual(
            sip,
            lambda p: cd((p - crpix) + poly(lsst.afw.geom.Point2D(p - crpix)))
        )

    def testConvertPolynomial(self):
        poly = makeRandomPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = lsst.meas.astrom.SipForwardTransform.convert(poly, crpix, cd)
        self.assertTransformsAlmostEqual(sip, poly)

    def testConvertScaledPolynomialManual(self):
        scaled = makeRandomScaledPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = lsst.meas.astrom.SipForwardTransform.convert(scaled, crpix, cd)
        self.assertTransformsAlmostEqual(sip, scaled)

    def testConvertScaledPolynomialAutomatic(self):
        scaled = makeRandomScaledPolynomialTransform(4)
        sip = lsst.meas.astrom.SipForwardTransform.convert(scaled)
        self.assertTransformsAlmostEqual(sip, scaled)

    def testTransformPixels(self):
        sip = makeRandomSipForwardTransform(4)
        affine = makeRandomAffineTransform()
        self.assertTransformsAlmostEqual(
            sip.transformPixels(affine),
            lambda p: sip(affine.invert()(p))
        )

    def testMakeWcs(self):
        """Test SipForwardTransform, SipReverseTransform and makeWcs
        """
        filename = os.path.join(os.path.dirname(__file__),
                                'imgCharSources-v85501867-R01-S00.sipheader')
        sipMetadata = readMetadata(filename)
        # We're building an ICRS-based TAN-SIP using coefficients read from metadata
        # so ignore the RADESYS in metadata (which is missing anyway, falling back to FK5)
        sipMetadata.set("RADESYS", "ICRS")
        crpix = lsst.afw.geom.Point2D(
            sipMetadata.get("CRPIX1") - 1,
            sipMetadata.get("CRPIX2") - 1,
        )
        crval = lsst.afw.geom.SpherePoint(
            sipMetadata.get("CRVAL1"),
            sipMetadata.get("CRVAL2"), lsst.afw.geom.degrees,
        )
        cdLinearTransform = lsst.afw.geom.LinearTransform(getCdMatrixFromMetadata(sipMetadata))
        aArr = getSipMatrixFromMetadata(sipMetadata, "A")
        bArr = getSipMatrixFromMetadata(sipMetadata, "B")
        apArr = getSipMatrixFromMetadata(sipMetadata, "AP")
        bpArr = getSipMatrixFromMetadata(sipMetadata, "BP")
        abPoly = PolynomialTransform(aArr, bArr)
        abRevPoly = PolynomialTransform(apArr, bpArr)
        fwd = SipForwardTransform(crpix, cdLinearTransform, abPoly)
        rev = SipReverseTransform(crpix, cdLinearTransform, abRevPoly)
        wcsFromMakeWcs = lsst.meas.astrom.makeWcs(fwd, rev, crval)
        wcsFromMetadata = lsst.afw.geom.makeSkyWcs(sipMetadata, strip=False)

        # Check SipForwardTransform against a local implementation
        localPixelToIwc = makeSipPixelToIwc(sipMetadata)
        self.assertTransformsAlmostEqual(fwd, localPixelToIwc.applyForward, maxval=2000)

        # Compare SipReverseTransform against a local implementation
        # Use the forward direction first to get sensible inputs
        localIwcToPixel = makeSipIwcToPixel(sipMetadata)

        def fwdThenRev(p):
            return rev(fwd(p))

        def fwdThenLocalRev(p):
            return localIwcToPixel.applyForward(fwd(p))

        self.assertTransformsAlmostEqual(fwdThenRev, fwdThenLocalRev, maxval=2000)

        # Check that SipReverseTransform is the inverse of SipForwardTransform;
        # this is not perfect because the coefficients don't define a perfect inverse
        def nullTransform(p):
            return p

        self.assertTransformsAlmostEqual(fwdThenRev, nullTransform, maxval=2000, atol=1e-3)

        # Check SipForwardTransform against the one contained in wcsFromMakeWcs
        # (Don't bother with the other direction because the WCS transform is iterative,
        # so it doesn't tell us anything useful about SipReverseTransform
        pixelToIwc = lsst.afw.geom.getPixelToIntermediateWorldCoords(wcsFromMetadata)
        self.assertTransformsAlmostEqual(fwd, pixelToIwc.applyForward, maxval=2000)

        # Check a WCS constructed from SipForwardTransform, SipReverseTransform
        # against one constructed directly from the metadata
        bbox = lsst.afw.geom.Box2D(lsst.afw.geom.Point2D(0, 0), lsst.afw.geom.Extent2D(2000, 2000))
        self.assertWcsAlmostEqualOverBBox(wcsFromMakeWcs, wcsFromMetadata, bbox)

    def testTransformWcsPixels(self):
        filename = os.path.join(os.path.dirname(__file__),
                                'imgCharSources-v85501867-R01-S00.sipheader')
        wcs1 = lsst.afw.geom.makeSkyWcs(readMetadata(filename))
        s = makeRandomAffineTransform()
        wcs2 = transformWcsPixels(wcs1, s)
        crvalDeg = wcs1.getSkyOrigin().getPosition(lsst.afw.geom.degrees)

        def t1a(p):
            raDeg, decDeg = crvalDeg + lsst.afw.geom.Extent2D(p)
            sky = lsst.afw.geom.SpherePoint(raDeg, decDeg, lsst.afw.geom.degrees)
            return s(wcs1.skyToPixel(sky))

        def t2a(p):
            raDeg, decDeg = crvalDeg + lsst.afw.geom.Extent2D(p)
            sky = lsst.afw.geom.SpherePoint(raDeg, decDeg, lsst.afw.geom.degrees)
            return wcs2.skyToPixel(sky)

        self.assertTransformsAlmostEqual(t1a, t2a)

        def t1b(p):
            sky = wcs1.pixelToSky(s.invert()(p))
            return sky.getPosition(lsst.afw.geom.degrees)

        def t2b(p):
            sky = wcs2.pixelToSky(p)
            return sky.getPosition(lsst.afw.geom.degrees)

        self.assertTransformsAlmostEqual(t1b, t2b)

    def testRotateWcsPixelsBy90(self):
        filename = os.path.join(os.path.dirname(__file__),
                                'imgCharSources-v85501867-R01-S00.sipheader')
        wcs0 = lsst.afw.geom.makeSkyWcs(readMetadata(filename))
        w, h = 11, 12
        image0 = lsst.afw.image.ImageD(w, h)
        x, y = np.meshgrid(np.arange(w), np.arange(h))
        # Make a slowly-varying image of an asymmetric function
        image0.getArray()[:, :] = (x/w)**2 + 0.5*(x/w)*(y/h) - 3.0*(y/h)**2
        dimensions = image0.getBBox().getDimensions()

        image1 = lsst.afw.math.rotateImageBy90(image0, 1)
        wcs1 = rotateWcsPixelsBy90(wcs0, 1, dimensions)
        image2 = lsst.afw.math.rotateImageBy90(image0, 2)
        wcs2 = rotateWcsPixelsBy90(wcs0, 2, dimensions)
        image3 = lsst.afw.math.rotateImageBy90(image0, 3)
        wcs3 = rotateWcsPixelsBy90(wcs0, 3, dimensions)

        bbox = image0.getBBox()
        image0r = lsst.afw.image.ImageD(bbox)
        image1r = lsst.afw.image.ImageD(bbox)
        image2r = lsst.afw.image.ImageD(bbox)
        image3r = lsst.afw.image.ImageD(bbox)

        ctrl = lsst.afw.math.WarpingControl("nearest")
        lsst.afw.math.warpImage(image0r, wcs0, image0, wcs0, ctrl)
        lsst.afw.math.warpImage(image1r, wcs0, image1, wcs1, ctrl)
        lsst.afw.math.warpImage(image2r, wcs0, image2, wcs2, ctrl)
        lsst.afw.math.warpImage(image3r, wcs0, image3, wcs3, ctrl)

        # warpImage doesn't seem to handle the first row and column,
        # even with nearest-neighbor interpolation, so we have to
        # ignore pixels it didn't know how to populate.
        def compareFinite(ref, target):
            finitPixels = np.isfinite(target.getArray())
            self.assertGreater(finitPixels.sum(), 0.7*target.getArray().size)
            self.assertFloatsAlmostEqual(
                ref.getArray()[finitPixels],
                target.getArray()[finitPixels],
                rtol=1E-6
            )

        compareFinite(image0, image0r)
        compareFinite(image0, image1r)
        compareFinite(image0, image2r)
        compareFinite(image0, image3r)


class SipReverseTransformTestCase(lsst.utils.tests.TestCase, TransformTestMixin):

    def setUp(self):
        np.random.seed(50)

    def makeRandom(self):
        return makeRandomSipReverseTransform(4)

    def testConstruction(self):
        poly = makeRandomPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = SipReverseTransform(crpix, cd, poly)
        offset = lsst.afw.geom.Extent2D(crpix)
        cdInverse = cd.invert()
        self.assertTransformsAlmostEqual(
            sip,
            lambda p: offset + lsst.afw.geom.Extent2D(cdInverse(p)) + poly(cdInverse(p))
        )

    def testConvertPolynomial(self):
        poly = makeRandomPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = lsst.meas.astrom.SipReverseTransform.convert(poly, crpix, cd)
        self.assertTransformsAlmostEqual(sip, poly)

    def testConvertScaledPolynomialManual(self):
        scaled = makeRandomScaledPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = lsst.meas.astrom.SipReverseTransform.convert(scaled, crpix, cd)
        self.assertTransformsAlmostEqual(sip, scaled)

    def testConvertScaledPolynomialAutomatic(self):
        scaled = makeRandomScaledPolynomialTransform(4)
        sip = lsst.meas.astrom.SipReverseTransform.convert(scaled)
        self.assertTransformsAlmostEqual(sip, scaled)

    def testTransformPixels(self):
        sip = makeRandomSipReverseTransform(4)
        affine = makeRandomAffineTransform()
        self.assertTransformsAlmostEqual(
            sip.transformPixels(affine),
            lambda p: affine(sip(p))
        )


class ScaledPolynomialTransformFitterTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        np.random.seed(50)

    def testFromMatches(self):
        # Setup artifical matches that correspond to a known (random) PolynomialTransform.
        order = 3
        truePoly = makeRandomPolynomialTransform(order)
        crval = lsst.afw.geom.SpherePoint(35.0, 10.0, lsst.afw.geom.degrees)
        crpix = lsst.afw.geom.Point2D(50, 50)
        cd = lsst.afw.geom.LinearTransform.makeScaling((0.2*lsst.afw.geom.arcseconds).asDegrees()).getMatrix()
        initialWcs = lsst.afw.geom.makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cd)
        bbox = lsst.afw.geom.Box2D(
            crval.getPosition(lsst.afw.geom.arcseconds) - lsst.afw.geom.Extent2D(20, 20),
            crval.getPosition(lsst.afw.geom.arcseconds) + lsst.afw.geom.Extent2D(20, 20),
        )
        srcSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        srcPosKey = lsst.afw.table.Point2DKey.addFields(srcSchema, "pos", "source position", "pix")
        srcErrKey = lsst.afw.table.CovarianceMatrix2fKey.addFields(srcSchema, "pos",
                                                                   ["x", "y"], ["pix", "pix"])
        srcSchema.getAliasMap().set("slot_Centroid", "pos")
        nPoints = 10
        trueSrc = lsst.afw.table.SourceCatalog(srcSchema)
        trueSrc.reserve(nPoints)
        measSrc = lsst.afw.table.SourceCatalog(srcSchema)
        measSrc.reserve(nPoints)
        ref = lsst.afw.table.SimpleCatalog(lsst.afw.table.SimpleTable.makeMinimalSchema())
        ref.reserve(nPoints)
        refCoordKey = ref.getCoordKey()
        errScaling = 1E-14
        matches = []
        initialIwcToSky = lsst.afw.geom.getIntermediateWorldCoordsToSky(initialWcs)
        for i in range(nPoints):
            refRec = ref.addNew()
            raDeg, decDeg = (
                np.random.uniform(low=bbox.getMinX(), high=bbox.getMaxX()),
                np.random.uniform(low=bbox.getMinY(), high=bbox.getMaxY()),
            )
            skyCoord = lsst.afw.geom.SpherePoint(raDeg, decDeg, lsst.afw.geom.arcseconds)
            refRec.set(refCoordKey, skyCoord)
            trueRec = trueSrc.addNew()
            truePos = truePoly(initialIwcToSky.applyInverse(skyCoord))
            trueRec.set(srcPosKey, truePos)
            measRec = measSrc.addNew()
            covSqrt = np.random.randn(3, 2)
            cov = (errScaling*(np.dot(covSqrt.transpose(), covSqrt) +
                   np.diag([1.0, 1.0]))).astype(np.float32)
            # We don't actually perturb positions according to noise level, as
            # this makes it much harder to test that the result agrees with
            # what we put in.
            measPos = truePos
            measRec.set(srcPosKey, measPos)
            measRec.set(srcErrKey, cov)
            match = lsst.afw.table.ReferenceMatch(refRec, measRec, (measPos - truePos).computeNorm())
            matches.append(match)
        # Construct a fitter, and verify that the internal catalog it constructs is what we expect.
        fitter = ScaledPolynomialTransformFitter.fromMatches(order, matches, initialWcs, 0.0)
        expected = lsst.meas.astrom.compose(
            fitter.getOutputScaling(),
            lsst.meas.astrom.compose(truePoly, fitter.getInputScaling().invert())
        )
        data = fitter.getData()
        dataOutKey = lsst.afw.table.Point2DKey(data.schema["src"])
        dataInKey = lsst.afw.table.Point2DKey(data.schema["intermediate"])
        dataErrKey = lsst.afw.table.CovarianceMatrix2fKey(data.schema["src"], ["x", "y"])
        scaledInBBox = lsst.afw.geom.Box2D()
        scaledOutBBox = lsst.afw.geom.Box2D()
        vandermonde = np.zeros((nPoints, (order + 1)*(order + 2)//2), dtype=float)
        for i, (match, dataRec, trueRec) in enumerate(zip(matches, data, trueSrc)):
            self.assertEqual(match.second.getX(), dataRec.get("src_x"))
            self.assertEqual(match.second.getY(), dataRec.get("src_y"))
            self.assertEqual(match.first.getId(), dataRec.get("ref_id"))
            self.assertEqual(match.second.getId(), dataRec.get("src_id"))
            refPos = initialIwcToSky.applyInverse(match.first.getCoord())
            self.assertEqual(refPos.getX(), dataRec.get("intermediate_x"))
            self.assertEqual(refPos.getY(), dataRec.get("intermediate_y"))
            self.assertFloatsAlmostEqual(match.second.get(srcErrKey), dataRec.get(dataErrKey), rtol=1E-7)
            scaledIn = fitter.getInputScaling()(dataRec.get(dataInKey))
            scaledOut = fitter.getOutputScaling()(dataRec.get(dataOutKey))
            scaledInBBox.include(scaledIn)
            scaledOutBBox.include(scaledOut)
            self.assertFloatsAlmostEqual(np.array(expected(scaledIn)), np.array(scaledOut), rtol=1E-7)
            j = 0
            for n in range(order + 1):
                for p in range(n + 1):
                    q = n - p
                    vandermonde[i, j] = scaledIn.getX()**p * scaledIn.getY()**q
                    j += 1
        # Verify that scaling transforms move inputs and outputs into [-1, 1]
        self.assertFloatsAlmostEqual(scaledInBBox.getMinX(), -1.0, rtol=1E-12)
        self.assertFloatsAlmostEqual(scaledInBBox.getMinY(), -1.0, rtol=1E-12)
        self.assertFloatsAlmostEqual(scaledInBBox.getMaxX(), 1.0, rtol=1E-12)
        self.assertFloatsAlmostEqual(scaledInBBox.getMaxY(), 1.0, rtol=1E-12)
        self.assertFloatsAlmostEqual(scaledOutBBox.getMinX(), -1.0, rtol=1E-12)
        self.assertFloatsAlmostEqual(scaledOutBBox.getMinY(), -1.0, rtol=1E-12)
        self.assertFloatsAlmostEqual(scaledOutBBox.getMaxX(), 1.0, rtol=1E-12)
        self.assertFloatsAlmostEqual(scaledOutBBox.getMaxY(), 1.0, rtol=1E-12)
        # Run the fitter, and check that we get out approximately what we put in.
        fitter.fit(order)
        fitter.updateModel()
        # Check the transformed input points.
        self.assertFloatsAlmostEqual(data.get("model_x"), trueSrc.getX(), rtol=1E-15)
        self.assertFloatsAlmostEqual(data.get("model_y"), trueSrc.getY(), rtol=1E-15)
        # Check the actual transform's coefficients (after composing in the scaling, which is
        # a lot of the reason we lose a lot of precision here).
        fittedPoly = lsst.meas.astrom.PolynomialTransform.convert(fitter.getTransform())
        self.assertFloatsAlmostEqual(fittedPoly.getXCoeffs(), truePoly.getXCoeffs(), rtol=1E-5, atol=1E-5)
        self.assertFloatsAlmostEqual(fittedPoly.getYCoeffs(), truePoly.getYCoeffs(), rtol=1E-5, atol=1E-5)

    def testFromGrid(self):
        outOrder = 8
        inOrder = 2
        toInvert = makeRandomScaledPolynomialTransform(inOrder)
        bbox = lsst.afw.geom.Box2D(lsst.afw.geom.Point2D(432, -671), lsst.afw.geom.Point2D(527, -463))
        fitter = ScaledPolynomialTransformFitter.fromGrid(outOrder, bbox, 50, 50, toInvert)
        fitter.fit(outOrder)
        fitter.updateModel()
        data = fitter.getData()
        result = fitter.getTransform()
        inputKey = lsst.afw.table.Point2DKey(data.schema["input"])
        outputKey = lsst.afw.table.Point2DKey(data.schema["output"])
        for record in data:
            self.assertFloatsAlmostEqual(np.array(record.get(inputKey)),
                                         np.array(toInvert(record.get(outputKey))))
            self.assertFloatsAlmostEqual(np.array(result(record.get(inputKey))),
                                         np.array(record.get(outputKey)),
                                         rtol=1E-2)  # even at much higher order, inverse can't be perfect.


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
