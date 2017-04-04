#!/usr/bin/env python
from __future__ import absolute_import, division, print_function
from builtins import range

#
# LSST Data Management System
# Copyright 2016 LSST/AURA.
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
import lsst.pex.exceptions
import lsst.afw.geom
from lsst.meas.astrom import (
    PolynomialTransform,
    ScaledPolynomialTransform,
    SipForwardTransform,
    SipReverseTransform,
    ScaledPolynomialTransformFitter
)


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

    def assertTransformsNearlyEqual(self, a, b, rtol=1E-8):
        for i in range(10):
            point = lsst.afw.geom.Point2D(*np.random.randn(2))
            aOut = a(point)
            bOut = b(point)
            self.assertClose(np.array(aOut), np.array(bOut), rtol=rtol)

    def testLinearize(self):
        """Test that the AffineTransform returned by linearize() is equivalent
        to the transform at the expansion point, and matches finite differences.
        """
        transform = self.makeRandom()
        point = lsst.afw.geom.Point2D(*np.random.randn(2))
        affine = transform.linearize(point)
        self.assertClose(np.array(transform(point)), np.array(affine(point)), rtol=1E-14)
        delta = 1E-4
        deltaX = lsst.afw.geom.Extent2D(delta, 0.0)
        deltaY = lsst.afw.geom.Extent2D(0.0, delta)
        dtdx = (transform(point + deltaX) - transform(point - deltaX)) / (2*delta)
        dtdy = (transform(point + deltaY) - transform(point - deltaY)) / (2*delta)
        self.assertClose(affine[affine.XX], dtdx.getX(), rtol=1E-6)
        self.assertClose(affine[affine.YX], dtdx.getY(), rtol=1E-6)
        self.assertClose(affine[affine.XY], dtdy.getX(), rtol=1E-6)
        self.assertClose(affine[affine.YY], dtdy.getY(), rtol=1E-6)


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
        self.assertClose(p.getXCoeffs(), xc, atol=0, rtol=0)
        self.assertClose(p.getYCoeffs(), yc, atol=0, rtol=0)
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
        self.assertTransformsNearlyEqual(scaled, converted)

    def testConvertSipForward(self):
        """Test that we can convert a SipForwardTransform to a PolynomialTransform.
        """
        sipForward = makeRandomSipForwardTransform(4)
        converted = PolynomialTransform.convert(sipForward)
        self.assertTransformsNearlyEqual(sipForward, converted)

    def testConvertSipReverse(self):
        """Test that we can convert a SipForwardTransform to a PolynomialTransform.
        """
        sipReverse = makeRandomSipReverseTransform(4)
        converted = PolynomialTransform.convert(sipReverse)
        self.assertTransformsNearlyEqual(sipReverse, converted)

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
        self.assertTransformsNearlyEqual(composed1, lambda p: poly(affine(p)))
        self.assertTransformsNearlyEqual(composed2, lambda p: affine(poly(p)))
        # Test that composition with an identity transform is a no-op
        composed3 = lsst.meas.astrom.compose(poly, lsst.afw.geom.AffineTransform())
        composed4 = lsst.meas.astrom.compose(lsst.afw.geom.AffineTransform(), poly)
        self.assertClose(composed3.getXCoeffs(), poly.getXCoeffs())
        self.assertClose(composed3.getYCoeffs(), poly.getYCoeffs())
        self.assertClose(composed4.getXCoeffs(), poly.getXCoeffs())
        self.assertClose(composed4.getYCoeffs(), poly.getYCoeffs())


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
        self.assertTransformsNearlyEqual(
            scaled,
            lambda p: outputScalingInverse(poly(inputScaling(p)))
        )

    def testConvertPolynomial(self):
        """Test that we can convert a PolynomialTransform to a ScaledPolynomialTransform.
        """
        poly = makeRandomPolynomialTransform(4)
        converted = ScaledPolynomialTransform.convert(poly)
        self.assertTransformsNearlyEqual(poly, converted)

    def testConvertSipForward(self):
        """Test that we can convert a SipForwardTransform to a ScaledPolynomialTransform.
        """
        sipForward = makeRandomSipForwardTransform(4)
        converted = ScaledPolynomialTransform.convert(sipForward)
        self.assertTransformsNearlyEqual(sipForward, converted)

    def testConvertSipReverse(self):
        """Test that we can convert a SipReverseTransform to a ScaledPolynomialTransform.
        """
        sipReverse = makeRandomSipReverseTransform(4)
        converted = ScaledPolynomialTransform.convert(sipReverse)
        self.assertTransformsNearlyEqual(sipReverse, converted)


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
        self.assertTransformsNearlyEqual(
            sip,
            lambda p: cd((p - crpix) + poly(lsst.afw.geom.Point2D(p - crpix)))
        )

    def testConvertPolynomial(self):
        poly = makeRandomPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = lsst.meas.astrom.SipForwardTransform.convert(poly, crpix, cd)
        self.assertTransformsNearlyEqual(sip, poly)

    def testConvertScaledPolynomialManual(self):
        scaled = makeRandomScaledPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = lsst.meas.astrom.SipForwardTransform.convert(scaled, crpix, cd)
        self.assertTransformsNearlyEqual(sip, scaled)

    def testConvertScaledPolynomialAutomatic(self):
        scaled = makeRandomScaledPolynomialTransform(4)
        sip = lsst.meas.astrom.SipForwardTransform.convert(scaled)
        self.assertTransformsNearlyEqual(sip, scaled)

    def testMakeWcs(self):
        fwd = self.makeRandom()
        rev = SipReverseTransform(  # this isn't actually the inverse of fwd, but that doesn't matter
            fwd.getPixelOrigin(),
            fwd.getCDMatrix(),
            makeRandomPolynomialTransform(4)
        )
        crval = lsst.afw.geom.Point2D(45.0, 45.0)
        wcs = lsst.meas.astrom.makeWcs(
            fwd, rev,
            lsst.afw.coord.IcrsCoord(crval, lsst.afw.geom.degrees)
        )
        # We can only test agreement with TanWcs in one direction, because
        # TanWcs doesn't provide an inverse to skyToIntermediateWorldCoord.

        def t1(p):
            sky = lsst.afw.coord.IcrsCoord(crval + lsst.afw.geom.Extent2D(p), lsst.afw.geom.degrees)
            return rev(wcs.skyToIntermediateWorldCoord(sky))

        def t2(p):
            sky = lsst.afw.coord.IcrsCoord(crval + lsst.afw.geom.Extent2D(p), lsst.afw.geom.degrees)
            return wcs.skyToPixel(sky)

        self.assertTransformsNearlyEqual(t1, t2)


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
        self.assertTransformsNearlyEqual(
            sip,
            lambda p: offset + lsst.afw.geom.Extent2D(cdInverse(p)) + poly(cdInverse(p))
        )

    def testConvertPolynomial(self):
        poly = makeRandomPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = lsst.meas.astrom.SipReverseTransform.convert(poly, crpix, cd)
        self.assertTransformsNearlyEqual(sip, poly)

    def testConvertScaledPolynomialManual(self):
        scaled = makeRandomScaledPolynomialTransform(4)
        cd = lsst.afw.geom.LinearTransform(np.random.randn(2, 2))
        crpix = lsst.afw.geom.Point2D(*np.random.randn(2))
        sip = lsst.meas.astrom.SipReverseTransform.convert(scaled, crpix, cd)
        self.assertTransformsNearlyEqual(sip, scaled)

    def testConvertScaledPolynomialAutomatic(self):
        scaled = makeRandomScaledPolynomialTransform(4)
        sip = lsst.meas.astrom.SipReverseTransform.convert(scaled)
        self.assertTransformsNearlyEqual(sip, scaled)


class ScaledPolynomialTransformFitterTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        np.random.seed(50)

    def testFromMatches(self):
        # Setup artifical matches that correspond to a known (random) PolynomialTransform.
        order = 3
        truePoly = makeRandomPolynomialTransform(order)
        crval = lsst.afw.coord.IcrsCoord(lsst.afw.geom.Point2D(35.0, 10.0), lsst.afw.geom.degrees)
        crpix = lsst.afw.geom.Point2D(50, 50)
        cd = lsst.afw.geom.LinearTransform.makeScaling((0.2*lsst.afw.geom.arcseconds).asDegrees())
        initialWcs = lsst.afw.image.makeWcs(crval, crpix, cd[cd.XX], cd[cd.XY], cd[cd.YX], cd[cd.YY])
        bbox = lsst.afw.geom.Box2D(
            (lsst.afw.geom.Point2D(crval.getPosition(lsst.afw.geom.arcseconds)) -
                lsst.afw.geom.Extent2D(20, 20)),
            (lsst.afw.geom.Point2D(crval.getPosition(lsst.afw.geom.arcseconds)) +
                lsst.afw.geom.Extent2D(20, 20))
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
        for i in range(nPoints):
            refRec = ref.addNew()
            skyPos = lsst.afw.geom.Point2D(
                np.random.uniform(low=bbox.getMinX(), high=bbox.getMaxX()),
                np.random.uniform(low=bbox.getMinY(), high=bbox.getMaxY())
            )
            skyCoord = lsst.afw.coord.IcrsCoord(skyPos, lsst.afw.geom.arcseconds)
            refRec.set(refCoordKey, skyCoord)
            trueRec = trueSrc.addNew()
            truePos = truePoly(initialWcs.skyToIntermediateWorldCoord(skyCoord))
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
            refPos = initialWcs.skyToIntermediateWorldCoord(match.first.getCoord())
            self.assertEqual(refPos.getX(), dataRec.get("intermediate_x"))
            self.assertEqual(refPos.getY(), dataRec.get("intermediate_y"))
            self.assertClose(match.second.get(srcErrKey), dataRec.get(dataErrKey), rtol=1E-7)
            scaledIn = fitter.getInputScaling()(dataRec.get(dataInKey))
            scaledOut = fitter.getOutputScaling()(dataRec.get(dataOutKey))
            scaledInBBox.include(scaledIn)
            scaledOutBBox.include(scaledOut)
            self.assertClose(np.array(expected(scaledIn)), np.array(scaledOut), rtol=1E-7)
            j = 0
            for n in range(order + 1):
                for p in range(n + 1):
                    q = n - p
                    vandermonde[i, j] = scaledIn.getX()**p * scaledIn.getY()**q
                    j += 1
        # Verify that scaling transforms move inputs and outputs into [-1, 1]
        self.assertClose(scaledInBBox.getMinX(), -1.0, rtol=1E-12)
        self.assertClose(scaledInBBox.getMinY(), -1.0, rtol=1E-12)
        self.assertClose(scaledInBBox.getMaxX(), 1.0, rtol=1E-12)
        self.assertClose(scaledInBBox.getMaxY(), 1.0, rtol=1E-12)
        self.assertClose(scaledOutBBox.getMinX(), -1.0, rtol=1E-12)
        self.assertClose(scaledOutBBox.getMinY(), -1.0, rtol=1E-12)
        self.assertClose(scaledOutBBox.getMaxX(), 1.0, rtol=1E-12)
        self.assertClose(scaledOutBBox.getMaxY(), 1.0, rtol=1E-12)
        # Run the fitter, and check that we get out approximately what we put in.
        fitter.fit(order)
        fitter.updateModel()
        # Check the transformed input points.
        self.assertClose(data.get("model_x"), trueSrc.getX(), rtol=1E-15)
        self.assertClose(data.get("model_y"), trueSrc.getY(), rtol=1E-15)
        # Check the actual transform's coefficients (after composing in the scaling, which is
        # a lot of the reason we lose a lot of precision here).
        fittedPoly = lsst.meas.astrom.PolynomialTransform.convert(fitter.getTransform())
        self.assertClose(fittedPoly.getXCoeffs(), truePoly.getXCoeffs(), rtol=1E-5, atol=1E-5)
        self.assertClose(fittedPoly.getYCoeffs(), truePoly.getYCoeffs(), rtol=1E-5, atol=1E-5)

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
            self.assertClose(np.array(record.get(inputKey)), np.array(toInvert(record.get(outputKey))))
            self.assertClose(
                np.array(result(record.get(inputKey))),
                np.array(record.get(outputKey)),
                rtol=1E-2  # even at much higher order, inverse can't be perfect.
            )


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
