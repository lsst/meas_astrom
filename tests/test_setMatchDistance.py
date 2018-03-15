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
from __future__ import absolute_import, division, print_function

from builtins import object
import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from lsst.meas.algorithms import LoadReferenceObjectsTask
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.astrom import setMatchDistance


class BaseTestCase(unittest.TestCase):

    """A test case for setMatchDistance

    Use involves setting one class attribute:
    * MatchClass: match class, e.g. ReferenceMatch or SourceMatch

    This test is a bit messy because it exercises two templates of MatchClass

    The test creates source and reference object catalogs that intentionally have
    some separation in on-sky coordinates. The reference catalog is set using
    a uniform grid of pixel positions and a simple WCS to compute on-sky coordinates.
    The source catalog is created by using a distorted version of the same grid
    of pixel positions, which is converted to on-sky coordinates using the same WCS.
    """
    MatchClass = None

    def setUp(self):
        crval = afwGeom.SpherePoint(44, 45, afwGeom.degrees)
        crpix = afwGeom.PointD(0, 0)

        scale = 1 * afwGeom.arcseconds
        self.tanWcs = afwGeom.makeSkyWcs(crpix=crpix, crval=crval,
                                         cdMatrix=afwGeom.makeCdMatrix(scale=scale))

        S = 300
        N = 5

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
        self.refCoordKey = afwTable.CoordKey(refSchema["coord"])
        self.srcCoordKey = afwTable.CoordKey(srcSchema["coord"])
        self.srcCentroidKey = afwTable.Point2DKey(srcSchema["slot_Centroid"])
        self.sourceCat = afwTable.SourceCatalog(srcSchema)
        self.origSourceCat = afwTable.SourceCatalog(srcSchema)  # undistorted copy
        self.matches = []

        for i in np.linspace(0., S, N):
            for j in np.linspace(0., S, N):
                src = self.sourceCat.addNew()
                refObj = self.refCat.addNew()

                src.set(self.srcCentroidKey, afwGeom.Point2D(i, j))

                c = self.tanWcs.pixelToSky(afwGeom.Point2D(i, j))
                refObj.setCoord(c)

                self.matches.append(self.MatchClass(refObj, src, 0.0))

    def tearDown(self):
        del self.refCat
        del self.origSourceCat
        del self.sourceCat
        del self.matches
        del self.tanWcs

    def doTest(self, name, func):
        """Apply func(x, y) to each source in self.sourceCat, then set coord, compute and check dist
        """
        for refObj, src, d in self.matches:
            origPos = src.get(self.srcCentroidKey)
            x, y = func(*origPos)
            distortedPos = afwGeom.Point2D(*func(*origPos))
            src.set(self.srcCentroidKey, distortedPos)
            src.set(self.srcCoordKey, self.tanWcs.pixelToSky(distortedPos))

        setMatchDistance(self.matches)
        maxDistErr = afwGeom.Angle(0)
        for refObj, source, distRad in self.matches:
            sourceCoord = source.get(self.srcCoordKey)
            refCoord = refObj.get(self.refCoordKey)
            predDist = sourceCoord.separation(refCoord)
            distErr = abs(predDist - distRad*afwGeom.radians)
            maxDistErr = max(distErr, maxDistErr)

        self.assertLess(maxDistErr.asArcseconds(), 1e-7)


class SideLoadTestCases(object):

    """Base class implementations of testing methods.

    Explicitly does not inherit from unittest.TestCase"""

    def testTrivial(self):
        """Add no distortion"""
        self.doTest("testTrivial", lambda x, y: (x, y))

    def testQuadraticX(self):
        """Add quadratic distortion in x"""
        self.doTest("testQuadraticX", lambda x, y: (x + 1e-4*x**2, y))

    def testRadial(self):
        """Add radial distortion"""
        radialTransform = afwGeom.makeRadialTransform([0, 1.02, 1e-6])

        def radialDistortion(x, y):
            x, y = radialTransform.applyForward(afwGeom.Point2D(x, y))
            return (x, y)
        self.doTest("testRadial", radialDistortion)

# The test classes inherit from two base classes and differ in the match
# class being used.


class SetMatchDistanceTestCaseReferenceMatch(BaseTestCase, SideLoadTestCases):
    MatchClass = afwTable.ReferenceMatch


class SetMatchDistanceTestCaseSourceMatch(BaseTestCase, SideLoadTestCases):
    MatchClass = afwTable.SourceMatch


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
