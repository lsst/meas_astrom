
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
import math
import os
import unittest

import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.utils.tests
import lsst.afw.image as afwImage
from lsst.meas.algorithms import LoadReferenceObjectsTask
import lsst.meas.astrom.sip.genDistortedImage as distort
import lsst.meas.astrom as measAstrom


class TestMatchOptimisticB(unittest.TestCase):

    def setUp(self):

        self.config = measAstrom.MatchPessimisticBTask.ConfigClass()
        self.MatchPessimisticB = measAstrom.MatchPessimisticBTask(
            config=self.config)

        metadata = dafBase.PropertySet()
        metadata.set("RADECSYS", "FK5")
        metadata.set("EQUINOX", 2000.0)
        metadata.set("CTYPE1", "RA---TAN")
        metadata.set("CTYPE2", "DEC--TAN")
        metadata.set("CUNIT1", "deg")
        metadata.set("CUNIT2", "deg")
        metadata.set("CRVAL1", 36.930640)
        metadata.set("CRVAL2", -4.939560)
        metadata.set("CRPIX1", 792.4)
        metadata.set("CRPIX2", 560.7)
        metadata.set("CD1_1", -5.17e-05)
        metadata.set("CD1_2", 0.0)
        metadata.set("CD2_2", 5.17e-05)
        metadata.set("CD2_1", 0.0)
        self.wcs = afwImage.makeWcs(metadata)
        self.distortedWcs = self.wcs

        self.filename = os.path.join(os.path.dirname(__file__), "cat.xy.fits")
        self.tolArcsec = .4
        self.tolPixel = .1

    def tearDown(self):
        del self.config
        del self.MatchPessimisticB
        del self.wcs
        del self.distortedWcs

    def testLinearXDistort(self):
        self.singleTestInstance(self.filename, distort.linearXDistort)

    def testLinearYDistort(self):
        self.singleTestInstance(self.filename, distort.linearYDistort)

    def testQuadraticDistort(self):
        self.singleTestInstance(self.filename, distort.quadraticDistort)

    def testLargeDistortion(self):
        # This transform is about as extreme as I can get:
        # using 0.0005 in the last value appears to produce numerical issues.

        # It produces a maximum deviation of 459 pixels, which should be
        # sufficient.
        pixelsToTanPixels = afwGeom.RadialXYTransform([0.0, 1.1, 0.0004])
        self.distortedWcs = afwImage.DistortedTanWcs(
            self.wcs, pixelsToTanPixels)


        def applyDistortion(src):
            out = src.table.copyRecord(src)
            out.set(out.table.getCentroidKey(),
                    pixelsToTanPixels.reverseTransform(src.getCentroid()))
            return out

        self.singleTestInstance(self.filename, applyDistortion)

    def singleTestInstance(self, filename, distortFunc, doPlot=False):
        sourceCat = self.loadSourceCatalog(self.filename)
        refCat = self.computePosRefCatalog(sourceCat)
        distortedCat = distort.distortList(sourceCat, distortFunc)

        if doPlot:
            import matplotlib.pyplot as plt

            undistorted = [self.wcs.skyToPixel(
                            self.distortedWcs.pixelToSky(ss.getCentroid()))
                           for ss in distortedCat]
            refs = [self.wcs.skyToPixel(ss.getCoord()) for ss in refCat]

            def plot(catalog, symbol):
                plt.plot([ss.getX() for ss in catalog],
                         [ss.getY() for ss in catalog], symbol)

            plot(distortedCat, 'b+')  # Distorted positions: blue +
            plot(undistorted, 'g+')  # Undistorted positions: green +
            plot(refs, 'rx')  # Reference catalog: red x
            # The green + should overlap with the red x, because that's how
            # MatchPessimisticB does it.

            plt.show()

        sourceCat = distortedCat

        matchRes = self.MatchPessimisticB.matchObjectsToSources(
            refCat=refCat,
            sourceCat=sourceCat,
            wcs=self.distortedWcs,
            refFluxField="r_flux",
        )
        matches = matchRes.matches
        if doPlot:
            measAstrom.plotAstrometry(matches=matches, refCat=refCat,
                                      sourceCat=sourceCat)
        if distortFunc == distort.quadraticDistort:
            # Quad distort finds 181 real matches for Pessimistic
            self.assertEqual(len(matches), 184)
        else:
            self.assertEqual(len(matches), 186)

        refCoordKey = afwTable.CoordKey(refCat.schema["coord"])
        srcCoordKey = afwTable.CoordKey(sourceCat.schema["coord"])
        refCentroidKey = afwTable.Point2DKey(refCat.getSchema()["centroid"])
        maxDistErr = afwGeom.Angle(0)

        for refObj, source, distRad in matches:
            sourceCoord = source.get(srcCoordKey)
            refCoord = refObj.get(refCoordKey)
            predDist = sourceCoord.angularSeparation(refCoord)
            distErr = abs(predDist - distRad*afwGeom.radians)
            maxDistErr = max(distErr, maxDistErr)

            # Current skip refObj 15 as it matches to source 12.
            # All other matches are correct.
            if refObj.getId() != source.getId() and refObj.getId() != 15:
                refCentroid = refObj.get(refCentroidKey)
                sourceCentroid = source.getCentroid()
                radius = math.hypot(*(refCentroid - sourceCentroid))
                self.fail(
                    "ID mismatch: %s at %s != %s at %s; error = %0.1f pix" %
                    (refObj.getId(), refCentroid, source.getId(),
                     sourceCentroid, radius))

        self.assertLess(maxDistErr.asArcseconds(), 1e-7)

    def computePosRefCatalog(self, sourceCat):
        """Generate a position reference catalog from a source catalog
        """
        minimalPosRefSchema = LoadReferenceObjectsTask.makeMinimalSchema(
            filterNameList=["r"],
            addFluxSigma=True,
        )
        refCat = afwTable.SimpleCatalog(minimalPosRefSchema)
        for source in sourceCat:
            refObj = refCat.addNew()
            refObj.setCoord(source.getCoord())
            refObj.set("centroid_x", source.getX())
            refObj.set("centroid_y", source.getY())
            refObj.set("hasCentroid", True)
            refObj.set("r_flux", source.get("slot_ApFlux_flux"))
            refObj.set("r_fluxSigma", source.get("slot_ApFlux_fluxSigma"))
            refObj.setId(source.getId())
        return refCat

    def loadSourceCatalog(self, filename):
        """Load a list of xy points from a file, set coord, and return a
        SourceSet of points

        """
        sourceCat = afwTable.SourceCatalog.readFits(filename)
        aliasMap = sourceCat.schema.getAliasMap()
        aliasMap.set("slot_ApFlux", "base_PsfFlux")
        fluxKey = sourceCat.schema["slot_ApFlux_flux"].asKey()
        fluxSigmaKey = sourceCat.schema["slot_ApFlux_fluxSigma"].asKey()

        # Source x,y positions are ~ (500,1500) x (500,1500)
        centroidKey = sourceCat.table.getCentroidKey()
        for src in sourceCat:
            adjCentroid = src.get(centroidKey) - afwGeom.Extent2D(500, 500)
            src.set(centroidKey, adjCentroid)
            src.set(fluxKey, 1000)
            src.set(fluxSigmaKey, 1)

        # Set catalog coord
        for src in sourceCat:
            src.updateCoord(self.wcs)
        return sourceCat


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":

    lsst.utils.tests.init()
    unittest.main()
