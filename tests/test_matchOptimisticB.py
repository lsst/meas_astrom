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
import pickle

import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.utils.tests
import lsst.pex.exceptions as pexExcept
from lsst.meas.algorithms import LoadReferenceObjectsTask
import lsst.meas.astrom.sip.genDistortedImage as distort
import lsst.meas.astrom as measAstrom
import lsst.meas.astrom.matchOptimisticB as matchOptimisticB


class TestMatchOptimisticB(unittest.TestCase):

    def setUp(self):

        self.config = measAstrom.MatchOptimisticBTask.ConfigClass()
        self.matchOptimisticB = measAstrom.MatchOptimisticBTask(config=self.config)
        self.wcs = afwGeom.makeSkyWcs(crpix=lsst.geom.Point2D(791.4, 559.7),
                                      crval=lsst.geom.SpherePoint(36.930640, -4.939560, lsst.geom.degrees),
                                      cdMatrix=afwGeom.makeCdMatrix(scale=5.17e-5*lsst.geom.degrees))
        self.distortedWcs = self.wcs

        self.filename = os.path.join(os.path.dirname(__file__), "cat.xy.fits")
        self.tolArcsec = .4
        self.tolPixel = .1

    def tearDown(self):
        del self.config
        del self.matchOptimisticB
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
        # It produces a maximum deviation of 459 pixels, which should be sufficient.
        pixelsToTanPixels = afwGeom.makeRadialTransform([0.0, 1.1, 0.0004])
        self.distortedWcs = afwGeom.makeModifiedWcs(pixelTransform=pixelsToTanPixels,
                                                    wcs=self.wcs,
                                                    modifyActualPixels=False)

        def applyDistortion(src):
            out = src.table.copyRecord(src)
            out.set(out.table.getCentroidSlot().getMeasKey(),
                    pixelsToTanPixels.applyInverse(src.getCentroid()))
            return out

        self.singleTestInstance(self.filename, applyDistortion)

    def singleTestInstance(self, filename, distortFunc, doPlot=False):
        sourceCat = self.loadSourceCatalog(self.filename)
        refCat = self.computePosRefCatalog(sourceCat)

        # Apply source selector to sourceCat, using the astrometry config defaults
        tempConfig = measAstrom.AstrometryTask.ConfigClass()
        tempConfig.matcher.retarget(measAstrom.MatchOptimisticBTask)
        tempConfig.sourceSelector["matcher"].excludePixelFlags = False
        tempSolver = measAstrom.AstrometryTask(config=tempConfig, refObjLoader=None)
        sourceSelection = tempSolver.sourceSelector.run(sourceCat)

        distortedCat = distort.distortList(sourceSelection.sourceCat, distortFunc)

        if doPlot:
            import matplotlib.pyplot as plt
            undistorted = [self.wcs.skyToPixel(self.distortedWcs.pixelToSky(ss.getCentroid())) for
                           ss in distortedCat]
            refs = [self.wcs.skyToPixel(ss.getCoord()) for ss in refCat]

            def plot(catalog, symbol):
                plt.plot([ss.getX() for ss in catalog], [ss.getY() for ss in catalog], symbol)

            # plot(sourceCat, 'k+') # Original positions: black +
            plot(distortedCat, 'b+')  # Distorted positions: blue +
            plot(undistorted, 'g+')  # Undistorted positions: green +
            plot(refs, 'rx')  # Reference catalog: red x
            # The green + should overlap with the red x, because that's how matchOptimisticB does it.
            # The black + happens to overlap with those also, but that's beside the point.
            plt.show()

        sourceCat = distortedCat

        matchRes = self.matchOptimisticB.matchObjectsToSources(
            refCat=refCat,
            sourceCat=sourceCat,
            wcs=self.distortedWcs,
            sourceFluxField='slot_ApFlux_instFlux',
            refFluxField="r_flux",
        )
        matches = matchRes.matches
        if doPlot:
            measAstrom.plotAstrometry(matches=matches, refCat=refCat, sourceCat=sourceCat)
        self.assertEqual(len(matches), 183)

        refCoordKey = afwTable.CoordKey(refCat.schema["coord"])
        srcCoordKey = afwTable.CoordKey(sourceCat.schema["coord"])
        refCentroidKey = afwTable.Point2DKey(refCat.getSchema()["centroid"])
        maxDistErr = 0*lsst.geom.radians
        for refObj, source, distRad in matches:
            sourceCoord = source.get(srcCoordKey)
            refCoord = refObj.get(refCoordKey)
            predDist = sourceCoord.separation(refCoord)
            distErr = abs(predDist - distRad*lsst.geom.radians)
            maxDistErr = max(distErr, maxDistErr)

            if refObj.getId() != source.getId():
                refCentroid = refObj.get(refCentroidKey)
                sourceCentroid = source.getCentroid()
                radius = math.hypot(*(refCentroid - sourceCentroid))
                self.fail("ID mismatch: %s at %s != %s at %s; error = %0.1f pix" %
                          (refObj.getId(), refCentroid, source.getId(), sourceCentroid, radius))

        self.assertLess(maxDistErr.asArcseconds(), 1e-7)

    def computePosRefCatalog(self, sourceCat):
        """Generate a position reference catalog from a source catalog
        """
        minimalPosRefSchema = LoadReferenceObjectsTask.makeMinimalSchema(filterNameList=["r"],
                                                                         addCentroid=True)
        refCat = afwTable.SimpleCatalog(minimalPosRefSchema)
        for source in sourceCat:
            refObj = refCat.addNew()
            refObj.setCoord(source.getCoord())
            refObj.set("centroid_x", source.getX())
            refObj.set("centroid_y", source.getY())
            refObj.set("hasCentroid", True)
            refObj.set("r_flux", source.get("slot_ApFlux_instFlux"))
            refObj.set("r_fluxErr", source.get("slot_ApFlux_instFluxErr"))
            refObj.setId(source.getId())
        return refCat

    def loadSourceCatalog(self, filename):
        """Load a list of xy points from a file, set coord, and return a SourceSet of points
        """
        sourceCat = afwTable.SourceCatalog.readFits(filename)
        aliasMap = sourceCat.schema.getAliasMap()
        aliasMap.set("slot_ApFlux", "base_PsfFlux")
        instFluxKey = sourceCat.schema["slot_ApFlux_instFlux"].asKey()
        instFluxErrKey = sourceCat.schema["slot_ApFlux_instFluxErr"].asKey()

        # print("schema=", sourceCat.schema)

        # Source x,y positions are ~ (500,1500) x (500,1500)
        centroidKey = sourceCat.table.getCentroidSlot().getMeasKey()
        for src in sourceCat:
            adjCentroid = src.get(centroidKey) - lsst.geom.Extent2D(500, 500)
            src.set(centroidKey, adjCentroid)
            src.set(instFluxKey, 1000)
            src.set(instFluxErrKey, 1)

        # Set catalog coord
        for src in sourceCat:
            src.updateCoord(self.wcs)
        return sourceCat

    def testArgumentErrors(self):
        """Test argument sanity checking in matchOptimisticB
        """
        matchControl = matchOptimisticB.MatchOptimisticBControl()

        sourceCat = self.loadSourceCatalog(self.filename)
        emptySourceCat = afwTable.SourceCatalog(sourceCat.schema)

        refCat = self.computePosRefCatalog(sourceCat)
        emptyRefCat = afwTable.SimpleCatalog(refCat.schema)

        with self.assertRaises(pexExcept.InvalidParameterError):
            matchOptimisticB.matchOptimisticB(
                emptyRefCat,
                sourceCat,
                matchControl,
                self.wcs,
                0,
            )
        with self.assertRaises(pexExcept.InvalidParameterError):
            matchOptimisticB.matchOptimisticB(
                refCat,
                emptySourceCat,
                matchControl,
                self.wcs,
                0,
            )
        with self.assertRaises(pexExcept.InvalidParameterError):
            matchOptimisticB.matchOptimisticB(
                refCat,
                sourceCat,
                matchControl,
                self.wcs,
                len(refCat),
            )
        with self.assertRaises(pexExcept.InvalidParameterError):
            matchOptimisticB.matchOptimisticB(
                refCat,
                sourceCat,
                matchControl,
                self.wcs,
                -1,
            )

    def testConfigPickle(self):
        """Test that we can pickle the Config

        This is required for use in singleFrameDriver.
        See DM-18314.
        """
        config = pickle.loads(pickle.dumps(self.config))
        self.assertEqual(config, self.config)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":

    lsst.utils.tests.init()
    unittest.main()
