#!/usr/bin/env python
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

import eups
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.utils.tests as utilsTests
import lsst.afw.image as afwImage
from lsst.meas.algorithms import LoadReferenceObjectsTask
import lsst.meas.astrom.sip.genDistortedImage as distort
import lsst.meas.astrom as measAstrom

############################
# Set up local astrometry_net_data
meas_astrom_dir = eups.productDir("meas_astrom")
datapath = os.path.join(meas_astrom_dir, 'tests', 'astrometry_net_data', 'cfhttemplate')
eupsObj = eups.Eups(root=datapath)
ok, version, reason = eupsObj.setup('astrometry_net_data')
if not ok:
    raise ValueError("Can't find astrometry_net_data version cfhttemplate (from path: %s): %s" %
                     (datapath, reason))

class TestMatchOptimisticB(unittest.TestCase):
    def setUp(self):

        self.config = measAstrom.MatchOptimisticBTask.ConfigClass()
        self.matchOptimisticB = measAstrom.MatchOptimisticBTask(config = self.config)

        metadata = dafBase.PropertySet()
        metadata.set("RADECSYS", "FK5")
        metadata.set("EQUINOX", 2000.0)
        metadata.set("CTYPE1", "RA---TAN")
        metadata.set("CTYPE2", "DEC--TAN")
        metadata.set("CUNIT1", "deg")
        metadata.set("CUNIT2", "deg")
        metadata.set("CRVAL1",  36.930640)
        metadata.set("CRVAL2",  -4.939560)
        metadata.set("CRPIX1", 792.4)
        metadata.set("CRPIX2", 560.7)
        metadata.set("CD1_1", -5.17e-05)
        metadata.set("CD1_2",  0.0)
        metadata.set("CD2_2",  5.17e-05)
        metadata.set("CD2_1",  0.0)
        self.wcs = afwImage.makeWcs(metadata)

        path=eups.productDir("meas_astrom")
        self.filename=os.path.join(path, "tests", "cat.xy.fits")
        self.tolArcsec = .4 
        self.tolPixel = .1

    def tearDown(self):
        del self.config
        del self.matchOptimisticB
        del self.wcs
        
    def testLinearXDistort(self):
        self.singleTestInstance(self.filename, distort.linearXDistort)

    def testLinearYDistort(self):
        self.singleTestInstance(self.filename, distort.linearYDistort)

    def testQuadraticDistort(self):
        self.singleTestInstance(self.filename, distort.quadraticDistort)
    
    def singleTestInstance(self, filename, distortFunc, doPlot=False):
        sourceCat = self.loadSourceCatalog(self.filename)
        refCat = self.computePosRefCatalog(sourceCat)
        sourceCat = distort.distortList(sourceCat, distortFunc)
        matchRes = self.matchOptimisticB.matchObjectsToSources(
            refCat = refCat,
            sourceCat = sourceCat,
            wcs = self.wcs,
            refFluxField = "r_flux",
        )
        matches = matchRes.matches
        if doPlot:
            measAstrom.plotMatches(matches=matches, refCat=refCat, sourceCat=sourceCat)
        self.assertEqual(len(matches), 183)

        refCoordKey = refCat.schema["coord"].asKey()
        srcCoordKey = sourceCat.schema["coord"].asKey()
        maxDistErr = afwGeom.Angle(0)
        for refObj, source, distRad in matches:
            sourceCoord = source.get(srcCoordKey)
            refCoord = refObj.get(refCoordKey)
            predDist = sourceCoord.angularSeparation(refCoord)
            distErr = abs(predDist - distRad*afwGeom.radians)
            maxDistErr = max(distErr, maxDistErr)

            if refObj.getId() != source.getId():
                refCentroid = refObj.get("centroid")
                sourceCentroid = source.getCentroid()
                radius = math.hypot(*(refCentroid - sourceCentroid))
                self.fail("ID mismatch: %s at %s != %s at %s; error = %0.1f pix" %
                    (refObj.getId(), refCentroid, source.getId(), sourceCentroid, radius))

        self.assertLess(maxDistErr.asArcseconds(), 1e-7)

    def computePosRefCatalog(self, sourceCat):
        """Generate a position reference catalog from a source catalog
        """
        minimalPosRefSchema = LoadReferenceObjectsTask.makeMinimalSchema(
            filterNameList = ["r"],
            addFluxSigma = True,
        )
        refCat = afwTable.SimpleCatalog(minimalPosRefSchema)
        for source in sourceCat:
            refObj = refCat.addNew()
            refObj.set("coord", source.get("coord"))
            refObj.set("centroid", source.getCentroid())
            refObj.set("hasCentroid", True)
            refObj.set("r_flux", source.get("slot_ApFlux_flux"))
            refObj.set("r_fluxSigma", source.get("slot_ApFlux_fluxSigma"))
            refObj.setId(source.getId())
        return refCat

    def loadSourceCatalog(self, filename):
        """Load a list of xy points from a file, set coord, and return a SourceSet of points
        """
        sourceCat = afwTable.SourceCatalog.readFits(filename)
        aliasMap = sourceCat.schema.getAliasMap()
        aliasMap.set("slot_ApFlux_flux", "base_PsfFlux_flux")
        aliasMap.set("slot_ApFlux_fluxSigma", "base_PsfFlux_fluxSigma")

        # print("schema=", sourceCat.schema)

        # Source x,y positions are ~ (500,1500) x (500,1500)
        centroidKey = sourceCat.table.getCentroidKey()
        for src in sourceCat:
            adjCentroid = src.get(centroidKey) - afwGeom.Extent2D(500, 500)
            src.set(centroidKey, adjCentroid)
            
        # Set catalog coord
        for src in sourceCat:
            src.updateCoord(self.wcs)
        return sourceCat

        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(TestMatchOptimisticB)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)



 
if __name__ == "__main__":
    run(True)
