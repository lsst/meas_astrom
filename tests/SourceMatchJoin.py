#!/usr/bin/env python

import os, sys
from math import *
import unittest
import eups
import random
import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils
import lsst.meas.astrom as measAstrom
from lsst.pex.logging import Log
#from lsst.afw.geom import Angle
import lsst.afw.geom as afwGeom
import lsst.afw.image              as afwImg

try:
    type(verbose)
except NameError:
    display = False
    verbose = 0

def roundTripSourceMatch(storagetype, filename, matchlist):
    pol = policy.Policy()
    additionalData = dafBase.PropertySet()

    loc = dafPersist.LogicalLocation(filename)
    persistence = dafPersist.Persistence.getPersistence(pol)
    storageList = dafPersist.StorageList()
    storage = persistence.getPersistStorage(storagetype, loc)
    storageList.append(storage)
    persistence.persist(matchlist, storageList, additionalData)

    storageList2 = dafPersist.StorageList()
    storage2 = persistence.getRetrieveStorage(storagetype, loc)
    storageList2.append(storage2)
    matchlistptr = persistence.unsafeRetrieve("PersistableSourceMatchVector", storageList2, additionalData)
    matchlist2 = afwDet.PersistableSourceMatchVector.swigConvert(matchlistptr)

    return matchlist2

class matchlistTestCase(unittest.TestCase):
    def setUp(self):
        # Load sample input from disk
        mypath = eups.productDir("meas_astrom")
        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        print 'Setting up astrometry_net_data:', datapath
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Couldn't set up local photocal version of astrometry_net_data (from path: %s): %s" % (datapath, reason))

        path = os.path.join(mypath, "examples")
        self.srcSet = afwTable.SourceCatalog.readFits(os.path.join(path, "v695833-e0-c000.xy.fits"))
        self.imageSize = (2048, 4612) # approximate
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci"))

        conf = measAstrom.MeasAstromConfig()
        loglvl = Log.INFO
        #loglvl = Log.DEBUG
        self.astrom = measAstrom.Astrometry(conf, logLevel=loglvl)

    def tearDown(self):
        del self.srcSet
        del self.imageSize
        del self.exposure
        del self.astrom
        
    def getAstrometrySolution(self):
        return self.astrom.determineWcs(self.srcSet, self.exposure, imageSize=self.imageSize)

    def testJoin(self):
        res = self.getAstrometrySolution()

        matches = res.getMatches()
        print 'Matches:', matches
        matchmeta = res.getMatchMetadata()
        print 'match meta', matchmeta

        smv = matches
        psmv = afwDet.PersistableSourceMatchVector(matches)
        psmv.setSourceMatchMetadata(matchmeta)

        psmv2 = roundTripSourceMatch('FitsStorage', 'tests/data/matchlist.fits', psmv)
        smv2 = psmv2.getSourceMatches()
        extra2 = psmv2.getSourceMatchMetadata()

        print 'Got SMV:', smv2
        print 'Got metadata:', extra2
        print extra2.toString()

        self.astrom.joinMatchListWithCatalog(smv2, extra2)

        self.assertEqual(len(smv2), len(smv))
        for i in xrange(len(smv)):
            self.assertEqual(smv2[i].first.getSourceId(), smv[i].first.getSourceId())
            self.assertEqual(smv2[i].second.getSourceId(), smv[i].second.getSourceId())
            self.assertEqual(smv2[i].first.getRa().asDegrees(), smv[i].first.getRa().asDegrees())
            self.assertEqual(smv2[i].first.getDec().asDegrees(), smv[i].first.getDec().asDegrees())
            self.assertEqual(smv2[i].first.getPsfFlux(), smv[i].first.getPsfFlux())

        return

            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= silly boilerplate -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(matchlistTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
