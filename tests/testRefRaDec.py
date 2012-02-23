#!/usr/bin/env python

# 2011-05-02, RHL reports:
'''
I cannot read the latest weekly run astrometry_net catalogues with the
latest meas_astrom, but I _can_ read them with the previous version
(20851:21214).  This is very odd, as I think that the code ran
correctly.

By "cannot read" I mean that all the XAstrom and YAstrom values in

    X = solver.getCatalogue(ra, dec, radius, filterName, idName, anid)
    self.ref = X.refsources

self.ref are 0.0
'''
import os
import unittest
import math

import eups
import sourceSetIO                 as ssi
import lsst.meas.astrom            as measAstrom
import lsst.afw.image              as afwImg
import lsst.afw.coord              as afwCoord
import lsst.utils.tests            as utilsTests
import lsst.pex.policy             as pexPolicy
from lsst.pex.logging import Log
import lsst.afw.geom  as afwGeom
                
class XYAstromTest(unittest.TestCase):

    def setUp(self):
        # Set up local astrometry_net_data
        mypath = eups.productDir("meas_astrom")
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))
        conf = measAstrom.MeasAstromConfig()
        loglvl = Log.INFO
        #loglvl = Log.DEBUG
        self.astrom = measAstrom.Astrometry(conf, logLevel=loglvl)

    def tearDown(self):
        del self.astrom

    def testXYAstrom(self):
        ra,dec = 215.505616846, 53.1874283624
        print 'ra,dec', ra,dec
        radius = 300.
        filterName = 'r'

        ref = self.astrom.getReferenceSources(ra * afwGeom.degrees, dec * afwGeom.degrees,
                                              radius * afwGeom.arcseconds, filterName)
        print 'Got', len(ref), 'ref sources'
        # Number of sources within the search radius.
        self.assertEqual(len(ref), 245)
        allzero = True
        for i in range(len(ref)):
            if i < 10:
                print ref[i]
                print 'r,d', i, ':', ref[i].getRa().asDegrees(), ref[i].getDec().asDegrees()
            ra = ref[i].getRa().asRadians()
            dec = ref[i].getDec().asRadians()
            if ra != 0 or dec != 0:
                allzero = False
            self.assertTrue(ra >= 0.)
            self.assertTrue(ra < 2.*math.pi)
            self.assertTrue(dec >= -math.pi)
            self.assertTrue(dec <=  math.pi)
            self.assertEqual(ra,  ref[i].getRaAstrom().asRadians())
            self.assertEqual(dec, ref[i].getDecAstrom().asRadians())
            self.assertEqual(ra,  ref[i].getRaFlux().asRadians())
            self.assertEqual(dec, ref[i].getDecFlux().asRadians())
            self.assertEqual(ra,  ref[i].getRaPeak().asRadians())
            self.assertEqual(dec, ref[i].getDecPeak().asRadians())
            self.assertEqual(ra,  ref[i].getRaObject().asRadians())
            self.assertEqual(dec, ref[i].getDecObject().asRadians())
        self.assertFalse(allzero)

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(XYAstromTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
