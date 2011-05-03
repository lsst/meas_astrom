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
                
class XYAstromTest(unittest.TestCase):

    def setUp(self):
        self.defaultPolicy = pexPolicy.Policy.createPolicy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>     
            matchThreshold: 22
            """
        ))
        # Set up local astrometry_net_data
        mypath = eups.productDir("meas_astrom")
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
		# Work around lame scons bug (doesn't pass HOME)
        os.environ['HOME'] = 'iswheretheheartis'
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))

    def tearDown(self):
        del self.defaultPolicy

    def testXYAstrom(self):
        log = Log.getDefaultLog()
        solver = measAstrom.createSolver(self.defaultPolicy, log)
        ra,dec = 215.505616846, 53.1874283624
        print 'ra,dec', ra,dec
        radius = 300.
        filterName = 'mag'
        idName = 'id'
        anid = 2033

        X = solver.getCatalogue(ra, dec, radius, filterName, idName, anid)
        ref = X.refsources
        print 'Got', len(ref), 'ref sources'
        # Number of sources within the search radius.
        self.assertEqual(len(ref), 143)
        allzero = True
        for i in range(10):
            print ref[i]
            print 'r,d', i, ':', ref[i].getRa(), ref[i].getDec()
            #print 'r,d', i, ':', ref[i].getRaAstrom(), ref[i].getDecAstrom()
            ra = ref[i].getRa()
            dec = ref[i].getDec()
            if ra != 0 or dec != 0:
                allzero = False
            self.assertTrue(ra >= 0.)
            self.assertTrue(ra < 2.*math.pi)
            self.assertTrue(dec >= -math.pi)
            self.assertTrue(dec <=  math.pi)
            self.assertEqual(ra,  ref[i].getRaAstrom())
            self.assertEqual(dec, ref[i].getDecAstrom())
            self.assertEqual(ra,  ref[i].getRaFlux())
            self.assertEqual(dec, ref[i].getDecFlux())
            self.assertEqual(ra,  ref[i].getRaPeak())
            self.assertEqual(dec, ref[i].getDecPeak())
            self.assertEqual(ra,  ref[i].getRaObject())
            self.assertEqual(dec, ref[i].getDecObject())
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
