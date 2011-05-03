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
        distanceForCatalogueMatchinArcsec: 1.0
        cleaningParameter: 3
        calculateSip: false
        numBrightStars: 75
        defaultFilterName: mag
        """
        ))
        
        #Load sample input from disk
        mypath = eups.productDir("meas_astrom")
        path = os.path.join(mypath, "examples")
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci"))
        self.srcSet = ssi.read(os.path.join(path, "v695833-e0-c000.xy.txt"))
        for s in self.srcSet:
            s.setApFlux(s.getPsfFlux())
        
        # The .xy.txt file has sources in the range ~ [0,2000],[0,4500], but
        # the exposure is only one amp -- 1024x1153.  Work around.
        print 'Exposure image size: %i x %i' % (self.exposure.getWidth(), self.exposure.getHeight())
        self.forceImageSize = (2048, 4612) # approximately; 2x4 x (1024 x 1153)
        print 'Forcing image size to %i x %i to match source list.' % (self.forceImageSize[0],
                                                                       self.forceImageSize[1])

        # Set up local astrometry_net_data
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
        del self.exposure
        del self.srcSet

    def testXYAstrom(self):
        log = Log.getDefaultLog()
        solver = measAstrom.createSolver(self.defaultPolicy, log)
        solver.setStarlist(self.srcSet)
        solver.setImageSize(*self.forceImageSize)
        key = 'pixelScaleUncertainty'
        policy = self.defaultPolicy
        if policy.exists(key):
            dscale = float(policy.get(key))
            solver.solve(self.exposure.getWcs(), dscale)
        else:
            solver.solve(self.exposure.getWcs())

        wcs = self.exposure.getWcs()
        radec = wcs.pixelToSky(0, 0)
        print 'ra,dec', radec
        ra,dec = (radec.toIcrs().getRa(afwCoord.DEGREES),
                  radec.toIcrs().getDec(afwCoord.DEGREES))
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
            print 'r,d', i, ':', ref[i].getRaAstrom(), ref[i].getDecAstrom()
            print 'x,y', i, ':', ref[i].getXAstrom(), ref[i].getYAstrom()
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
