#!/usr/bin/env python

import os
import sys
import unittest

import lsst.afw.image as afwImage
import lsst.utils.tests as utilsTests
from lsst.afw.coord import DEGREES



class SipTransformationTest(unittest.TestCase):

    def setUp(self):
        basefn = os.path.join(os.path.dirname(__file__), 'imgCharSources-v85501867-R01-S00')
        wcsfn = basefn + '.header' #'.wcs'
        print 'WCS filename:', wcsfn
        metaData = afwImage.readMetadata(wcsfn)
        self.wcs = afwImage.makeWcs(metaData)
        print 'WCS object:', self.wcs

        # Using Astrometry.net:
        # > wcs-rd2xy -w $V.wcs -r 1.5 -d 3.3
        # RA,Dec (1.500000, 3.300000) -> pixel (3711.012884, 3134.440250)
        # > wcs-xy2rd -w $V.wcs -x 3711.012884 -y 3134.440250
        # Pixel (3711.012884, 3134.440250) -> RA,Dec (1.500000, 3.300000)

        # If only TAN is used:
        # > wcs-rd2xy -t -w $V.wcs -r 1.5 -d 3.3
        # RA,Dec (1.500000, 3.300000) -> pixel (3711.902259, 3134.017979)
        # > wcs-xy2rd -w $V.wcs -x 3711.012884 -y 3134.440250 -t
        # Pixel (3711.012884, 3134.440250) -> RA,Dec (1.500016, 3.300052)

        # These are 1-indexed pixels, hence the '-1's below.

        self.rdxy = [
            (1.500000, 3.300000, 3711.012884 - 1, 3134.440250 - 1),
            ]

    # UGH, the coord interface is nasty.
    def pixelToRaDec(self, wcs, xx, yy):
            rd = wcs.pixelToSky(xx, yy)
            print 'rd is', rd
            #rr = rd.getRa(DEGREES)
            #dd = rd.getDec(DEGREES)
            rr = rd.getLongitude(DEGREES)
            dd = rd.getLatitude(DEGREES)
            return (rr, dd)

    def testRoundTrip(self):
        for (ra,dec,x,y) in self.rdxy:
            xx,yy = self.wcs.skyToPixel(ra, dec)
            rr,dd = self.pixelToRaDec(self.wcs, xx, yy)
            self.assertAlmostEqual(rr, ra, 5)
            self.assertAlmostEqual(dd, dec, 5)

            ra,dec = self.pixelToRaDec(self.wcs, x, y)
            xx,yy = self.wcs.skyToPixel(ra, dec)
            self.assertAlmostEqual(x, xx, 5)
            self.assertAlmostEqual(y, yy, 5)


    def testSip1(self):
        for (ra, dec, x, y) in self.rdxy:
            xx,yy = self.wcs.skyToPixel(ra, dec)
            self.assertAlmostEqual(x, xx, 5)
            self.assertAlmostEqual(y, yy, 5)
            rr,dd = self.pixelToRaDec(self.wcs, x, y)
            self.assertAlmostEqual(r, rr, 5)
            self.assertAlmostEqual(d, dd, 5)


# UGH boilerplate.
def suite():
    """Returns a suite containing all the test cases in this module."""
    suites = []
    suites += unittest.makeSuite(SipTransformationTest)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
