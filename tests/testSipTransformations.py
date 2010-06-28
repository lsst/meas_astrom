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

        sipfn = basefn + '.header' #'.wcs'
        tanfn = basefn + '.tanheader' #'.wcs'
        print 'TAN WCS filename:', tanfn
        print 'SIP WCS filename:', sipfn

        metaData = afwImage.readMetadata(tanfn)
        self.tan = afwImage.makeWcs(metaData)
        metaData = afwImage.readMetadata(sipfn)
        self.sip = afwImage.makeWcs(metaData)

        # Using Astrometry.net on SIP solution:

        '''
        wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.header -x 2168.54521667 -y 2020.40323873

        wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.header -x 2168.54521667 -y 2020.40323873

        wcs-rd2xy -w tests/imgCharSources-v85501867-R01-S00.header -r 1.5 -d 3.3
        '''
        
        self.sip_rdxy = [
            (1.42667846826, 3.37583321746, 2167.54521667 - 1, 2020.40323873 - 1),

            (1.4266863759, 3.3757783481,  2168.5452166700 - 1, 2020.4032387300 - 1),

            (1.5000000000, 3.3000000000, 3711.0128841126 - 1, 3134.4402504066 - 1),
            ]


        # If only TAN is used:

        # > wcs-rd2xy -w tests/imgCharSources-v85501867-R01-S00.tanheader -r 1.5 -d 3.3
        # RA,Dec (1.5000000000, 3.3000000000) -> pixel (3711.9022585704, 3134.0179793251)
        #  (wcslib gives same)

        # > wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.tanheader -x 3711.9022585704 -y 3134.0179793251
        # Pixel (3711.9022585704, 3134.0179793251) -> RA,Dec (1.5000000000, 3.3000000000)


        # These are 1-indexed pixels, hence the '-1's below.

        # crpix0 2167.54521667
        # crpix1 2020.40323873
        # crval0 1.42667846826
        # crval1 3.37583321746

        #  (crpix.x, crpix.y)
        # > wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.tanheader -x 2167.54521667 -y 2020.40323873
        # Pixel (2167.545217, 2020.403239) -> RA,Dec (1.426678, 3.375833)

        #  (crpix.x + 1, crpix.y)
        # > wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.tanheader -x 2168.54521667 -y 2020.40323873
        # Pixel (2168.5452166700, 2020.4032387300) -> RA,Dec (1.4266863759, 3.3757783481)

        #  (crpix.x + 1, crpix.y) with WCSLib:
        # wcs-xy2rd -L -w tests/imgCharSources-v85501867-R01-S00.tanheader -x 2168.54521667 -y 2020.40323873
        # Pixel (2168.5452166700, 2020.4032387300) -> RA,Dec (1.4266863759, 3.3757783481)

        # inverse:
        # > wcs-rd2xy -L -w tests/imgCharSources-v85501867-R01-S00.tanheader -r 1.4266863759 -d 3.3757783481
        # RA,Dec (1.4266863759, 3.3757783481) -> pixel (2168.5452162485, 2020.4032394061)


        # This is the SIP position of 1.5,3.3:
        # > wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.tanheader -x 3711.0128841126 -y 3134.4402504066
        # Pixel (3711.0128841126, 3134.4402504066) -> RA,Dec (1.5000161440, 3.3000521757)

        # --> good to about 5 decimal digits in pixels.

        self.tan_rdxy = [
            (1.42667846826, 3.37583321746, 2167.54521667 - 1,   2020.40323873 - 1),
            (1.4266863759,  3.3757783481,  2168.5452166700 - 1, 2020.4032387300 - 1),
            (1.5000000000, 3.3000000000,  3711.9022585704 - 1, 3134.0179793251 - 1),
            (1.5000161440, 3.3000521757, 3711.0128841126 - 1, 3134.4402504066 - 1),
            ]

    # UGH, the coord interface is nasty.
    def pixelToRaDec(self, wcs, xx, yy):
            rd = wcs.pixelToSky(xx, yy)
            #print 'rd is', rd
            #rr = rd.getRa(DEGREES)
            #dd = rd.getDec(DEGREES)
            rr = rd.getLongitude(DEGREES)
            dd = rd.getLatitude(DEGREES)
            return (rr, dd)

    def roundTrip(self, wcs, rdxy):
        for (ra,dec,x,y) in rdxy:
            xx,yy = wcs.skyToPixel(ra, dec)
            rr,dd = self.pixelToRaDec(wcs, xx, yy)
            print
            print 'RA,Dec %-14.10g, %-14.10g --> pixel %g, %g -->' % (ra, dec, xx, yy)
            print 'RA,Dec %-14.10g, %-14.10g' % (rr, dd)

            #print 'pixels are', type(xx), type(yy)
            #print 'ra,dec are', type(rr), type(rr)
            
            self.assertAlmostEqual(rr, ra, 5)
            self.assertAlmostEqual(dd, dec, 5)

            ra,dec = self.pixelToRaDec(wcs, x, y)
            xx,yy = wcs.skyToPixel(ra, dec)
            print
            print 'Pixel %-14.10g, %-14.10g --> RA,Dec %-14.10g, %-14.10g -->' %  (x, y, ra, dec)
            print 'Pixel %-14.10g, %-14.10g' % (xx, yy)
            self.assertAlmostEqual(x, xx, 5)
            self.assertAlmostEqual(y, yy, 5)
        

    def testRoundTripTAN(self):
        print
        print 'Round trip TAN'
        self.roundTrip(self.tan, self.tan_rdxy)

    def testRoundTripSIP(self):
        print
        print 'Round trip SIP'
        self.roundTrip(self.sip, self.sip_rdxy)

    def againstReality(self, wcs, rdxy):
        for (ra, dec, x, y) in rdxy:
            xx,yy = wcs.skyToPixel(ra, dec)
            print 'RA,Dec %-14.10g, %-14.10g --> x,y %-14.10g, %-14.10g' % (ra, dec, xx, yy)
            print '  Expected:                                   %-14.10g, %-14.10g' % (x,y)
            self.assertAlmostEqual(x, xx, 5)
            self.assertAlmostEqual(y, yy, 5)
            rr,dd = self.pixelToRaDec(wcs, x, y)
            self.assertAlmostEqual(ra, rr, 5)
            self.assertAlmostEqual(dec, dd, 5)

    def testTan1(self):
        print
        print 'TAN against reality:'
        self.againstReality(self.tan, self.tan_rdxy)

    def testSip1(self):
        print
        print 'SIP against reality:'
        self.againstReality(self.sip, self.sip_rdxy)



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
