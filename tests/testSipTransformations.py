#!/usr/bin/env python

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


import os
import sys
import unittest

import lsst.afw.image as afwImage
import lsst.utils.tests as utilsTests
import lsst.afw.geom as afwGeom



class SipTransformationTest(unittest.TestCase):

    def tearDown(self):
        import lsst.meas.astrom.astrometry_net as an
        an.finalize()
    
    def setUp(self):
        basefn = os.path.join(os.path.dirname(__file__), 'imgCharSources-v85501867-R01-S00')

        sipfn = basefn + '.sipheader'
        tanfn = basefn + '.tanheader'
        print 'TAN WCS filename:', tanfn
        print 'SIP WCS filename:', sipfn

        metaData = afwImage.readMetadata(tanfn)
        self.tan = afwImage.makeWcs(metaData)
        metaData = afwImage.readMetadata(sipfn)
        self.sip = afwImage.makeWcs(metaData)

        #print 'TAN is', self.tan
        #print 'SIP is', self.sip

        # Using Astrometry.net on SIP solution:

        '''
        wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.sipheader -x 2168.54521667 -y 2020.40323873

        wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.sipheader -x 2168.54521667 -y 2020.40323873

        wcs-rd2xy -w tests/imgCharSources-v85501867-R01-S00.sipheader -r 1.5 -d 3.3

        wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.sipheader -x 0 -y 0

        Roundtrip:
        > wcs-rd2xy -w tests/imgCharSources-v85501867-R01-S00.sipheader -r 1.5 -d 3.3
        RA,Dec (1.5000000000, 3.3000000000) -> pixel (3711.0128841126, 3134.4402504066)
        > wcs-xy2rd -w tests/imgCharSources-v85501867-R01-S00.sipheader -x 3711.0128841126 -y 3134.4402504066
        Pixel (3711.0128841126, 3134.4402504066) -> RA,Dec (1.5000000064, 3.3000000037)
        > wcs-rd2xy -w tests/imgCharSources-v85501867-R01-S00.sipheader -r 1.5000000064 -d 3.3000000037
        RA,Dec (1.5000000064, 3.3000000037) -> pixel (3711.0128347543, 3134.4403740779)

        --> good to about 8 digits in RA,Dec in degrees.

        (3711.0128841126, 3134.4402504066)
        (3711.0128347543, 3134.4403740779)

        --> and about 3-4 digits in pixels.
        '''

        ''' WCSTools 3.8.0:
        > sky2xy -v -j tests/imgCharSources-v85501867-R01-S00.wcs 1.5 3.3
        1.5 3.3 J2000 -> 3711.013 3134.440

        > xy2sky -v -j -d tests/imgCharSources-v85501867-R01-S00.wcs 3711.013 3134.440
        XY2SKY WCSTools 3.8.0, 12 November 2009, Doug Mink (dmink@cfa.harvard.edu)
        Print sky coordinates from tests/imgCharSources-v85501867-R01-S00.wcs image coordinates
        RA           Dec       Sys          X        Y
        1.50000   3.30000 J2000  <- 3711.013 3134.440

        '''

        
        self.sip_rdxy = [(r * afwGeom.degrees, d * afwGeom.degrees, x-1,y-1) for (r,d,x,y) in [
            (1.42667846826, 3.37583321746, 2167.54521667, 2020.40323873),
            (1.4266863759, 3.3757783481,  2168.5452166700, 2020.4032387300),
            (1.5000000000, 3.3000000000, 3711.0128841126, 3134.4402504066),
            (1.2986490876, 3.4785952816, 0.0000000000, 0.0000000000),
            ]]


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

        self.tan_rdxy = [(r * afwGeom.degrees, d * afwGeom.degrees, x-1, y-1) for (r,d,x,y) in [
            (1.42667846826, 3.37583321746, 2167.54521667,   2020.40323873),
            (1.4266863759,  3.3757783481,  2168.5452166700, 2020.4032387300),
            (1.5000000000, 3.3000000000,  3711.9022585704, 3134.0179793251),
            (1.5000161440, 3.3000521757, 3711.0128841126, 3134.4402504066),
            (1.2986453989, 3.4785959230, 0.0000000000, 0.0000000000),
            ]]

    # UGH, the coord interface is nasty.
    def pixelToRaDec(self, wcs, xx, yy):
        rd = wcs.pixelToSky(xx, yy)
        rr = rd.getLongitude().asDegrees()
        dd = rd.getLatitude().asDegrees()
        return (rr, dd)

    def roundTrip(self, wcs, rdxy):
        for (ra,dec,x,y) in rdxy:
            xx,yy = wcs.skyToPixel(ra, dec)
            rr,dd = self.pixelToRaDec(wcs, xx, yy)
            print
            print 'RA,Dec %-14.12g, %-14.12g --> pixel %g, %g -->' % (ra, dec, xx, yy)
            print 'RA,Dec %-14.12g, %-14.12g' % (rr, dd)

            #print 'pixels are', type(xx), type(yy)
            #print 'ra,dec are', type(rr), type(rr)
            
            self.assertAlmostEqual(rr, ra.asDegrees(), 5)
            self.assertAlmostEqual(dd, dec.asDegrees(), 5)

            ra,dec = self.pixelToRaDec(wcs, x, y)
            xx,yy = wcs.skyToPixel(ra * afwGeom.degrees, dec * afwGeom.degrees)
            print
            print 'Pixel %-14.12g, %-14.12g --> RA,Dec %-14.12g, %-14.12g -->' %  (x, y, ra, dec)
            print 'Pixel %-14.12g, %-14.12g' % (xx, yy)
            self.assertAlmostEqual(x, xx, 3)
            self.assertAlmostEqual(y, yy, 3)
        

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
            print 'RA,Dec %-14.12g, %-14.12g --> x,y %-14.12g, %-14.12g' % (ra, dec, xx, yy)
            print '  Expected:                                   %-14.12g, %-14.12g' % (x,y)
            self.assertAlmostEqual(x, xx, 3)
            self.assertAlmostEqual(y, yy, 3)
            rr,dd = self.pixelToRaDec(wcs, x, y)
            self.assertAlmostEqual(ra.asDegrees(), rr, 5)
            self.assertAlmostEqual(dec.asDegrees(), dd, 5)
            print 'x,y %-14.12g, %-14.12g --> ra,dec %-14.12g, %-14.12g' % (x, y, ra, dec)
            print '  Expected:                                   %-14.12g, %-14.12g' % (rr,dd)

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
