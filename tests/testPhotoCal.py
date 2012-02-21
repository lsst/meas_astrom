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

import re
import os
import sys
import glob
import math
import unittest

import matplotlib
matplotlib.use('Agg')
import pylab as plt

import numpy as np

import eups
import lsst.meas.astrom            as measAstrom
import lsst.meas.algorithms.utils  as measAlgUtil
import lsst.afw.detection          as afwDet
import lsst.afw.math               as afwMath
import lsst.afw.image              as afwImg
import lsst.utils.tests            as utilsTests
import lsst.pex.policy             as pexPolicy
from lsst.pex.logging import Log
import lsst.meas.photocal          as photocal

from lsst.pex.exceptions import LsstCppException

import sourceSetIO                 as ssi


class PhotoCalTest(unittest.TestCase):

    def setUp(self):
        self.conf = measAstrom.MeasAstromConfig()
        
        #Load sample input from disk
        mypath = eups.productDir("meas_astrom")
        path = os.path.join(mypath, "examples")
        self.srcSet = ssi.read(os.path.join(path, "v695833-e0-c000.xy.txt"))
        for s in self.srcSet:
            s.setApFlux(s.getPsfFlux())
        
        # The .xy.txt file has sources in the range ~ [0,2000],[0,4500]
        self.imageSize = (2048, 4612) # approximate
        self.exposure = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci"))

        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))

    def tearDown(self):
        del self.srcSet
        del self.conf
        del self.exposure

    def getAstrometrySolution(self, loglvl = Log.INFO):
        astrom = measAstrom.Astrometry(self.conf, logLevel=loglvl)
        print 'Calling determineWcs...'
        res = astrom.determineWcs(self.srcSet, self.exposure, imageSize=self.imageSize)
        return res

    def testGetSourceMatch(self):
        res = self.getAstrometrySolution(loglvl=Log.DEBUG)
        print 'Result:', res
        M = res.getMatches()
        print 'Matches:', M
        print 'N matches:', len(M)
        assert(len(M) > 50)

        refflux = np.array([m.first.getPsfFlux() for m in M])
        srcflux = np.array([m.second.getPsfFlux() for m in M])
        I = np.logical_and(refflux > 0, srcflux > 0)
        refflux = refflux[I]
        srcflux = srcflux[I]
        refmag = -2.5 * np.log10(refflux)
        srcmag = -2.5 * np.log10(srcflux)
        plt.clf()
        plt.plot(srcmag, refmag, 'r.')
        plt.xlabel('src mag')
        plt.ylabel('ref mag')
        plt.savefig('mags1.png')

        plt.clf()
        plt.plot(refmag, srcmag - refmag, 'r.')
        plt.xlabel('ref mag')
        plt.ylabel('src mag - ref mag')
        plt.savefig('mags2.png')



    def test1(self):
        res = self.getAstrometrySolution()
        matches = res.getMatches()
           
        print 'Test1'
        #for m in matches:
        #    print '  ', m

        logLevel = Log.DEBUG
        log = Log(Log.getDefaultLog(),
                  'meas.astrom',
                  logLevel)

        pCal = photocal.calcPhotoCal(matches, log=log)
        print pCal

        diff=[]
        for m in matches:
            catFlux = m[0].getPsfFlux()     #Catalogue flux
            catMag = -2.5*np.log10(catFlux) #Cat mag
            instFlux = m[1].getPsfFlux()    #Instrumental Flux
            mag = pCal.getMag(instFlux)     #Instrumental mag
            diff.append(mag-catMag)

        #A very loose test, but the input data has a lot of scatter

        diff = np.array(diff)
        self.assertAlmostEqual(np.mean(diff), 0, 0)
        
        
    def test2(self):
        """Check that negative fluxes dealt with properly"""
        
        astrom = measAstrom.determineWcs(self.defaultPolicy, self.exposure,
                                         self.srcSet, forceImageSize=self.forceImageSize)
        matches = astrom.getMatches()
        
        matches[0].first.setPsfFlux(0)
        matches[1].first.setPsfFlux(-1)
        matches[2].second.setPsfFlux(0)
        
        pCal = photocal.calcPhotoCal(matches)
        print pCal


        diff=[]
        for m in matches:
            catFlux = m[0].getPsfFlux()     #Catalogue flux
            instFlux = m[1].getPsfFlux()    #Instrumental Flux

            if catFlux > 0 and instFlux > 0:
                catMag = -2.5*np.log10(catFlux) #Cat mag
                mag = pCal.getMag(instFlux)     #Instrumental mag
                diff.append(mag-catMag)


        #A very loose test, but the input data has a lot of scatter

        diff = np.array(diff)
        self.assertAlmostEqual(np.mean(diff), 0, 0)


    def testKnownZP(self):
        """Verify we recover a known zeropoint

        N.b. This is a very bad test as there's no scatter in the zeropoint measurements, but it does
        test this edge case of the photo calibration code
        """
        
        nS = 20
        magLo = 10.0
        magHi = 20.0
        dmag = (magHi - magLo)/nS

        # This test only works if zpCat == 0
        # - calcPhotoCal() assumes the getPsfFlux() values for catalog Sources
        #   were computed as f = 10**(-(mag-zp)/2.5), with zp=0, so we must do exactly that.
        zpCat = 0
        zpSrc = 20

        flags = measAlgUtil.getDetectionFlags()
        def fluxToMag(flux, zp=0):
            return zp - 2.5*math.log10(flux)
        def magToFlux(mag, zp=0):
            return 10**(-(mag-zp)/2.5)

        
        matchList = []
        for i in range(nS + 1):
            s1 = afwDet.Source()
            s2 = afwDet.Source()
            s1.setFlagForDetection(flags["BINNED1"] | flags["STAR"])
            s2.setFlagForDetection(flags["BINNED1"] | flags["STAR"])
            
            m1 = magLo + i*dmag
            f1 = magToFlux(m1, zp=zpCat) # the catalog mag

            # set the instrument flux for a different zeropoint
            f2 = magToFlux(m1, zp=zpSrc) 
            
            s1.setPsfFlux(f1)
            s2.setPsfFlux(f2)
            s1.setApFlux(f1)
            s2.setApFlux(f2)

            if False:
                print "flux1 =", f1, "flux2 =", f2
            matchList.append(afwDet.SourceMatch(s1, s2, 0.0))

        # do the cal
        pCal = photocal.calcPhotoCal(matchList, goodFlagValue=flags["BINNED1"])

        print "ZP_known = ", zpSrc, "ZP_cat", zpCat, "ZP = ", pCal.zeroMag, "Zflux = ", pCal.zeroFlux
        self.assertAlmostEqual(zpSrc, pCal.zeroMag)

        for m in matchList:
            s1, s2 = m.first, m.second
            mag1Known = fluxToMag(s1.getPsfFlux(), zpCat)  # catalog
            mag2Known = fluxToMag(s2.getPsfFlux(), zpSrc)  # inst

            # calibrate the fluxes and see if we get back the mags we put in.
            mag2 = pCal(s2.getPsfFlux())                    # inst calib
            
            print mag1Known, mag2Known, mag2
            self.assertAlmostEqual(mag2Known, mag2)
            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PhotoCalTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
