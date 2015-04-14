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
import unittest

try:                                    # used in plotPhotoCal
    import matplotlib
    matplotlib.use('Agg')
    import pylab as plt
except ImportError:
    plt = None

import numpy as np

import eups
import lsst.meas.astrom            as measAstrom
import lsst.afw.geom               as afwGeom
import lsst.afw.table              as afwTable
import lsst.afw.image              as afwImage
import lsst.utils.tests            as utilsTests
from lsst.pex.logging import Log
import lsst.meas.photocal          as photocal

class PhotoCalTest(unittest.TestCase):

    def setUp(self):
        self.conf = measAstrom.AstrometryConfig()
        
        # Load sample input from disk
        mypath = eups.productDir("meas_astrom")
        path = os.path.join(mypath, "examples")
        self.srcCat = afwTable.SourceCatalog.readFits(os.path.join(path, "v695833-e0-c000.xy.fits"))
        self.srcCat.table.defineApFlux("flux.psf")
        
        # The .xy.fits file has sources in the range ~ [0,2000],[0,4500]
        # which is bigger than the exposure
        self.bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(2048, 4612))
        self.exposure = afwImage.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci.fits"))

        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))

    def tearDown(self):
        del self.srcCat
        del self.conf
        del self.exposure
        import lsst.meas.astrom.astrometry_net as an
        an.finalize()
        
    def getAstrometrySolution(self, loglvl = Log.INFO):
        astromConfig = measAstrom.AstrometryTask.ConfigClass()
        astrom = measAstrom.AstrometryTask(config=astromConfig)
        # use solve instead of run because the exposure has the wrong bbox
        return astrom.solve(
            sourceCat = self.srcCat,
            bbox = self.bbox,
            initWcs = self.exposure.getWcs(),
            filterName = self.exposure.getFilter().getName(),
            calib = self.exposure.getCalib(),
        )

    def testGetSolution(self):
        res = self.getAstrometrySolution(loglvl=Log.DEBUG)
        self.assertTrue(res is not None)
        self.assertTrue(len(res.matches) > 50)

    def plotPhotoCal(self):
        res = self.getAstrometrySolution(loglvl=Log.DEBUG)
        print 'Result:', res
        M = res.matches
        #print 'Matches:', M
        print 'N matches:', len(M)
        assert(len(M) > 50)

        logLevel = Log.DEBUG
        log = Log(Log.getDefaultLog(), 'meas.astrom', logLevel)
        pcal = photocal.calcPhotoCal(M, log=log)
        print 'PhotoCal:', pcal
        zp = pcal.getMag(1.)
        print 'zeropoint:', zp

        refflux = np.array([m.first .getPsfFlux() for m in M])
        srcflux = np.array([m.second.getPsfFlux() for m in M])
        refferr = np.array([m.first .getPsfFluxErr() for m in M])
        srcferr = np.array([m.second.getPsfFluxErr() for m in M])

        refmag = -2.5 * np.log10(refflux)
        srcmag = -2.5 * np.log10(srcflux)
        refmagerr = refferr / refflux / np.log(10.)
        srcmagerr = srcferr / srcflux / np.log(10.)

        plt.clf()
        plt.errorbar(srcmag, refmag, xerr=srcmagerr, yerr=refmagerr, fmt='.',
                     ecolor='r', mec='r', mfc='r')
        ax = plt.axis()
        plt.plot([ax[0], ax[1]], [ax[0]+zp, ax[1]+zp], 'r-', lw=3, alpha=0.3)
        plt.axis(ax)
        plt.xlabel('src mag')
        plt.ylabel('ref mag')
        plt.savefig('mags1.png')

        def myerrbars(I, c):
            return plt.errorbar(XX[I], YY[I], xerr=XE[I], yerr=YE[I],
                                fmt='.', ecolor=c, mec=c, mfc=c)

        plt.clf()
        XX = refmag
        YY = srcmag-refmag
        XE = refmagerr
        YE = srcmagerr
        p1 = plt.errorbar(XX, YY, xerr=XE, yerr=YE, fmt='.',
                          ecolor='k', mec='k', mfc='k')
        LP = [p1[0]]
        LT = ['all']

        print 'Initial match list:', len(M)
        print 'Iflags:', len(pcal.Iflags)
        print 'Istar:', len(pcal.Istar)
        print 'Iflux:', len(pcal.Iflux)
        print 'Ibright:', len(pcal.Ibright)

        Ibadflags = np.ones(len(M), bool)
        Ibadflags[pcal.Iflags] = False
        p2 = myerrbars(Ibadflags, 'm')
        LP.append(p2[0])
        LT.append('bad flags')

        Inotstars = np.zeros(len(M), bool)
        Inotstars[pcal.Iflags] = True
        Inotstars[pcal.Istar] = False
        if sum(Inotstars):
            p2 = myerrbars(Inotstars, 'g')
            LP.append(p2[0])
            LT.append('not STARs')

        Ifaint = np.zeros(len(M), bool)
        Ifaint[pcal.Iflux] = True
        Ifaint[pcal.Ibright] = False
        if sum(Ifaint):
            p2 = myerrbars(Ifaint, 'y')
            LP.append(p2[0])
            LT.append('too faint')

        Igood = np.zeros(len(M), bool)
        Igood[pcal.Ibright] = True
        if sum(Igood):
            p2 = myerrbars(Igood, 'r')
            LP.append(p2[0])
            LT.append('good')

        plt.legend(LP, LT)
        #plt.plot(refmag, srcmag - refmag, 'r.')
        plt.axhline(-zp, color='r', lw=3, alpha=0.3)
        plt.xlabel('ref mag')
        plt.ylabel('src mag - ref mag')
        plt.savefig('mags2.png')



    def test1(self):
        res = self.getAstrometrySolution()
        matches = res.matches

        print 'Test1'

        logLevel = Log.DEBUG
        log = Log(Log.getDefaultLog(),
                  'meas.astrom',
                  logLevel)

        schema = matches[0].second.schema

        config = photocal.PhotoCalConfig()
        config.doWriteOutput = False    # schema is fixed because we already loaded the data
        task = photocal.PhotoCalTask(config=config, schema=schema)
        pCal = task.run(exposure=self.exposure, matches=matches)
        print "Ref flux fields list =", pCal.arrays.refFluxFieldList
        refFluxField = pCal.arrays.refFluxFieldList[0]

        # These are *all* the matches; we don't really expect to do that well.
        diff=[]
        for m in matches:
            refFlux = m[0].get(refFluxField) # reference catalog flux
            if refFlux <= 0:
                continue
            refMag = afwImage.abMagFromFlux(refFlux) # reference catalog mag
            instFlux = m[1].getPsfFlux()    #Instrumental Flux
            if instFlux <= 0:
                continue
            instMag = pCal.calib.getMagnitude(instFlux)     #Instrumental mag
            diff.append(instMag - refMag)
        diff = np.array(diff)

        self.assertTrue(len(diff) > 50)
        log.info('%i magnitude differences; mean difference %g; mean abs diff %g' %
                 (len(diff), np.mean(diff), np.mean(np.abs(diff))))
        self.assertLess(np.mean(diff), 0.6)

        # Differences of matched objects that were used in the fit.
        zp = pCal.calib.getMagnitude(1.)
        log.logdebug('zeropoint: %g' % zp)
        fitdiff = pCal.arrays.srcMag + zp - pCal.arrays.refMag
        log.logdebug('number of sources used in fit: %i' % len(fitdiff))
        log.logdebug('median diff: %g' % np.median(fitdiff))
        log.logdebug('mean diff: %g' % np.mean(fitdiff))
        log.logdebug('median abs(diff): %g' % np.median(np.abs(fitdiff)))
        log.logdebug('mean abs(diff): %g' % np.mean(np.abs(fitdiff)))

        # zeropoint: 31.3145
        # number of sources used in fit: 65
        # median diff: -0.009681
        # mean diff: 0.00331871
        # median abs(diff): 0.0368904
        # mean abs(diff): 0.0516589

        self.assertTrue(abs(zp - 31.3145) < 0.05)

        self.assertTrue(len(fitdiff) > 50)
        # These are kind of arbitrary
        self.assertTrue(abs(np.median(fitdiff)) < 0.02)
        self.assertTrue(abs(np.mean(fitdiff)) < 0.004)
        #
        self.assertTrue(np.median(np.abs(fitdiff)) < 0.04)
        self.assertTrue(np.mean(np.abs(fitdiff)) < 0.06)
            
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
