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
import tempfile
import numpy as np
from math import ceil

import eups
import lsst.meas.astrom as measAstrom
import lsst.meas.astrom.net as net
import lsst.afw.detection as det
import lsst.afw.math as afwMath
import lsst.afw.image as afwImg
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.utils.tests as utilsTests
import lsst.pex.policy as pexPolicy
import lsst.meas.photocal as photocal
from astrometry.util.pyfits_utils import fits_table

from lsst.pex.exceptions import LsstCppException

import sourceSetIO as ssi


class TagAlongTest2(unittest.TestCase):

    def setUp(self):
        self.defaultPolicy = pexPolicy.Policy.createPolicy(pexPolicy.PolicyString(
        """#<?cfg paf policy?>     
        inputExposureKey: visitExposure
        inputSourceSetKey: sourceSet
        allowDistortion: true
        matchThreshold: 22
        blindSolve: false
        outputWcsKey: measuredWcs
        outputMatchListKey: matchList
        distanceForCatalogueMatchinArcsec: 1.0
        cleaningParameter: 3
        calculateSip: false
        numBrightStars: 75
        defaultFilterName: rmag
        sipOrder: 4
        """
        ))

        # Load fake image header from disk
        mypath = eups.productDir('meas_astrom')
        tests = os.path.join(mypath, 'tests')

        # Create fake blank image using a FITS image; this is only required
        # because measAstrom requires an Exposure.
        #self.tempimage = tempfile.NamedTemporaryFile(suffix='_img.fits')
        # delete=True (only available in python >= 2.6)

        (fd,self.tempimage) = tempfile.mkstemp(suffix='_img.fits')
        os.close(fd)
        f = open(self.tempimage, 'wb')

        hdr = open(os.path.join(tests, 'fpobjc-0745-3-0564-cut.hdr')).read()
        f.write(hdr)
        W,H = (2048,1489)
        bpp = 1
        #for i in range(H):
        #   self.tempimage.write(chr(0)*W)
        fitsbytes = 2880 * int(ceil(W*H*bpp / 2880.))
        f.write(chr(0)*fitsbytes)
        f.close()
        print 'Wrote fake image to temp file', self.tempimage
        base = self.tempimage.replace('_img.fits', '')
        #self.exposure = afwImg.ExposureF(base)

        self.exposure = afwImg.ExposureF(W, H)

        # Grab source list
        T = fits_table(os.path.join(tests, 'fpobjc-0745-3-0564-cut.fits'))
        sourceSet = det.SourceSet()
        for i,(x,y,f) in enumerate(zip(T.x, T.y, T.flux)):
            s = det.Source()
            s.setId(i)
            s.setFlagForDetection(8320)
            s.setRa(0.)
            s.setDec(0.)
            s.setXAstrom(x)
            s.setYAstrom(y)
            s.setPsfFlux(f)
            sourceSet.append(s)
        self.srcSet = sourceSet

        print 'Exposure image size: %i x %i' % (self.exposure.getWidth(), self.exposure.getHeight())
        print '%i sources' % len(sourceSet)

        # Set up local astrometry_net_data
        datapath = os.path.join(tests, 'astrometry_net_data', 'tagalong')
        # Work around lame scons bug (doesn't pass HOME)
        os.environ['HOME'] = 'iswheretheheartis'
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Couldn't set up local 'tagalong' version of astrometry_net_data (from path: %s): %s" % (datapath, reason))

    def tearDown(self):
        del self.defaultPolicy
        del self.exposure
        del self.srcSet

                        
    def test1(self):
		# this is how determineWcs.py initializes it...
        path = os.path.join(os.environ['ASTROMETRY_NET_DATA_DIR'], "metadata.paf")
        gas = measAstrom.net.GlobalAstrometrySolution(path)
        matchThreshold = self.defaultPolicy.get('matchThreshold')
        gas.setMatchThreshold(matchThreshold)

        matches, wcs = measAstrom.determineWcs(self.defaultPolicy, self.exposure,
                                               self.srcSet, solver=gas)

        print 'Grabbing index stars inside the solved field...'
        (xyz,radec,inds,tag) = gas.getIndexStarsInSolvedField(10.)
        print 'Found %i index stars in field' % len(xyz)
        print 'Found tag-along columns:', tag.keys()

        #(ra,dec,radius) = (-145, 53, 0.02)
        #print 'Searching RA,Dec %g,%g, radius %g deg' % (ra,dec,radius)
        #(xyz,radec,inds,tag) = gas.getIndexStars(ra,dec,radius)
        #print 'Found %i index stars' % len(xyz)
        #print 'Found tag-along columns:', tag.keys()
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TagAlongTest2)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
