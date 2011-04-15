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
from astrometry.util.pyfits_utils import fits_table
from lsst.pex.logging import Log

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

        mypath = eups.productDir('meas_astrom')
        tests = os.path.join(mypath, 'tests')
        self.testdir = tests

        # Create a fake blank image of the appropriate size
        (W,H) = (2048,1489)
        self.exposure = afwImg.ExposureF(afwGeom.Extent2I(W, H))

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
        log = Log.getDefaultLog()
        gas = measAstrom.createSolver(self.defaultPolicy, log)

        astrom = measAstrom.determineWcs(self.defaultPolicy, self.exposure,
                                         self.srcSet, solver=gas)

        print 'Grabbing index stars inside the solved field...'
        (xyz,radec,inds,tag) = gas.getIndexStarsInSolvedField(10.)
        print 'Found %i index stars in field' % len(xyz)
        print 'Found tag-along columns:', tag.keys()

        self.assertEqual(123, len(xyz))
        self.assertEqual(145, len(tag.keys()))
        for k,v in tag.items():
            self.assertEqual(len(v), len(xyz))

        # This is from:
        # liststruc tests/astrometry_net_data/tagalong/tsobj-0745-3-40-0564.index+13 | awk '{printf("\"%s\", ", $2);}' > cols
        keys = ["run", "camCol", "rerun", "field", "parent", "id", "nchild", "objc_type", "objc_prob_psf", "catID", "objc_flags", "objc_flags2", "objc_rowc", "objc_rowcErr", "objc_colc", "objc_colcErr", "rowv", "rowvErr", "colv", "colvErr", "rowc", "rowcErr", "colc", "colcErr", "sky", "skyErr", "psfCounts", "psfCountsErr", "fiberCounts", "fiberCountsErr", "petroCounts", "petroCountsErr", "petroRad", "petroRadErr", "petroR50", "petroR50Err", "petroR90", "petroR90Err", "Q", "QErr", "U", "UErr", "M_e1", "M_e2", "M_e1e1Err", "M_e1e2Err", "M_e2e2Err", "M_rr_cc", "M_rr_ccErr", "M_cr4", "M_e1_psf", "M_e2_psf", "M_rr_cc_psf", "M_cr4_psf", "iso_rowc", "iso_rowcErr", "iso_rowcGrad", "iso_colc", "iso_colcErr", "iso_colcGrad", "iso_a", "iso_aErr", "iso_aGrad", "iso_b", "iso_bErr", "iso_bGrad", "iso_phi", "iso_phiErr", "iso_phiGrad", "r_deV", "r_deVErr", "ab_deV", "ab_deVErr", "phi_deV", "phi_deVErr", "counts_deV", "counts_deVErr", "r_exp", "r_expErr", "ab_exp", "ab_expErr", "phi_exp", "phi_expErr", "counts_exp", "counts_expErr", "counts_model", "counts_modelErr", "texture", "star_L", "star_lnL", "exp_L", "exp_lnL", "deV_L", "deV_lnL", "fracPSF", "flags", "flags2", "type", "prob_psf", "nprof", "profMean", "profErr", "status", "lambda", "eta", "l", "b", "offsetRa", "offsetDec", "primTarget", "secTarget", "reddening", "propermotionmatch", "propermotiondelta", "propermotion", "propermotionangle", "usnoBlue", "usnoRed", "firstMatch", "firstId", "firstLambda", "firstEta", "firstDelta", "firstPeak", "firstInt", "firstRms", "firstMajor", "firstMinor", "firstPa", "rosatMatch", "rosatDelta", "rosatPosErr", "rosatCps", "rosatCpsErr", "rosatHr1", "rosatHr1Err", "rosatHr2", "rosatHr2Err", "rosatExt", "rosatExtLike", "rosatDetectLike", "rosatExposure", "priority", "matchid", "rmag"]
        for k in keys:
            self.assertTrue(k in tag)

        # Check a single row... id=3
        #X = (tag['id'])[:]
        #X.sort()
        #print X
        I = tag['id'].index(3)
        print I
        # I got these values from:
        #    sdss_das.py -r 745 -f 564 -c 3 -b r tsObj
        #    python tests/testTagAlong2.py truetags > tests/truetags.py
        # And then I checked them against:
        #    tablist tsObj-000745-3-40-0564.fit"[id==3]"

        # By the way, the row id==3 is row 63 in the index's tag-along table:
        #    tablist tests/astrometry_net_data/tagalong/tsobj-0745-3-40-0564.index+13"[#row==63]"

        # Column "priority":
        # the tsobj file has the value 0x80000000 and type J (signed int);
        # it seems that both pyfits and tablist handle this incorrectly and produce 0.
        # This code, on the other hand, correctly produces -2147483648

        # To see this: the tsobj file has 2732 bytes/row, and the data chunk starts at 40320.
        #   "priority" is the second-last column, format 1J, followed by format 50J.
        #  id==3 is in row 3, thus is at byte offset  '%x'%(40320 + (3*2732) - (51*4)) --> 0xbcb8
        # hexedit tsObj-000745-3-40-0564.fit
        # od -N 4 -j 48312 -t xC tsObj-000745-3-40-0564.fit
        # --> 0136270    80  00  00  00

        # In the index file, it starts at 123840, is followed by 50J + 1D,
        # and rows are 2724 bytes, and it's in row 63:
        # 123840 + (63*2724) - (51*4 + 1*8)   = 295240
        # fitsgetext -i tests/astrometry_net_data/tagalong/tsobj-0745-3-40-0564.index
        # --> Extension 13 : header start 86400 , length 37440 ; data start 123840 , length 336960 .
        # modhead tests/astrometry_net_data/tagalong/tsobj-0745-3-40-0564.index+13 NAXIS1
        # --> NAXIS1  =                 2724 / Bytes in row
        # od -N 4 -j 295240 -t xC tests/astrometry_net_data/tagalong/tsobj-0745-3-40-0564.index
        # --> 1100510 80 00 00 00

        sys.path.append(self.testdir)
        from truetags import truetags

        for k,truev in truetags.items():
            self.assertTrue(k in tag)
            v = (tag[k])[I]
            #print k, type(v), type(truev)
            if k == 'priority':
                print 'Skipping bad column "priority"'
                continue
            if type(v) is float:
                if truev == 0.0:
                    self.assertEqual(v, truev)
                else:
                    #print k, v, truev, abs(v - truev) / abs(truev)
                    self.assertTrue(abs(v - truev) / abs(truev) < 1e-15)
            else:
                self.assertEqual(v, truev)
        
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
    if 'truetags' in sys.argv:
        import numpy
        T = fits_table('tsObj-000745-3-40-0564.fit')
        T = T[T.id == 3]
        T = T[0]
        print 'truetags={'
        # this is same as above, but minus "rmag", since we hacked that in.
        keys = ["run", "camCol", "rerun", "field", "parent", "id", "nchild", "objc_type", "objc_prob_psf", "catID", "objc_flags", "objc_flags2", "objc_rowc", "objc_rowcErr", "objc_colc", "objc_colcErr", "rowv", "rowvErr", "colv", "colvErr", "rowc", "rowcErr", "colc", "colcErr", "sky", "skyErr", "psfCounts", "psfCountsErr", "fiberCounts", "fiberCountsErr", "petroCounts", "petroCountsErr", "petroRad", "petroRadErr", "petroR50", "petroR50Err", "petroR90", "petroR90Err", "Q", "QErr", "U", "UErr", "M_e1", "M_e2", "M_e1e1Err", "M_e1e2Err", "M_e2e2Err", "M_rr_cc", "M_rr_ccErr", "M_cr4", "M_e1_psf", "M_e2_psf", "M_rr_cc_psf", "M_cr4_psf", "iso_rowc", "iso_rowcErr", "iso_rowcGrad", "iso_colc", "iso_colcErr", "iso_colcGrad", "iso_a", "iso_aErr", "iso_aGrad", "iso_b", "iso_bErr", "iso_bGrad", "iso_phi", "iso_phiErr", "iso_phiGrad", "r_deV", "r_deVErr", "ab_deV", "ab_deVErr", "phi_deV", "phi_deVErr", "counts_deV", "counts_deVErr", "r_exp", "r_expErr", "ab_exp", "ab_expErr", "phi_exp", "phi_expErr", "counts_exp", "counts_expErr", "counts_model", "counts_modelErr", "texture", "star_L", "star_lnL", "exp_L", "exp_lnL", "deV_L", "deV_lnL", "fracPSF", "flags", "flags2", "type", "prob_psf", "nprof", "profMean", "profErr", "status", "lambda", "eta", "l", "b", "offsetRa", "offsetDec", "primTarget", "secTarget", "reddening", "propermotionmatch", "propermotiondelta", "propermotion", "propermotionangle", "usnoBlue", "usnoRed", "firstMatch", "firstId", "firstLambda", "firstEta", "firstDelta", "firstPeak", "firstInt", "firstRms", "firstMajor", "firstMinor", "firstPa", "rosatMatch", "rosatDelta", "rosatPosErr", "rosatCps", "rosatCpsErr", "rosatHr1", "rosatHr1Err", "rosatHr2", "rosatHr2Err", "rosatExt", "rosatExtLike", "rosatDetectLike", "rosatExposure", "priority", "matchid"]
        #for k in T.columns():
        for k in keys:
            v = T.getcolumn(k)
            print "'%s':" % k,
            #print type(v)
            if isinstance(v, numpy.ndarray):
                print v.tolist(),
            elif type(v) in [numpy.float32, numpy.float64]:
                print '%.16g' % v,
            else:
                print v,
            print ','
            #print ' # ', type(v)
        print '}'
        sys.exit(0)

    run(True)
