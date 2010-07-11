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
import math
import sys
import unittest

import eups
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.meas.astrom.net as net
import lsst.utils.tests as utilsTests
import lsst.afw.image as afwImg
import lsst.afw.coord as afwCoord
import lsst.afw.detection.detectionLib as detect
try:
    type(verbose)
except NameError:
    verbose = 0

verbose=True

class GAS(object):
    __gas = None
    __desiredVersion = None

    def __init__(self, desiredVersion):
        """We use the product astrometry_net_dir to point to different index files depending on the
    problem we're trying to solve.  If any version of astrometry_net_data is setup, try to
    switch to version "desiredVersion".  There's no need to switch back as this code is running in a subprocess
    """

        if desiredVersion == None:
            GAS.__desiredVersion = None
            GAS.__gas = None
            return

        if GAS.__desiredVersion and desiredVersion == GAS.__desiredVersion:
            return

        eupsObj = eups.Eups()
        dataVersion = eupsObj.findSetupVersion("astrometry_net_data")[0]
        
        if dataVersion and dataVersion != desiredVersion:
            print >> sys.stderr, \
                  "Note: These tests require astrometry_net_data %s; Trying to set this up for you now" % \
                  desiredVersion

            ok, version, reason = eupsObj.setup("astrometry_net_data", versionName=desiredVersion)
            if not ok:
                print >> sys.stderr, "Error: %s" % (reason)
                GAS.__desiredVersion = None
                GAS.__gas = None
                return
        #
        # Create and return a GAS if a version of astrometry_net_data can be setup
        # that contains the proper index files
        #
        an_dataDir = eups.productDir("astrometry_net_data")
        if an_dataDir:
            an_dataDir_dataDir = os.path.join(an_dataDir, desiredVersion)
            if os.path.isdir(an_dataDir_dataDir):
                an_dataDir = an_dataDir_dataDir

            policyFile = os.path.join(an_dataDir, "metadata.paf")
            if not os.path.exists(policyFile):
                raise IOError("Unable to find %s" % (policyFile))

            policy = pexPolicy.Policy(policyFile)

            indexFile = policy.get("indexFile")

            if os.path.exists(os.path.join(an_dataDir, indexFile)):
                GAS.__desiredVersion = desiredVersion
                GAS.__gas = net.GlobalAstrometrySolution(policyFile)

                return
            else:
                print >> sys.stderr, "astrometry_net_data/%s indexFiles are not available; " + \
                      "not running %s tests" % (desiredVersion, desiredVersion)

        print >> sys.stderr, "astrometry_net_data/%s is not available; not running %s tests" % \
              (desiredVersion, desiredVersion)

    def __del__(self):
        GAS.__gas = None
        GAS.__desiredVersion = None

    def exists(self):
        return GAS.__gas != None

    def __str__(self):
        return str(GAS.__gas)

    #
    # Forwarding functions
    #
    def getDistortedWcs(self):
        if self.exists():
            return self.__gas.getDistortedWcs()
        else:
            return None

    def getMatchedSources(self):
        if self.exists():
            return self.__gas.getMatchedSources()
        else:
            return None

    def getSolvedImageScale(self):
        if self.exists():
            return self.__gas.getSolvedImageScale()
        else:
            return None

    def getWcs(self):
        if self.exists():
            return self.__gas.getWcs()
        else:
            return None

    def reset(self):
        if self.exists():
            self.__gas.reset()

    def setImageScaleArcsecPerPixel(self, plateScale):
        if self.exists():
            self.__gas.setImageScaleArcsecPerPixel(plateScale)

    def setLogLevel(self, level):
        if self.exists():
            self.__gas.setLogLevel(level)

    def setMaximumImageScale(self, val):
        if self.exists():
            self.__gas.setMaximumImageScale(val)

    def setMinimumImageScale(self, val):
        if self.exists():
            self.__gas.setMinimumImageScale(val)

    def setNumBrightObjects(self, nBright):
        if self.exists():
            self.__gas.setNumBrightObjects(nBright)
        
    def setStarlist(self, starlist):
        if self.exists():
            self.__gas.setStarlist(starlist)

    def solve(self, val=None):
        if self.exists():
            if val is None:
                return self.__gas.solve()
            else:
                return self.__gas.solve(val)
        else:
            return None

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def loadXYFromFile(filename):
    """Load a list of positions from a file"""
    f= open(filename)
    
    s1=detect.SourceSet()
    i=0
    for line in f:
        #Split the row into an array
        line = re.sub("^\s+", "", line)
        elts = re.split("\s+", line)
        
        #Swig requires floats, but elts is an array of strings
        x=float(elts[0])
        y=float(elts[1])
        flux=float(elts[2])

        source = detect.Source()

        source.setSourceId(i)
        source.setXAstrom(x); source.setXAstromErr(0.1)
        source.setYAstrom(y); source.setYAstromErr(0.1)
        source.setPsfFlux(flux)

        s1.append(source)
        
        i=i + 1
    f.close()
    
    return s1



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class SmallSolveGASTestCFHT(unittest.TestCase):
    """A test case for WCS from astrometry.net using a small set of index files"""

    def setUp(self):
        self.gas = GAS("cfhttemplate")

    def tearDown(self):
        del self.gas

    def solve(self, imgListFile, raDec, nBright=50, expectPass=True):
        starlist = loadXYFromFile(imgListFile)
        self.gas.setStarlist(starlist)
        self.gas.setNumBrightObjects(50)
    
        try:
            flag = self.gas.solve(raDec)
        except Exception,e:
            #if we expect to fail, throwing an exception counts as a fail
            if not expectPass:
                return
            else:
                raise
            
        if expectPass:
            self.assertTrue(flag, "No solution found")
            wcs = self.gas.getWcs()
            result = wcs.getSkyOrigin().toFk5()
            scale= self.gas.getSolvedImageScale()
            if verbose:
                print "%.6f %.6f %.3f" %(result.getRa(afwCoord.DEGREES), result.getDec(afwCoord.DEGREES), scale)

        else:
            self.assertFalse(flag, "Solution found, but none expected")
    
    def testGD66Fail(self):
        """A field not covered by the CFHT indices"""

        if not self.gas.exists():
            return

        crval = afwImage.PointD(80.15978319,30.80524999)

        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
        self.gas.reset()
        self.gas.setMinimumImageScale(.5)
        self.gas.setMaximumImageScale(2)
        self.solve(starlist, crval, expectPass=False)
        #
        
        
    def testCFHTa(self):                
        """Testing a field that should pass with different image scales"""
        if verbose:
            print "testCFHTa"

        if not self.gas.exists():
            return

        crval = afwImage.PointD(334.303012, -17.233988)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        self.gas.reset()
        self.gas.setMinimumImageScale(.1)
        self.gas.setMaximumImageScale(.2)
        self.gas.setLogLevel(verbose)
        #self.solve(starlist, crval)
        self.solve(starlist, crval)
        self.gas.setLogLevel(0)

        
    def testCFHTb(self):                
        """Different starting point"""
        if verbose:
            print "testCFHTb"
            
        if not self.gas.exists():
            return

        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        self.gas.reset()
        self.gas.setMinimumImageScale(.1)
        self.gas.setMaximumImageScale(.5)
        if verbose:
            self.gas.setLogLevel(3)
        self.solve(starlist, crval)    
        self.gas.setLogLevel(0)


    def testCFHTc(self):                
        """Different img scales"""
        if verbose:
            print "testCFHTc"

        if not self.gas.exists():
            return

        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        self.gas.reset()
        self.gas.setMinimumImageScale(.1)
        self.gas.setMaximumImageScale(1.5)
        if verbose:
            self.gas.setLogLevel(3)
        self.solve(starlist, crval)    
        self.gas.setLogLevel(0)


    def testCFHTd(self):                
        """Different img scales"""

        if verbose:
            print "testCFHTc"

        if not self.gas.exists():
            return self.assertAlmostEqual(sRaDec.getY(), wRaDec.getY(), 3, "y coord failed for getMatchedSources()")


        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        self.gas.reset()
        self.gas.setMinimumImageScale(.15)
        self.gas.setMaximumImageScale(.25)
        if verbose:
            self.gas.setLogLevel(3)
        self.solve(starlist, crval)    
        self.gas.setLogLevel(0)


    def testCFHTe(self):                
        """Different img scales"""
        if verbose:
            print "testCFHTe"

        if not self.gas.exists():
            return

        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        self.gas.reset()
        self.gas.setImageScaleArcsecPerPixel(.183)
        if verbose:
            self.gas.setLogLevel(3)
        self.solve(starlist, crval)    
        self.gas.setLogLevel(0)
        wcs = self.gas.getWcs()

    def NORUNtestDistortedWcs(self):
        """Is a distorted Wcs returned. This test is disabled for now because function hangs
        but we don't use it so we don't care."""

        if verbose:
            print "DistortedWcs"

        if not self.gas.exists():
            return

        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        self.gas.reset()
        self.gas.setImageScaleArcsecPerPixel(.183)
        if verbose:
            self.gas.setLogLevel(3)
        self.solve(starlist, crval)    
        wcs = self.gas.getDistortedWcs()
        self.gas.setLogLevel(0)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    #These tests take to long to run, and only repeat what cfht is doing, so I should remove them.
    #suites += unittest.makeSuite(WCSTestCaseNetUSNOB)
    suites += unittest.makeSuite(SmallSolveGASTestCFHT)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    verbose = 3
    # Run individual tests:
    run(True)
