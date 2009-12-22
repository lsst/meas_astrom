
#!/usr/bin/env python
import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.afw.image as afwImage
import lsst.meas.astrom.net as net
import lsst.utils.tests as utilsTests
import lsst.afw.image.imageLib as img
import lsst.afw.detection.detectionLib as detect
try:
    type(verbose)
except NameError:
    verbose = 0


eupsObj = eups.Eups()
dataVersion=eupsObj.findSetupVersion("astrometry_net_data")[0]
if dataVersion != "cfhttemplate":
    print "Warning: These tests require astrometry_net_data cfhttemplate"
    print "Setting this up for you now"
    try:
        eups.setup(eupsObj, "astrometry_net_data", version="cfhttemplate")
    except RuntimeError, e:
        print e
        raise RuntimeError("Failed to set up astrometry_net_data cfhttemplate")
    
dataDir = eups.productDir("astrometry_net_data")


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

class SmallSolveGASTest(unittest.TestCase):
    """A test case for WCS from astrometry.net"""


    def setUp(self):
        pass

    def tearDown(self):
        gas.reset()

    def solve(self, imgListFile, raDec, nBright=50, expectPass=True):
        starlist = loadXYFromFile(imgListFile)
        gas.setStarlist(starlist)
        gas.setNumBrightObjects(50)
    
        flag = gas.solve(raDec)
        if expectPass:
            self.assertTrue(flag, "No solution found")
            wcs = gas.getWcs()
            result = wcs.getOriginRaDec()
            scale= gas.getSolvedImageScale()
            if verbose:
                print "%.6f %.6f %.3f" %(result.getX(), result.getY(), scale)

        else:
            self.assertFalse(flag, "Solution found, but none expected")
    
    def testGD66Fail(self):
        """A field not covered by the CFHT indices"""
        crval = afwImage.PointD(80.15978319,30.80524999)

        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
        gas.reset()
        gas.setMinimumImageScale(.5)
        gas.setMaximumImageScale(2)
        self.solve(starlist, crval, expectPass=False)
        #
        
        
    def testCFHTa(self):                
        """Testing a field that should pass with different image scales"""
        if verbose:
            print "testCFHTa"
        crval = afwImage.PointD(334.303012, -17.233988)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        gas.reset()
        gas.setMinimumImageScale(.1)
        gas.setMaximumImageScale(.2)
        if verbose:
            gas.setLogLevel(3)
        self.solve(starlist, crval)    
        gas.setLogLevel(0)

        #
    def testCFHTb(self):                
        """Different starting point"""
        if verbose:
            print "testCFHTb"
            
        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        gas.reset()
        gas.setMinimumImageScale(.1)
        gas.setMaximumImageScale(.5)
        if verbose:
            gas.setLogLevel(3)
        self.solve(starlist, crval)    
        gas.setLogLevel(0)


    def testCFHTc(self):                
        """Different img scales"""
        if verbose:
            print "testCFHTc"

        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        gas.reset()
        gas.setMinimumImageScale(.1)
        gas.setMaximumImageScale(1.5)
        if verbose:
            gas.setLogLevel(3)
        self.solve(starlist, crval)    
        gas.setLogLevel(0)


    def testCFHTd(self):                
        if verbose:
            print "testCFHTc"
        """Different img scales"""
        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        gas.reset()
        gas.setMinimumImageScale(.15)
        gas.setMaximumImageScale(.25)
        if verbose:
            gas.setLogLevel(3)
        self.solve(starlist, crval)    
        gas.setLogLevel(0)


    def testCFHTe(self):                
        if verbose:
            print "testCFHTe"
        """Different img scales"""
        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        gas.reset()
        gas.setImageScaleArcsecPerPixel(.183)
        if verbose:
            gas.setLogLevel(3)
        self.solve(starlist, crval)    
        gas.setLogLevel(0)
        wcs = gas.getWcs()

    def testDistortedWcs(self):
        """Is a distorted Wcs returned"""
        if verbose:
            print "DistortedWcs"
        """Different img scales"""
        crval = afwImage.PointD(334.303215, -17.329315)
        #Set starlist    
        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")

        gas.reset()
        gas.setImageScaleArcsecPerPixel(.183)
        if verbose:
            gas.setLogLevel(3)
        self.solve(starlist, crval)    
        wcs = gas.getDistortedWcs()
        gas.setLogLevel(0)


    def testSolveWcs(self):
        """Can't be tested until I figure out the LSST way of interfacing Eigen Matrices"""
        pass
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(SmallSolveGASTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)


#Create a globally accessible instance of a GAS
policyFile=eups.productDir("astrometry_net_data")
policyFile=os.path.join(policyFile, "metadata.paf")
gas = net.GlobalAstrometrySolution(policyFile)
 
if __name__ == "__main__":
    run(True)
    if dataVersion != "cfhttemplate":
        eups.setup(eupsObj, "astrometry_net_data", dataVersion)  #Restore old a_n_data package
