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
import lsst.meas.astrom.net as net
import lsst.afw.image.imageLib as img
import lsst.afw.detection.detectionLib as detect
try:
    type(verbose)
except NameError:
    verbose = 0

if False:
    dataDir = eups.productDir("afwdata")
    if not dataDir:
	raise RuntimeError("Must set up afwdata to run these tests")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class WCSTestCaseNet(unittest.TestCase):
    """A test case for WCS from astrometry.net"""

    def setUp(self):
	self.gas = net.GlobalAstrometrySolution()

    def tearDown(self):
        del self.gas

    def testCleanup(self):
        """Test that tearDown does"""
        pass

    def testFindGD66(self):
	"""Pass the positions of objects near the white dwarf GD66 and test that the correct position is returned
	"""
	gas = net.GlobalAstrometrySolution();

	#Read in the indices (i.e the files containing the positions of known asterisms
	#and add them to the gas object
        print "Loading indices..."
        indices=glob.glob( os.path.join(eups.productDir("astrometry_net_data"), "index-2*.fits") )
        gas.setLogLevel(2)
        for f in indices:
            gas.addIndexFile(f)
        print gas.getNumIndices()
	
	#Read in a list of object positions in an image
	starlist = loadXYFromFile(os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt"))
	gas.setStarlist(starlist)

	#To speed the test, tell the GAS what the size of the image is
	#The image is 1780 pixels on a side and covers half a square degree on the sky
	gas.setImageScaleArcsecPerPixel(.5*3600/1780.)

	flag = gas.blindSolve()

	if flag:
	    radec = gas.xy2RaDec(890, 890)

	    self.assertAlmostEqual(radec[0], 80.15978319, 7, "Ra doesn't match")
            self.assertAlmostEqual(radec[1], 30.80524999, 7, "Dec doesn't match")
	else:
	    #If we didn't get a match, that's a failure
            self.assertEqual(flag, 1, "Failed to find a match")
	    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def loadXYFromFile(filename):
    """Load a list of positions from a file"""
    f= open(filename)
    
    s1=detect.SourceContainer()
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
        source.setPsfMag(flux)

        s1.append(source)
        
        i=i + 1
    f.close()
    
    return s1

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(WCSTestCaseNet)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
