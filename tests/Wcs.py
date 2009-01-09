#!/usr/bin/env python
import os
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.afw.image as afwImage
import lsst.meas.astrom.net as net
import lsst.utils.tests as utilsTests

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
	"""Pass the positions of objects near the white dwarf GD66 and test that the correct
	position is returned

	"""

	gas = net.GlobalAstrometrySolution("backend.cfg")
	starlist = loadXYFromFile("gd66xy.txt")
	gas.setStarList(starlist)

	#To speed the test, tell the GAS what the size of the image is
	#The image is 1780 pixels on a side and covers half a square degree on the sky
	gas.setImageScaleArcsecPerPixel(.5*3600/1780.)

	flag = gas.blindSolve()

	if flag:
	    radec = gas.xy2RaDec(890, 890)

	    #RA
	    delta = math.fabs(radec[0]-80.159783)
	    if delta > 1e-6:
		return False

	    #Dec
	    delta = math.fabs(radec[1]-30.805249)
	    if delta > 1e-6:
		return False

	    return True
	else:
	    #If we didn't get a match, that's a failure
	    return False
	    
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
