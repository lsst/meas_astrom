#!/usr/bin/env python
import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.afw.image as afwImage
import lsst.utils.tests as utilsTests
import lsst.afw.image.imageLib as img
import lsst.afw.detection.detectionLib as detect
import lsst.pex.exceptions.exceptionsLib as exception

import lsst.meas.astrom.sip as sip


class Lsf2dTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testBadArgs1(self):
        """Check that code fails when order is too low"""
        
        x = [630, 1000, 1205, 1427]
        y = [535, 929, 623, 602]
        z = [.622, 1.109, 1.198, 1.429]
        s = [1,1,1,1]
        
        order=0

        self.assertRaises(exception.LsstCppException,
            sip.LeastSqFitter1dPoly, x, y,s,order )
        

    def testBadArgs2(self):
        """Check that code fails when not len(x) != len(y)"""
        
        x = [630, 1000, 1205, 1427]
        y = [535, 929, 623]
        z = [.622, 1.109, 1.198, 1.429]
        s = [1,1,1,1]
        
        order=2

        self.assertRaises(exception.LsstCppException,
            sip.LeastSqFitter1dPoly, x, y,s,order )

    def testBadArgs3(self):
        """Check that code fails when not len(x) != len(s)"""
        
        x = [630, 1000, 1205, 1427]
        y = [535, 929, 623, 602]
        z = [.622, 1.109, 1.198, 1.429]
        s = [1,1,1,1,1]
        
        order=0

        self.assertRaises(exception.LsstCppException,
            sip.LeastSqFitter1dPoly, x, y,s,order )


    def testBadArgs4(self):
        """Check that code fails when not order > number data points"""
        
        x = [630, 1000, 1205, 1427]
        y = [535, 929, 623, 602]
        z = [.622, 1.109, 1.198, 1.429]
        s = [1,1,1,1]
        
        order=5

        self.assertRaises(exception.LsstCppException,
            sip.LeastSqFitter1dPoly, x, y,s,order )

    def testFitLinearXSurface2(self):
        """Python equivalent of C++ test case"""
        
        x=[599.59899999999993, 1172.7726619709097, 512.51199999999994, 1083.6078436082901]
        y=[512.0, 539.77214401699996, 541.0, 562.09371856300004]
        z=[0.5989999999999327, 1.1716010609097793, 0.51199999999994361, 1.0825253182899814]
        s=[.1,.1,.1,.1]
        
        print " "
        lsf = sip.LeastSqFitter2dPoly(x, y, z, s, 2)
        lsf.printParams()
        
        print "Output:"
        for i in range(len(x)):
            print x[i], y[i], z[i], lsf.valueAt(x[i], y[i])
            self.assertAlmostEqual(z[i], lsf.valueAt(x[i], y[i]))
                

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(Lsf2dTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
