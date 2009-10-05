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


class WCSTestCaseNet(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testBadArgs1(self):
        """Check that code fails when order is too low"""
        
        x = [630, 1000, 1205, 1427]
        y = [.622, 1.109, 1.198, 1.429]
        s = [1,1,1,1]
        
        order=0

        self.assertRaises(exception.LsstCppException,
            sip.LeastSqFitter1dPoly, x, y,s,order )
        

    def testBadArgs2(self):
        """Check that code fails when not len(x) != len(y)"""
        
        x = [630, 1000, 1205, 1427]
        y = [.622, 1.109, 1.198]
        s = [1,1,1,1]
        
        order=0

        self.assertRaises(exception.LsstCppException,
            sip.LeastSqFitter1dPoly, x, y,s,order )

    def testBadArgs3(self):
        """Check that code fails when not len(x) != len(s)"""
        
        x = [630, 1000, 1205, 1427]
        y = [.622, 1.109, 1.198, 1.429]
        s = [1,1,1,1,1]
        
        order=0

        self.assertRaises(exception.LsstCppException,
            sip.LeastSqFitter1dPoly, x, y,s,order )


    def testBadArgs4(self):
        """Check that code fails when not order > number data points"""
        
        x = [630, 1000, 1205, 1427]
        y = [.622, 1.109, 1.198, 1.429]
        s = [1,1,1,1]
        
        order=5

        self.assertRaises(exception.LsstCppException,
            sip.LeastSqFitter1dPoly, x, y,s,order )

    #def testConst1(self):
        #"""Check that code fits a dc offset correctly (1)"""
        #
        #x = [630, 1000, 1205, 1427]
        #y = [.622, 1.109, 1.198, 1.429]
        #s = [1,1,1,1]
        #
        #order=1
#
        #lsf = sip.LeastSqFitter1dPoly(x, y,s,order)
        #
        #for i in range(1, len(x)):
            #self.assertAlmostEqual(lsf.valueAt(x[0]), lsf.valueAt(x[i]), 4)
            ##print x[i], y[i], lsf.valueAt(x[i])
 #       

    #def testCompareToCpp1(self):
        #"""Confirm that I get the same behaviour in Python as I 
        #do in the C++ test fitLinear3
        #"""
        #
        ##This test arises because I was having trouble using Eigen's
        ##svd inversion. fitLinear3 in C++ worked fine, but this test
        ##behaved differently, and I tracked the result down to
        ##the svd routines. 
        #
        #x = [689.301136505, 1112.8573687, 1386.67168477]
        #y = [0.66911456573, 1.1147439759, 1.39597284177]
        #s = [1,1,1]
        #
        #order=2
        #print
        #lsf = sip.LeastSqFitter1dPoly(x, y,s,order)
#
        #self.assertAlmostEqual(lsf.valueAt(x[0]), 0.670187891961, 4 )
        #self.assertAlmostEqual(lsf.valueAt(x[1]), 1.11201034929, 4 )
        #self.assertAlmostEqual(lsf.valueAt(x[2]), 1.39763314215, 4 )
        

    def testCompareToCpp2(self):
        """Confirm that I get the same behaviour in Python as I 
        do in the C++ test fitLinear3
        """
        
        #This test arises because I was having trouble using Eigen's
        #svd inversion. fitLinear3 in C++ worked fine, but this test
        #behaved differently, and I tracked the result down to
        #the svd routines. 
        
        x = [628.857680996, 995.008255088, 1203.39412154, 1425.1404727]
        y = [0.616728229875, 1.01887634344, 1.19679830465, 1.42873084062]
        s = [1,1,1, 1]
        
        order=3
        print
        lsf = sip.LeastSqFitter1dPoly(x, y,s,order)
        
        for i in range(len(x)):
            f = lsf.valueAt(x[i])
            print "%.1f %.3f %.3f" % (x[i], y[i], f)

        self.assertAlmostEqual(lsf.valueAt(x[0]), y[0], 1 )
        self.assertAlmostEqual(lsf.valueAt(x[1]), y[1], 1 )
        self.assertAlmostEqual(lsf.valueAt(x[2]), y[2], 1 )
        
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
