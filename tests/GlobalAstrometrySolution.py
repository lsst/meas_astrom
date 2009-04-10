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

dataDir = eups.productDir("astrometry_net_data")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run these tests")

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


class WCSTestCaseNet(unittest.TestCase):
    """A test case for WCS from astrometry.net"""


    def setUp(self):
        pass

    def tearDown(self):
        gas.reset()


    def solveOrVerify(self, starlist, crval, crpix,  plateScale=0, verify=False):
        """Test the solve() function
        
            Input: 
            starlist    List of objects as returned by loadXYFromFile
            crval       lsst.afw.image.PointD ra/dec of a known position
                        on the image
            crpix       lsst.afw.image.PointD xy of a known position
                        on the image
            plateScale  Size of image in arcsec/pixel. Specifing this 
                        dramatically improves search time
            verify      If True, crval is passed to to solve() to speed up the match,
                        otherwise a solution is found with no inital guess at the position
        """
        
        #Set plate scale
        if plateScale > 0:
            gas.setImageScaleArcsecPerPixel(plateScale)
        
        #Set starlist    
        starlist = loadXYFromFile(starlist)
        gas.setStarlist(starlist)
        
        #Run solver
        if verify:
            flag = gas.solve(crval)
        else:
            flag = gas.solve()

        if flag:
            #Test xy->radec
            radec = gas.xyToRaDec(crpix.getX(), crpix.getY())
            self.assertAlmostEqual(radec.getX(), crval.getX(), 6, "Ra doesn't match")
            self.assertAlmostEqual(radec.getY(), crval.getY(), 6, "Dec doesn't match")

            #Test the reverse operation
            xy = gas.raDecToXY(crval.getX(), crval.getY())
            self.assertAlmostEqual(xy.getX(), crpix.getX(), 2, "X pos doesn't match")
            self.assertAlmostEqual(xy.getY(), crpix.getY(), 2, "Y pos doesn't match")

        else:
            #If we didn't get a match, that's a failure
            self.assertEqual(flag, 1, "Failed to find a match")
        

    def solveWcs(self, wcsPtr, starlist):
        """Test that the solve(wcs) function works correctly.
        """

        gas.setStarlist(starlist)
        gas.setLogLevel(2)
        return gas.solve(wcsPtr)


    def testSolveGD66Wcs(self):
        """Run solveWcs on GD66"""

        #This Wcs solution is taken from the header of the image.
        crval = afwImage.PointD(80.1601717554, 30.8060971075)
        crpix = afwImage.PointD(890,890)
        wcsPtr = afwImage.createWcs(crval, crpix, -0.0002802350, -0.0000021800, -0.0000022507, 0.0002796878)

        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
        starlist = loadXYFromFile(starlist)

        flag = self.solveWcs(wcsPtr, starlist)

        if flag:
            radec = gas.xyToRaDec(890,890)
            self.assertAlmostEqual(crval.getX(), radec.getX(), 6, "Ra doesn't match")
            self.assertAlmostEqual(crval.getY(), radec.getY(), 6, "Dec doesn't match")
        else:
            self.assertEqual(flag, 1, "Failed to find a match")
        
        gas.reset()
        
    def testSolveGD66(self):
        """Pass the positions of objects near the white dwarf GD66 and test that the correct position is returned
    """
        crval = afwImage.PointD(80.15978319, 30.80524999)
        crpix = afwImage.PointD(890,890)
        listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
        #To speed the test, tell the GAS what the size of the image is
        #The image is 1780 pixels on a side and covers half a square degree 
        #on the sky
        plateScale = .5*3600/1780.
        self.solveOrVerify(listFile, crval, crpix, plateScale)
        gas.reset()

    def testSolveG117(self):
        crval = afwImage.PointD(141.063590, +35.280919)
        crpix = afwImage.PointD(446, 447)
        listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "g117.xy.txt")
        #To speed the test, tell the GAS what the size of the image is
        #The image is 1780 pixels on a side and covers half a square degree 
        #on the sky
        plateScale = .5*3600/1780.
        self.solveOrVerify(listFile, crval, crpix, plateScale)
    #
    def testVerifyCFHTField(self):
        crval = afwImage.PointD(334.303012, -17.233988)
        crpix = afwImage.PointD(512,512)
        listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")
        #To speed the test, tell the GAS what the size of the image is
        #The image is 1780 pixels on a side and covers half a square degree 
        #on the sky
        plateScale = .185
        self.solveOrVerify(listFile, crval, crpix, plateScale, verify=True)

#
    def testMultiple(self):
        """Test that solver can handle doing two solves in a row"""
        
        
        #GD66
        crval = afwImage.PointD(80.15978319,30.80524999)
        crpix = afwImage.PointD(890,890)
        listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
        #To speed the test, tell the GAS what the size of the image is
        #The image is 1780 pixels on a side and covers half a square degree 
        #on the sky
        plateScale = .5*3600/1780.
        self.solveOrVerify(listFile, crval, crpix, plateScale)

        gas.reset()
        
        #G117
        crval = afwImage.PointD(141.063590, +35.280919)
        crpix = afwImage.PointD(446, 447)
        listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "g117.xy.txt")
        #To speed the test, tell the GAS what the size of the image is
        #The image is 1780 pixels on a side and covers half a square degree 
        #on the sky
        plateScale = .5*3600/1780.
        self.solveOrVerify(listFile, crval, crpix, plateScale)

    
        def testGetWcs(self):
            """Test the functions that return wcs structures for a field"""
     
            crval = afwImage.PointD(80.15978319,30.80524999)
            crpix = afwImage.PointD(890,890)
            listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
            plateScale = .5*3600/1780.
            
            flag = self.solveOrVerify(listFile, crval, crpix, plateScale)
            
            if flag:
                wcs1 = gas.getWcs();
                radec = wcs1.xyToRaDec(crpix.getX(), crpix.getY())
                
                #Test xy->radec
                radec = gas.xyToRaDec(crpix.getX(), crpix.getY())
                self.assertAlmostEqual(radec.getX(), crval.getX(), 6, "Ra doesn't match")
                self.assertAlmostEqual(radec.getY(), crval.getY(), 6, "Dec doesn't match")
    
                #Test the reverse operation
                xy = gas.raDecToXY(crval.getX(), crval.getY())
                self.assertAlmostEqual(xy.getX(), crpix.getX(), 2, "X pos doesn't match")
                self.assertAlmostEqual(xy.getY(), crpix.getY(), 2, "Y pos doesn't match")
    
            else:
                #If we didn't get a match, that's a failure
                self.assertEqual(flag, 1, "Failed to find a match")
            

        def testGetDistortedWcs(self):
            """Test the functions that return wcs structures for a field"""
     
            crval = afwImage.PointD(80.15978319,30.80524999)
            crpix = afwImage.PointD(890,890)
            listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
            plateScale = .5*3600/1780.
            
            flag = self.solveOrVerify(listFile, crval, crpix, plateScale)
            
            if flag:
                wcs1 = gas.getDistortedWcs();
                radec = wcs1.xyToRaDec(crpix.getX(), crpix.getY())
                print radec
                
                #Test xy->radec
                radec = gas.xyToRaDec(crpix.getX(), crpix.getY())
                self.assertAlmostEqual(radec.getX(), crval.getX(), 6, "Ra doesn't match")
                self.assertAlmostEqual(radec.getY(), crval.getY(), 6, "Dec doesn't match")
    
                #Test the reverse operation
                xy = gas.raDecToXY(crval.getX(), crval.getY())
                self.assertAlmostEqual(xy.getX(), crpix.getX(), 2, "X pos doesn't match")
                self.assertAlmostEqual(xy.getY(), crpix.getY(), 2, "Y pos doesn't match")
    
            else:
                #If we didn't get a match, that's a failure
                self.assertEqual(flag, 1, "Failed to find a match")


        def testWcsSinglePixelOffset(self):

            mastrom = eups.productDir("meas_astrom")
            imageFilename = "gd66.fits"
            starlist = os.path.join(mastrom, "tests", "gd66.xy.txt")


            #filename = eups.productDir("afwdata")
            #filename = os.path.join(filename, "CFHT", "D4", imageFilename)
            filename = os.path.join(mastrom, "tests", "gd66.fits")

            starlist = loadXYFromFile(starlist)

            #Get Wcs from image header
            exposure = readExposure(filename)
            origWcs = exposure.getWcs()

            #Get Wcs from astrometry.net
            gasWcs = gas.solve(starlist, origWcs)

            #Pick an radec. The xy values corresponding to this radec should
            #differ by sqrt(2) between the two wcs'. Also, the values for
            #gasWcs should be larger in both axes
            radec = afwImage.PointD(80.139800, +30.7864306)
            origPix = origWcs.raDecToXY(radec)
            gasPix = gasWcs.raDecToXY(radec)

            
            self.assertTrue(origPix.getX() <= gasPix.getX(), "GAS Wcs moved in wrong direction in X")
            self.assertTrue(origPix.getY() <= gasPix.getY(), "GAS Wcs moved in wrong direction in Y")
            
            ds = origPix-gasPix
            ds = math.hypot(ds.getX(), ds.getY() )
            self.assertAlmostEqual(ds, Math.sqrt(2), 1, "Distance moved not 1 pixel")    
    
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


#Create a globally accessible instance of a GAS
policyFile=eups.productDir("astrometry_net_data")
policyFile=os.path.join(policyFile, "metadata.paf")
gas = net.GlobalAstrometrySolution(policyFile)
 
if __name__ == "__main__":
    #print "Warning, tests turned off"
    #return 0
    run(True)
