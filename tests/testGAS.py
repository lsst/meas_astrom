#!/usr/bin/env python
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
import lsst.afw.detection.detectionLib as detect
try:
    type(verbose)
except NameError:
    verbose = 0

class GAS(object):
    __gas = None
    __desiredVersion = None

    def __init__(self, desiredVersion=None):
        """We use the product astrometry_net_dir to point to different index files depending on the
    problem we're trying to solve.  If any version of astrometry_net_data is setup, try to
    switch to version "desiredVersion".  There's no need to switch back as this code is running in a subprocess
    """

        if desiredVersion == None:
            desiredVersion = GAS.__desiredVersion

        if GAS.__desiredVersion and desiredVersion == GAS.__desiredVersion:
            return

        eupsObj = eups.Eups()
        dataVersion = eupsObj.findSetupVersion("astrometry_net_data")[0]
        
        if dataVersion and dataVersion != desiredVersion:
            print >> sys.stderr, \
                  "Note: These tests require astrometry_net_data %s; Trying to set this up for you now" % \
                  desiredVersion

            try:
                eupsObj.setup("astrometry_net_data", versionName=desiredVersion)
            except Exception, e:
                print >> sys.stderr, "Error setting up %s: %s" % (desiredVersion, e)
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
                raise RuntimeError("astrometry_net_data/%s indexFiles are not available; " + 
                                   "not running %s tests" % (desiredVersion, desiredVersion))

        raise RuntimeError, ("astrometry_net_data/%s is not available; not running %s tests" %
                             (desiredVersion, desiredVersion))

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


class WCSTestCaseNetUSNOB(unittest.TestCase):
    """A test case for WCS from astrometry.net"""

    def setUp(self):
        self.gas = GAS("usnob")

    def tearDown(self):
        self.gas.reset()

    def solveOrVerify(self, starlist, crval, crpix,  plateScale=0, nBright=50, verify=False):
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

        #Set starlist    
        starlist = loadXYFromFile(starlist)
        self.gas.setStarlist(starlist)
        
        #Set plate scale
        if plateScale > 0:
            self.gas.setImageScaleArcsecPerPixel(plateScale)
        
        self.gas.setNumBrightObjects(nBright)
        
        
        #Run solver
        if verify:
            flag = self.gas.solve(crval)
        else:
            flag = self.gas.solve()

        if flag:
            #Test xy->radec
            wcs = self.gas.getWcs()
            radec = wcs.xyToRaDec(crpix.getX(), crpix.getY())
            self.assertAlmostEqual(radec.getX(), crval.getX(), 6, "Ra doesn't match")
            self.assertAlmostEqual(radec.getY(), crval.getY(), 6, "Dec doesn't match")

            #Test the reverse operation
            xy = wcs.raDecToXY(crval.getX(), crval.getY())
            self.assertAlmostEqual(xy.getX(), crpix.getX(), 2, "X pos doesn't match")
            self.assertAlmostEqual(xy.getY(), crpix.getY(), 2, "Y pos doesn't match")

        else:
            #If we didn't get a match, that's a failure
            self.assertTrue(flag, "Failed to find a match")
        

    def solveWcs(self, wcsPtr, starlist):
        """Test that the solve(wcs) function works correctly.
        """
        self.gas.setStarlist(starlist)
        self.gas.setLogLevel(verbose)
        return self.gas.solve(wcsPtr)

    def testSolveGD66(self):
        """Pass the positions of objects near the white dwarf GD66 and test that the correct position is returned
    """
        if verbose:
            print "testSolveGD66Wcs..."

        if not self.gas.exists():
            return

        crval = afwImage.PointD(80.15978319, 30.80524999)
        crpix = afwImage.PointD(890,890)
        listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
        #To speed the test, tell the GAS what the size of the image is
        #The image is 1780 pixels on a side and covers half a square degree 
        #on the sky
        plateScale = .5*3600/1780.
        self.solveOrVerify(listFile, crval, crpix, plateScale)
        self.gas.reset()

#
    def testSolveGD66Wcs(self):
        """Run solveWcs on GD66. Also does a sanity check on the list
        returned by getMatchedSources()"""
        if verbose:
            print "testSolveGD66Wcs..."

        if not self.gas.exists():
            return

        crval = afwImage.PointD(80.15978319, 30.80524999)
        crpix = afwImage.PointD(890,890)
        wcsPtr = afwImage.createWcs(crval, crpix, -0.0002802350, -0.0000021800, -0.0000022507, 0.0002796878)

        starlist = os.path.join(eups.productDir("meas_astrom"), "tests", "gd66.xy.txt")
        starlist = loadXYFromFile(starlist)

        flag = self.solveWcs(wcsPtr, starlist)

        solvedWcs = self.gas.getWcs()
        if flag:
            radec = solvedWcs.xyToRaDec(890,890)
            strr="rd= (%.6f, %.6f): crval=(%.6f %.6f)" %(radec.getX(), radec.getY(), crval.getX(), crval.getY())
            self.assertAlmostEqual(crval.getX(), radec.getX(), 6, "Ra doesn't match: %s" %(strr))
            self.assertAlmostEqual(crval.getY(), radec.getY(), 6, "Dec doesn't match: %s" %(strr))
        else:
            self.assertEqual(flag, 1, "Failed to find a match")
            
        sourceSet = self.gas.getMatchedSources()
        for i in range(len(sourceSet)):
            x = sourceSet[i].getXAstrom()
            y = sourceSet[i].getYAstrom()
            sXY = afwImage.PointD(x,y)
            
            x = sourceSet[i].getRa()
            y = sourceSet[i].getDec()
            sRaDec  = afwImage.PointD(x,y)
            
            wRaDec = solvedWcs.xyToRaDec(sXY)
            
            self.assertAlmostEqual(sRaDec.getX(), wRaDec.getX(), 3, "x coord failed for getMatchedSources()")
            self.assertAlmostEqual(sRaDec.getY(), wRaDec.getY(), 3, "y coord failed for getMatchedSources()")
        self.gas.reset()
        

    def testVerifyG117(self):
        if verbose:
            print "testVerifyG117..."

        if not self.gas.exists():
            return

        crval = afwImage.PointD(141.063590, +35.280919)
        crpix = afwImage.PointD(446, 447)
        listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "g117.xy.txt")
        #To speed the test, tell the GAS what the size of the image is
        #The image is 1780 pixels on a side and covers half a square degree 
        #on the sky
        plateScale = .5*3600/1780.
        self.gas.setLogLevel(verbose)
        self.gas.reset()
        self.solveOrVerify(listFile, crval, crpix, plateScale, verify=True)
        self.gas.setLogLevel(0)
        self.gas.reset()
        
        
    def testVerifyCFHTField(self):
        if verbose:
            print "testVerifyCFHTField..."

        if not self.gas.exists():
            return
        
        crval = afwImage.PointD(334.303012, -17.233988)
        
        crpix = afwImage.PointD(512,512)
        listFile = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")
        #To speed the test, tell the GAS what the size of the image is
        #The image is 1780 pixels on a side and covers half a square degree 
        #on the sky
        plateScale = .185
        self.gas.reset()
        #self.gas.setParity(net.UNKNOWN_PARITY)
        self.gas.setLogLevel(verbose)
        self.solveOrVerify(listFile, crval, crpix, plateScale=plateScale, verify=True)
        self.gas.setLogLevel(0)


    if False:
        def testWcsSinglePixelOffset(self):
            if not self.gas.exists():
                return

            mastrom = eups.productDir("meas_astrom")
            imageFilename = "gd66.fits"
            starlist = os.path.join(mastrom, "tests", "gd66.xy.txt")


            filename = os.path.join(mastrom, "tests", "gd66.fits")

            starlist = loadXYFromFile(starlist)

            #Get Wcs from image header
            exposure = afwImg.ExposureF(filename)
            origWcs = exposure.getWcs()

            #Get Wcs from astrometry.net
            gasWcs = self.gas.solveWcs(starlist, origWcs)

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

class SmallSolveGASTestCFHT(unittest.TestCase):
    """A test case for WCS from astrometry.net using a small set of index files"""

    def setUp(self):
        self.gas = GAS("cfhttemplate")

    def tearDown(self):
        self.gas.reset()

    def solve(self, imgListFile, raDec, nBright=50, expectPass=True):
        starlist = loadXYFromFile(imgListFile)
        self.gas.setStarlist(starlist)
        self.gas.setNumBrightObjects(50)
    
        flag = self.gas.solve(raDec)
        if expectPass:
            self.assertTrue(flag, "No solution found")
            wcs = self.gas.getWcs()
            result = wcs.getOriginRaDec()
            scale= self.gas.getSolvedImageScale()
            if verbose:
                print "%.6f %.6f %.3f" %(result.getX(), result.getY(), scale)

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
        self.solve(starlist, crval)
        self.gas.setLogLevel(0)

        #
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
            return

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

    def testDistortedWcs(self):
        """Is a distorted Wcs returned"""

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


    def testSolveWcs(self):
        """Can't be tested until I figure out the LSST way of interfacing Eigen Matrices"""
        pass

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(WCSTestCaseNetUSNOB)
    suites += unittest.makeSuite(SmallSolveGASTestCFHT)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
