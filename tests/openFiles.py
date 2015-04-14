#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2013 Dustin Lang.
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
import os
import unittest

import resource

import lsst.meas.astrom as measAstrom
import lsst.utils.tests as utilsTests
import lsst.afw.table as afwTable


# http://stackoverflow.com/a/7142094/834250
def getOpenFiles():
    '''
    return the filenames for open file descriptors for current process

    .. warning: will only work on UNIX-like os-es.
    '''
    import subprocess
    import os

    pid = os.getpid()
    output = subprocess.check_output(["/usr/sbin/lsof", '-w', '-Ffn', "-p", str(pid)])

    procs = [out for out in output.split('\n') if out and out[0] in "fn"]
    assert len(procs) % 2 == 0 # Expect an even number of lines: f and n pairs
    procs = [(procs[2*i], procs[2*i+1]) for i in range(len(procs)//2)]
    nameList = [name[1:] for desc,name in procs if desc[0] == 'f' and desc[1:].isdigit()]
    return nameList

def printOpenFiles():
    names = getOpenFiles()
    print "Open files (%d): %s" % (len(names), names)


class OpenFilesTest(unittest.TestCase):
    """This tests that the astrometry functions can run with a greatly reduced open file limit

    There is no specific assert; we're just testing that things can run.  If they can't, then
    the code will fail (often catastrophically with a segfault), so there's no need to check
    specific success conditions.
    """

    def setUp(self):
        limits = resource.getrlimit(resource.RLIMIT_NOFILE)
        print 'NOFILE rlimit:', limits
        resource.setrlimit(resource.RLIMIT_NOFILE, (10, limits[1]))
        print 'NOFILE rlimit:', resource.getrlimit(resource.RLIMIT_NOFILE)

        self.mypath = os.path.dirname(os.path.dirname(__file__))
        path = os.path.join(self.mypath, "examples")
        self.srcCat = afwTable.SourceCatalog.readFits(os.path.join(path, "v695833-e0-c000.xy.fits"))
        self.srcCat.table.defineApFlux("flux.psf")
        # The .xy.fits file has sources in the range ~ [0,2000],[0,4500]
        self.imageSize = (2048, 4612) # approximate

    def getAstrom(self):
        andpath = os.path.join(self.mypath, 'tests', 'astrometry_net_data', 'photocal')
        os.environ['ASTROMETRY_NET_DATA_DIR'] = andpath
        andcfn = os.path.join(andpath, 'andConfigOpenFiles.py')

        andconfig = measAstrom.AstrometryNetDataConfig()
        andconfig.load(andcfn)

        conf = measAstrom.ANetBasicAstrometryConfig()
        return measAstrom.ANetBasicAstrometryTask(config=conf, andConfig=andconfig,)
                                            #logLevel=pexLog.Log.DEBUG)

    def tearDown(self):
        del self.imageSize
        del self.srcCat


    def runDetermineWcs(self):
        astrom = self.getAstrom()
        result = astrom.determineWcs2(self.srcCat, imageSize=self.imageSize, filterName='i')
        print 'Got result from determineWcs:', result
        #printOpenFiles()
        return result.wcs

    def runUseKnownWcs(self, wcs):
        astrom = self.getAstrom()
        result = astrom.useKnownWcs(self.srcCat, wcs=wcs, filterName='i', imageSize=self.imageSize)
        print "Got result from useKnownWcs:", result
        #printOpenFiles()

    def testDetermineWcs(self):
        self.runDetermineWcs()
        self.runDetermineWcs()
        self.runDetermineWcs()

    def testUseKnownWcs(self):
        wcs = self.runDetermineWcs()
        self.runUseKnownWcs(wcs)
        self.runUseKnownWcs(wcs)
        self.runUseKnownWcs(wcs)

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(OpenFilesTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
    


