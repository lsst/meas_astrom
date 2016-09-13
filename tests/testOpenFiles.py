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
from __future__ import print_function
import os
import unittest
import resource

import lsst.meas.astrom as measAstrom
import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable

from testFindAstrometryNetDataDir import setupAstrometryNetDataDir


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
    assert len(procs) % 2 == 0  # Expect an even number of lines: f and n pairs
    procs = [(procs[2*i], procs[2*i+1]) for i in range(len(procs)//2)]
    nameList = [name[1:] for desc, name in procs if desc[0] == 'f' and desc[1:].isdigit()]
    return nameList


def printOpenFiles():
    names = getOpenFiles()
    print("Open files (%d): %s" % (len(names), names))


class OpenFilesTest(unittest.TestCase):
    """This tests that the astrometry functions can run with a greatly reduced open file limit

    There is no specific assert; we're just testing that things can run.  If they can't, then
    the code will fail (often catastrophically with a segfault), so there's no need to check
    specific success conditions.
    """

    def setUp(self):
        self.originalLimits = resource.getrlimit(resource.RLIMIT_NOFILE)
        print('NOFILE rlimit:', self.originalLimits)
        resource.setrlimit(resource.RLIMIT_NOFILE, (10, self.originalLimits[1]))
        print('NOFILE rlimit:', resource.getrlimit(resource.RLIMIT_NOFILE))

        mypath = os.path.dirname(__file__)
        self.andpath = setupAstrometryNetDataDir('photocal', rootDir=mypath)
        self.srcCat = afwTable.SourceCatalog.readFits(
            os.path.join(mypath, "v695833-e0-c000.xy.fits"))
        # The .xy.fits file has sources in the range ~ [0,2000],[0,4500]
        self.bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(2048, 4612))  # approximate

    def tearDown(self):
        resource.setrlimit(resource.RLIMIT_NOFILE, (self.originalLimits[0], self.originalLimits[1]))
        del self.bbox
        del self.srcCat

    def getAstrom(self):
        andcfn = os.path.join(self.andpath, 'andConfigOpenFiles.py')

        andconfig = measAstrom.AstrometryNetDataConfig()
        andconfig.load(andcfn)

        conf = measAstrom.ANetBasicAstrometryConfig()
        return measAstrom.ANetBasicAstrometryTask(config=conf, andConfig=andconfig,)

    def runDetermineWcs(self):
        astrom = self.getAstrom()
        result = astrom.determineWcs2(self.srcCat, bbox=self.bbox, filterName='i')
        print('Got result from determineWcs:', result)
        # printOpenFiles()
        return result.wcs

    def runUseKnownWcs(self, wcs):
        astrom = self.getAstrom()
        result = astrom.useKnownWcs(self.srcCat, wcs=wcs, filterName='i', bbox=self.bbox)
        print("Got result from useKnownWcs:", result)
        # printOpenFiles()

    def testDetermineWcs(self):
        self.runDetermineWcs()
        self.runDetermineWcs()
        self.runDetermineWcs()

    def testUseKnownWcs(self):
        wcs = self.runDetermineWcs()
        self.runUseKnownWcs(wcs)
        self.runUseKnownWcs(wcs)
        self.runUseKnownWcs(wcs)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
