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
import sys
import glob
import math
import unittest

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

import eups
import lsst.meas.astrom            as measAstrom
import lsst.utils.tests            as utilsTests
import lsst.afw.geom as afwGeom
from lsst.pex.logging import Log

class MultipleCatalogStarsTest(unittest.TestCase):

    def setUp(self):
        # Set up local astrometry_net_data
        mypath = eups.productDir("meas_astrom")
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))

        self.conf = measAstrom.MeasAstromConfig()
        # Load andConfig2.py rather than the default.
        confpath = os.path.join(datapath, 'andConfig2.py')
        self.andconf = measAstrom.AstrometryNetDataConfig()
        self.andconf.load(confpath)

    def tearDown(self):
        del self.conf
        import lsst.meas.astrom.astrometry_net as an
        an.finalize()

    def testGetCatalog(self, loglvl=Log.DEBUG):
        astrom = measAstrom.Astrometry(self.conf, andConfig=self.andconf, logLevel=loglvl)

        cat = astrom.getReferenceSources(215.6 * afwGeom.degrees,
                                         53.0 * afwGeom.degrees,
                                         0.1 * afwGeom.degrees,
                                         'z')
        print 'Got', len(cat), 'reference sources'

        ids = set(s.getId() for s in cat)
        print len(ids), 'unique IDs'
        ras = set(s.getRa() for s in cat)
        print len(ras), 'unique RAs'
        

        ids = set()
        for src in cat:
            sid = src.getId()
            if sid in ids:
                print 'Source id', sid, 'is duplicated'
            self.assertFalse(sid in ids)
            ids.add(sid)
            
            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(MultipleCatalogStarsTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
