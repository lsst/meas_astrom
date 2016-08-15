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

import os
import unittest

from lsst.afw.coord import IcrsCoord
import lsst.meas.astrom as measAstrom
import lsst.utils.tests
import lsst.afw.geom as afwGeom
from lsst.pex.logging import Log
import testFindAstrometryNetDataDir as helper


class MultipleCatalogStarsTest(unittest.TestCase):

    def setUp(self):
        # Set up local astrometry_net_data

        datapath = helper.setupAstrometryNetDataDir('photocal')
        self.conf = measAstrom.ANetBasicAstrometryConfig()
        # Load andConfig2.py rather than the default.
        confpath = os.path.join(datapath, 'andConfig2.py')
        self.andConfig = measAstrom.AstrometryNetDataConfig()
        self.andConfig.load(confpath)

    def tearDown(self):
        del self.conf
        del self.andConfig

    def testGetCatalog(self, logLevel=Log.DEBUG):
        astrom = measAstrom.ANetBasicAstrometryTask(self.conf, andConfig=self.andConfig)
        astrom.log.setThreshold(logLevel)

        ctrCoord = IcrsCoord(
            215.6 * afwGeom.degrees,
            53.0 * afwGeom.degrees,
        )
        cat = astrom.refObjLoader.loadSkyCircle(
            ctrCoord=ctrCoord,
            radius=0.1 * afwGeom.degrees,
            filterName='z',
        ).refCat
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
            self.assertNotIn(sid, ids)
            ids.add(sid)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(MultipleCatalogStarsTest)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)

    return unittest.TestSuite(suites)


def run(exit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
