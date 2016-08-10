#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

#
# LSST Data Management System
# Copyright 2008-2016 LSST Corporation.
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
import numpy

import lsst.utils.tests
import lsst.meas.astrom
import lsst.afw.geom

from lsst.utils import getPackageDir
from lsst.daf.persistence import Butler
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask

RefCatDir = os.path.join(getPackageDir("meas_astrom"), "tests", "data", "sdssrefcat")

numpy.random.seed(12345)


class MatchTestCase(lsst.utils.tests.TestCase):
    """Tests for lsst.meas.astrom.MatchTask"""

    def setUp(self):
        self.butler = Butler(RefCatDir)
        refObjLoader = LoadIndexedReferenceObjectsTask(butler=self.butler)
        center = lsst.afw.coord.IcrsCoord(215.5*lsst.afw.geom.degrees, 53.0*lsst.afw.geom.degrees)
        radius = 0.5*lsst.afw.geom.degrees
        self.filter = "r"
        self.references = refObjLoader.loadSkyCircle(center, radius, self.filter).refCat

    def tearDown(self):
        del self.butler
        del self.references

    def checkMatching(self, catalog):
        config = lsst.meas.astrom.MatchConfig()
        config.refObjLoader.retarget(LoadIndexedReferenceObjectsTask)
        task = lsst.meas.astrom.MatchTask(config=config, butler=self.butler)
        results = task.run(catalog, self.filter)

        self.assertEqual(len(results.matches), len(catalog))
        for match in results.matches:
            self.assertEqual(match.first.getId(), match.second.getId())
        maxDistance = max(match.distance for match in results.matches)
        self.assertLess(maxDistance, config.matchRadius)  # match.distance is in arcsec

        self.assertIsNotNone(results.matchMeta)
        names = results.matchMeta.names()
        for key in ("RA", "DEC", "RADIUS", "SMATCHV", "FILTER"):
            self.assertIn(key, names)

    def testWithoutNoise(self):
        """Match the reference catalog against itself"""
        self.checkMatching(self.references)

    def testWithNoise(self):
        """Match the reference catalog against a noised version of itself"""
        references = self.references.copy(True)
        offset = (0.1*lsst.afw.geom.arcseconds).asRadians()
        num = len(references)
        ra, dec = references["coord_ra"], references["coord_dec"]
        cosDec = numpy.cos(dec.mean())
        ra += offset/cosDec*numpy.random.uniform(-1.0, 1.0, num)
        dec += offset*numpy.random.uniform(-1.0, 1.0, num)
        self.checkMatching(references)


class MatchMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
