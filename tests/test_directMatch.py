
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
import sys
import unittest
import glob

import numpy as np

import lsst.utils.tests
import lsst.meas.astrom
import lsst.geom

from lsst.meas.algorithms.testUtils import MockReferenceObjectLoaderFromFiles


RefCatDir = os.path.join(os.path.dirname(__file__), "data", "sdssrefcat")


class DirectMatchTestCase(lsst.utils.tests.TestCase):
    """Tests for lsst.meas.astrom.DirectMatchTask"""

    def setUp(self):
        np.random.seed(12345)
        filenames = sorted(glob.glob(os.path.join(RefCatDir, 'ref_cats', 'cal_ref_cat', '??????.fits')))
        self.refObjLoader = MockReferenceObjectLoaderFromFiles(filenames, htmLevel=8)
        center = lsst.geom.SpherePoint(215.5, 53.0, lsst.geom.degrees)
        radius = 0.5*lsst.geom.degrees
        self.filter = "r"
        self.references = self.refObjLoader.loadSkyCircle(center, radius, self.filter).refCat

    def tearDown(self):
        del self.refObjLoader
        del self.references

    def checkMatching(self, catalog, config=None):
        if config is None:
            config = lsst.meas.astrom.DirectMatchConfig()
        task = lsst.meas.astrom.DirectMatchTask(config=config, refObjLoader=self.refObjLoader)
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
        offset = (0.1*lsst.geom.arcseconds).asRadians()
        num = len(references)
        ra, dec = references["coord_ra"], references["coord_dec"]
        cosDec = np.cos(dec.mean())
        ra += offset/cosDec*np.random.uniform(-1.0, 1.0, num)
        dec += offset*np.random.uniform(-1.0, 1.0, num)
        self.checkMatching(references)

    def testNoSourceSelection(self):
        """Same results with source selector disabled."""
        config = lsst.meas.astrom.DirectMatchConfig()
        config.doSourceSelection = False
        self.checkMatching(self.references, config=config)


class DirectMatchMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules[__name__])
    unittest.main()
