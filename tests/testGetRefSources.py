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
import eups

import numpy as np

import lsst.meas.astrom            as measAstrom
import lsst.meas.algorithms.utils  as measAlgUtil
import lsst.afw.detection          as afwDet
import lsst.afw.geom               as afwGeom
import lsst.afw.table              as afwTable
import lsst.afw.math               as afwMath
import lsst.afw.image              as afwImg
import lsst.utils.tests            as utilsTests
import lsst.pex.policy             as pexPolicy
from lsst.pex.logging import Log
import lsst.meas.photocal          as photocal

from lsst.pex.exceptions import LsstCppException

class GetRefSources(unittest.TestCase):

    def setUp(self):
        self.conf = measAstrom.MeasAstromConfig()
        mypath = eups.productDir("meas_astrom")
        # Set up local astrometry_net_data
        datapath = os.path.join(mypath, 'tests', 'astrometry_net_data', 'photocal')
        eupsObj = eups.Eups(root=datapath)
        ok, version, reason = eupsObj.setup('astrometry_net_data')
        if not ok:
            raise ValueError("Need photocal version of astrometry_net_data (from path: %s): %s" %
                             (datapath, reason))
        loglvl = Log.INFO
        self.astrom = measAstrom.Astrometry(self.conf, logLevel=loglvl)
        self.datapath = datapath

    def tearDown(self):
        del self.astrom
        import lsst.meas.astrom.astrometry_net as an
        an.finalize()

    def testGetSources(self):
        ra, dec, rad = (214.87 * afwGeom.degrees,
                        52.68 * afwGeom.degrees,
                        0.15 * afwGeom.degrees)

        cat = self.astrom.getReferenceSources(ra, dec, rad, "r")
        print 'cat', cat
        print len(cat)
        self.assertEqual(len(cat), 13)
        schema = cat.getSchema()
        print 'schema', schema
        # These must run without error.
        rkey = schema.find('r').key
        rekey = schema.find('r.err').key
        fkey = schema.find('flux').key
        fekey = schema.find('flux.err').key
        idkey = schema.find('id').key
        print '% 20s  % 20s  % 20s  % 20s  % 20s' % ('id', 'r flux', 'r flux err', 'r mag', 'r mag err')
        # flux == r
        for i,src in enumerate(cat):
            sid = src.get(idkey)
            r = src.get(rkey)
            f = src.get(fkey)
            self.assertEqual(r, f)
            re = src.get(rekey)
            fe = src.get(fekey)
            self.assertEqual(re, fe)
            rm = np.log10(r) / -0.4
            rme = np.abs(re / (-0.4 * r * np.log(10.)))
            print '% 20d  % 20.5g  % 20.5g  % 20.5f  % 20.5f' % (sid, r, re, rm, rme)

            if i == 0:
                self.assertLess(np.abs(r  - 1.6289e-10), 1e-13)
                self.assertLess(np.abs(re - 8.5042e-11), 1e-14)
                self.assertLess(np.abs(rm - 24.47027),   1e-4)
                self.assertLess(np.abs(rme- 0.56685),    1e-5)


        cat2 = self.astrom.getReferenceSources(ra, dec, rad, "g")
        self.assertEqual(len(cat2), len(cat))

        schema = cat2.getSchema()
        print 'schema', schema

        # These must run without error.
        bands = ['u','g','r','i','z', 'flux']
        for band in bands:
            fkey = schema.find(band).key
            fekey = schema.find(band + '.err').key
            print '% 10s  % 10s' % (band, band+'.err'),
        print
        nchecked = 0
        for i,src in enumerate(cat2):
            for band in bands:
                fkey = schema.find(band).key
                fekey = schema.find(band + '.err').key
                f = src.get(fkey)
                fe = src.get(fekey)
                m = np.log10(f) / -0.4
                me = np.abs(fe / (-0.4 * f * np.log(10.)))
                print '% 10.5f  % 10.5f' % (m, me),
                #print '% 10.5g  % 10.5g' % (f, fe),

                if i == 0:
                    # From Astrometry.net search-index command-line program
                    truemag = {'u': 24.1175, 'g': 24.0543, 'r': 24.4703,
                               'i': 22.1171, 'z': 21.1417 }
                    trueerr = {'u': 0.700183, 'g': 0.403793, 'r': 0.566848,
                               'i': 0.156623, 'z': 0.234956 }
                    truemag['flux'] = truemag['g']
                    trueerr['flux'] = trueerr['g']
                    if not band in truemag:
                        continue
                    self.assertLess(np.abs(m  - truemag[band]), 1e-4)
                    self.assertLess(np.abs(me - trueerr[band]), 1e-6)
                    nchecked += 1

            print
        self.assertEqual(nchecked, 6)

    def testNotAllFluxes(self):
        """Check that not all fluxes are returned when self.config.allFluxes=False.

        The default behavior (self.config.allFluxes=True) is checked by testGetSources.
        """
        ra, dec, rad = (214.87 * afwGeom.degrees,
                        52.68 * afwGeom.degrees,
                        0.15 * afwGeom.degrees)
        self.astrom.config.allFluxes = False
        cat = self.astrom.getReferenceSources(ra, dec, rad, "r")
        self.assertEqual(len(cat), 13)
        schema = cat.getSchema()

        # All these should pass without errors
        for name in ("r", "r.err"):
            schema.find(name)

        # All these should raise an exception
        for name in ("u", "g", "i", "z"):
            self.assertRaises(KeyError, schema.find, name)
            self.assertRaises(KeyError, schema.find, name + ".err")

    def testNoMagErrs(self):
        andconfig = measAstrom.AstrometryNetDataConfig()
        andconfig.load(os.path.join(self.datapath, 'andConfig2.py'))
        andconfig.magErrorColumnMap = {}
        astrom = measAstrom.Astrometry(self.conf, andConfig=andconfig)

        ra, dec, rad = (214.87 * afwGeom.degrees,
                        52.68 * afwGeom.degrees,
                        0.15 * afwGeom.degrees)

        cat = astrom.getReferenceSources(ra, dec, rad, "r")
        print 'cat', cat
        print len(cat)
        self.assertEqual(len(cat), 13)
        schema = cat.getSchema()
        print 'schema', schema

        for band in ['r', 'flux']:
            key = schema.find(band)
            with self.assertRaises(KeyError):
                ekey = schema.find(band + '.err')

        cat = astrom.getReferenceSources(ra, dec, rad, "r")
        print 'cat', cat
        print len(cat)
        self.assertEqual(len(cat), 13)
        schema = cat.getSchema()
        print 'schema', schema

        for band in ['u', 'g', 'r', 'i', 'z', 'flux']:
            key = schema.find(band)
            with self.assertRaises(KeyError):
                ekey = schema.find(band + '.err')

    def testRequestForeignFilter(self):
        """The user requests a filter not in the astrometry.net catalog.

        In that case, he's got a mapping in the MeasAstromConfig to point
        to an alternative filter he'll take instead (e.g., g instead of B).
        He should therefore expect the returned catalog to contain references
        to the bands that are in the catalog.
        """
        bands = ['u','g','r','i','z']
        andconfig = measAstrom.AstrometryNetDataConfig()
        andconfig.load(os.path.join(self.datapath, 'andConfig2.py'))
        self.conf.filterMap = dict(('my_'+b, b) for b in bands)
        astrom = measAstrom.Astrometry(self.conf, andConfig=andconfig)

        ra, dec, rad = (214.87 * afwGeom.degrees,
                        52.68 * afwGeom.degrees,
                        0.15 * afwGeom.degrees)

        cat = astrom.getReferenceSources(ra, dec, rad, "my_r")
        print 'cat', cat
        print len(cat)
        self.assertEqual(len(cat), 13)
        schema = cat.getSchema()
        print 'schema', schema

        for nm in bands + ['flux']:
            key = schema.find(nm)
            ekey = schema.find(nm + '.err')


    def testDifferentMagNames(self):
        """The astrometry.net catalog's magnitude columns are not named after filters.

        In that case, the AstrometryNetDataConfig has a mapping to point to the correct
        columns.  We should expect that the returned catalog refers to the bands
        requested (not the implementation-dependent column names).
        """
        andconfig = measAstrom.AstrometryNetDataConfig()
        andconfig.load(os.path.join(self.datapath, 'andConfig2.py'))
        bands = ['u','g','r','i','z']
        andconfig.magColumnMap = dict([('my_'+b, b) for b in bands])
        andconfig.magErrorColumnMap = dict([('my_'+b, b + "_err") for b in bands])
        astrom = measAstrom.Astrometry(self.conf, andConfig=andconfig)

        ra, dec, rad = (214.87 * afwGeom.degrees,
                        52.68 * afwGeom.degrees,
                        0.15 * afwGeom.degrees)

        cat = astrom.getReferenceSources(ra, dec, rad, "my_r")
        print 'cat', cat
        print len(cat)
        self.assertEqual(len(cat), 13)
        schema = cat.getSchema()
        print 'schema', schema

        for nm in ['my_' + b for b in bands] + ['flux']:
            key = schema.find(nm)
            ekey = schema.find(nm + '.err')


            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(GetRefSources)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
