#! /usr/bin/env python

import unittest
import lsst.utils.tests as utilsTests

import lsst.meas.astrom.sip as astromSip
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImg
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord

class removeOneToManyTest(unittest.TestCase):

	def test1(self):
		crval = afwGeom.Point2D()
		crpix = afwGeom.Point2D()
		pixscale = 1./3600 # deg/pix
		wcs = afwImg.createWcs(crval, crpix, pixscale, 0, 0, pixscale)

		sid = 1

		cat = afwDet.SourceSet()
		rd = wcs.pixelToSky(0,0)
		r,d = rd.getLongitude(afwCoord.DEGREES),rd.getLatitude(afwCoord.DEGREES)
		for i in range(3):
			s = afwDet.Source()
			s.setRa(r), s.setDec(d)
			s.setId(sid)
			sid += 1
			cat.append(s)
		rd = wcs.pixelToSky(50,0)
		r,d = rd.getLongitude(afwCoord.DEGREES),rd.getLatitude(afwCoord.DEGREES)
		for i in range(3):
			s = afwDet.Source()
			s.setRa(r), s.setDec(d)
			s.setId(sid)
			sid += 1
			cat.append(s)
		rd = wcs.pixelToSky(0.5,0)
		r,d = rd.getLongitude(afwCoord.DEGREES),rd.getLatitude(afwCoord.DEGREES)
		for i in range(3):
			s = afwDet.Source()
			s.setRa(r), s.setDec(d)
			s.setId(sid)
			sid += 1
			cat.append(s)
		rd = wcs.pixelToSky(0,0)
		r,d = rd.getLongitude(afwCoord.DEGREES),rd.getLatitude(afwCoord.DEGREES)
		for i in range(3):
			s = afwDet.Source()
			s.setRa(r), s.setDec(d)
			s.setId(sid)
			sid += 1
			cat.append(s)

		img = afwDet.SourceSet()
		for i in range(3):
			s = afwDet.Source()
			s.setXAstrom(0), s.setYAstrom(0)
			s.setId(sid)
			sid += 1
			img.append(s)
		for i in range(3):
			s = afwDet.Source()
			s.setXAstrom(50), s.setYAstrom(0)
			s.setId(sid)
			sid += 1
			img.append(s)
		for i in range(3):
			s = afwDet.Source()
			s.setXAstrom(0), s.setYAstrom(0)
			s.setId(sid)
			sid += 1
			img.append(s)

		dist = 1. # arcsec

		m = astromSip.MatchSrcToCatalogue(cat, img, wcs, dist)
		matches = m.getMatches()
		#print 'Got matches:', matches
		print 'Matches:'
		for m in matches:
			print '  ', m.first, '--', m.second
		self.assertEqual(len(matches), 2)

# ----------- ridiculous boilerplate:

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(removeOneToManyTest)
    return unittest.TestSuite(suites)

def run(exit=False):
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
