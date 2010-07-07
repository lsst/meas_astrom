import matplotlib
matplotlib.use('Agg')
from matplotlib.font_manager import FontProperties

from pylab import *
from numpy import array

#import lsst.afw.image.imageLib as afwImage
import lsst.afw.geom.geomLib as afwGeom

'''
matches
----------------------
Produced by afw::detection::matchRaDec(_catSet, _imgSet, _distInArcsec);
in afw/src/detection/SourceMatch.cc
is a vector<det::SourceMatch>
SourceMatch: Source first, Source second, double distance
Source.{getXAstrom(), getYAstrom, getRaObject(), getDecObject(), etc}

refsources
----------------------
Produced by solver.getCatalogue()
is an afw::detection::SourceSet
Source.getRa(), getDec()

--- both img sources and ref sources should have x,y,ra,dec set
 by meas_astrom:sip:MatchSrcToCatalogue.cc : findMatches()

'''
def wcsPlots(wcs, imgsources, refsources, matches, W, H, prefix, titleprefix):
	print 'WCS plots'

	clf()

	# Image sources
	ix = array([s.getXAstrom() for s in imgsources])
	iy = array([s.getYAstrom() for s in imgsources])
	iflux = array([s.getPsfFlux() for s in imgsources])
	I = argsort(-iflux)
	# First 200: red dots
	II = I[:200]
	p1 = plot(ix[II], iy[II], 'r.', zorder=10)
	# Rest: tiny dots
	II = I[200:]
	p2 = plot(ix[II], iy[II], 'r.', markersize=1, zorder=9)

	# Ref sources:
	# Only getRa() (not getRaAstrom(), getRaObject()) is non-zero.
	#print r0.getRa(), r0.getRaAstrom(), r0.getRaObject()

	rx,ry = [],[]
	for r in refsources:
		xy = wcs.skyToPixel(r.getRa(), r.getDec())
		rx.append(xy[0])
		ry.append(xy[1])
	p3 = plot(rx, ry, 'bo', mec='b', mfc='none', markersize=6, zorder=20)

	x,y = [],[]
	dx,dy = [],[]
	for m in matches:
		x0,x1 = m.first.getXAstrom(), m.second.getXAstrom()
		y0,y1 = m.first.getYAstrom(), m.second.getYAstrom()
		#plot([x0, x1], [y0, y1], 'g.-')
		x.append(x0)
		y.append(y0)
		dx.append(x1-x0)
		dy.append(y1-y0)
	#plot(x, y, 's', mec='g', mfc='none', markersize=5)
	p4 = plot(x, y, 'o', mec='g', mfc='g', alpha=0.5, markersize=8, zorder=5)
	p5 = quiver(x, y, dx, dy, angles='xy', scale=30., zorder=30)

	axis('scaled')
	axis([0, W, 0, H])

	print p1, p2, p3, p4, p5

	figlegend((p1, p2, p3, p4), #, p5),
			  ('Image sources (brightest 200)',
			   'Image sources (rest)',
			   'Reference sources',
			   'Matches',),
			   #'Match difference'),
			  #'upper right')
			  'center right',
			  numpoints=1,
			  prop=FontProperties(size='small'))
	  
	fn = prefix + '-matches.png'
	print 'Saving', fn
	savefig(fn)
	
def plotDistortion(sip, W, H, ncells, prefix, title, exaggerate=1.):
	print 'SIP:', sip
	pix = afwGeom.makePointD(W/2, H/2)
	print 'Pixel position', pix
	distpix = sip.distortPixel(pix)
	print 'Distorted:', distpix
	dx,dy = distpix[0]-pix[0], distpix[1]-pix[1]
	print 'dx,dy', dx, dy

	ncells = float(ncells)
	cellsize = sqrt(W * H / ncells)
	nw = int(floor(W / cellsize))
	nh = int(floor(H / cellsize))
	print 'Grid cell size', cellsize
	print 'N cells', nw, 'x', nh
	cx = arange(nw+1) * cellsize + ((W - (nw*cellsize))/2.)
	cy = arange(nh+1) * cellsize + ((H - (nh*cellsize))/2.)

	# pixel step size for grid lines
	step = 50

	xx = arange(-step, W+2*step, step)
	yy = arange(-step, H+2*step, step)

	clf()

	for y in cy:
		dx,dy = [],[]
		for x in xx:
			pix = afwGeom.makePointD(x, y)
			distpix = sip.distortPixel(pix)
			dx.append(distpix[0])
			dy.append(distpix[1])
		plot(xx, y*ones_like(xx), 'k-', zorder=10)
		dx = array(dx)
		dy = array(dy)
		if exaggerate != 1:
			dx += (exaggerate * (dx - xx))
			dy += (exaggerate * (dy - y))
		plot(dx, dy, 'r-', zorder=20)

	for x in cx:
		dx,dy = [],[]
		for y in yy:
			pix = afwGeom.makePointD(x, y)
			distpix = sip.distortPixel(pix)
			dx.append(distpix[0])
			dy.append(distpix[1])
		plot(x*ones_like(yy), yy, 'k-', zorder=10)
		dx = array(dx)
		dy = array(dy)
		if exaggerate != 1:
			dx += (exaggerate * (dx - x))
			dy += (exaggerate * (dy - yy))
		plot(dx, dy, 'r-', zorder=20)

	
	axis('scaled')
	axis([0, W, 0, H])

	fn = prefix + '-distort.png'
	print 'Saving', fn
	savefig(fn)
