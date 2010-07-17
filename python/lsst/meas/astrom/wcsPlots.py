import matplotlib
matplotlib.use('Agg')
from matplotlib.font_manager import FontProperties

from pylab import *
from numpy import array

#import lsst.afw.image.imageLib as afwImage
import lsst.afw.geom.geomLib as afwGeom
import lsst.afw.coord.coordLib as afwCoord

from astrometry.libkd import spherematch

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
    rx = array(rx)
    ry = array(ry)
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

    # All the id fields are zero, so I guess we have to do it the hard way...

    # NOTE that the reference source list here can contain duplicate
    # RA,Dec entries from each Astrometry.net index!
    # Also, ref sources can be *far* outside the image bounds.
    # --> not any more.

    # one milli-arcsec in degrees
    onemas = 1./(3600.*1000.)

    print '%i ref sources' % len(refsources)

    # Probably no longer necessary (and could be made *much* faster)...
    uniqrefsources = []
    for i,r1 in enumerate(refsources):
        c1 = afwCoord.Coord(r1.getRa(), r1.getDec(), 2000.0)
        duplicate = False
        for r2 in uniqrefsources:
            c2 = afwCoord.Coord(r2.getRa(), r2.getDec(), 2000.0)
            if c1.angularSeparation(c2, afwCoord.DEGREES) <= onemas:
                duplicate = True
                break
        if not duplicate:
            uniqrefsources.append(r1)
    print 'Trimmed reference sources from %i to %i\n' % (len(refsources), len(uniqrefsources))
    origrefsources = refsources
    refsources = uniqrefsources


    matchinds = []
    # match order is (cat,img).
    for m in matches:
        mcat = m.first
        mradec = afwCoord.Coord(mcat.getRa(), mcat.getDec(), 2000.0)
        cati = -1
        for i,s in enumerate(refsources):
            sradec = afwCoord.Coord(s.getRa(), s.getDec(), 2000.0)
            sep = mradec.angularSeparation(sradec, afwCoord.DEGREES)
            if sep < onemas:
                cati = i
                break
        mimg = m.second

        imgi = -1
        for i,s in enumerate(imgsources):
            sep = hypot(mimg.getXAstrom() - s.getXAstrom(),
                        mimg.getYAstrom() - s.getYAstrom())
            # 1 milli-pixel
            if sep < 1e-3:
                imgi = i
                break
        matchinds.append((cati, imgi))
    matchinds = array(matchinds)

    #print 'Match indices:', matchinds

    def flux2mag(f):
        return -2.5*log10(f)

    refmags = array([flux2mag(s.getPsfFlux()) for s in refsources])
    imgfluxes = array([s.getPsfFlux() for s in imgsources])

    # Matched mags:
    matchrefi = matchinds[:,0]
    matchimgi = matchinds[:,1]

    mimgflux = imgfluxes[matchimgi]
    okflux = (mimgflux > 1)
    mimgmag = flux2mag(mimgflux[okflux])
    mrefmag  = (refmags[matchrefi])[okflux]

    # unmatched:
    Uimg = ones(len(imgfluxes)).astype(bool)
    Uimg[matchimgi] = False
    Uref = ones(len(refmags)).astype(bool)
    Uimg[matchrefi] = False

    uimgflux = imgfluxes[Uimg]
    okflux = (uimgflux > 1)
    uimgmag = flux2mag(uimgflux[okflux])
    urefmag = refmags[Uref]

    clf()
    p1 = plot(mimgmag, mrefmag, 'b.')
    imag = append(mimgmag, uimgmag)
    axis([floor(min(imag))-0.5, ceil(max(imag)), floor(min(refmags))-0.5, ceil(max(refmags))])
    ax = axis()

    dy = (ax[3]-ax[2]) * 0.05
    y0 = ones_like(uimgmag) * ax[2]
    p2 = plot(vstack((uimgmag, uimgmag)), vstack((y0, y0+dy)), 'r-', alpha=0.5)
    p2 = p2[0]
    y0 = ones_like(mimgmag) * ax[2]
    p3 = plot(vstack((mimgmag, mimgmag)), vstack((y0+(0.25*dy), y0+(1.25*dy))), 'b-', alpha=0.5)
    p3 = p3[0]

    dx = (ax[1]-ax[0]) * 0.05
    x0 = ones_like(urefmag) * ax[0]
    p4 = plot(vstack((x0, x0+dx)), vstack((urefmag, urefmag)), 'r-', alpha=0.5)
    p4 = p4[0]
    x0 = ones_like(mrefmag) * ax[0]
    p5 = plot(vstack((x0+(0.25*dx), x0+(1.25*dx))), vstack((mrefmag, mrefmag)), 'b-', alpha=0.5)
    p5 = p5[0]

    figlegend((p1, p3, p2), ('Matched sources', 'Matched sources', 'Unmatched sources',),
              'center right', numpoints=1, prop=FontProperties(size='small'))

    axis(ax)
    xlabel('Image instrumental mag')
    ylabel('Reference catalog mag')
    
    fn = prefix + '-photom.png'
    print 'Saving', fn
    savefig(fn)


    # correspondences we could have hit...
    ixy = vstack((ix, iy)).T
    rxy = vstack((rx, ry)).T
    dcell = 50.
    radius = dcell * sqrt(2.)
    print 'ixy', ixy.shape
    print 'rxy', rxy.shape

    if False:
        (inds,dists) = spherematch.match(rxy, ixy, radius)
        mi = inds[:,0]
        ii = inds[:,1]
        matchx = rx[mi]
        matchy = ry[mi]
        matchdx = ix[ii] - matchx
        matchdy = iy[ii] - matchy
        ok = (matchdx >= -dcell) * (matchdx <= dcell) * (matchdy >= -dcell) * (matchdy <= dcell)
        matchx = matchx[ok]
        matchy = matchy[ok]
        matchdx = matchdx[ok]
        matchdy = matchdy[ok]
        mi = mi[ok]
        ii = ii[ok]
        print 'Found %i matches within %g pixels' % (len(dists), radius)

    ncells = 18.
    cellsize = sqrt(W * H / ncells)
    nw = int(round(W / cellsize))
    nh = int(round(H / cellsize))
    print 'Grid cell size', cellsize
    print 'N cells', nw, 'x', nh
    edgesx = linspace(0, W, nw+1)
    edgesy = linspace(0, H, nh+1)

    binx = digitize(rx, edgesx)
    biny = digitize(ry, edgesy)
    binx = clip(binx - 1, 0, nw-1)
    biny = clip(biny - 1, 0, nh-1)

    bin = biny * nw + binx
    
    clf()

    for i in range(nh):
        for j in range(nw):
            thisbin = i * nw + j
            R = (bin == thisbin)
            print 'cell %i, %i' % (j, i)
            print '%i ref sources' % sum(R)
            if sum(R) == 0:
                continue
            (inds,dists) = spherematch.match(rxy[R,:], ixy, radius)
            print 'Found %i matches within %g pixels' % (len(dists), radius)
            ri = inds[:,0]
            # un-cut ref inds...
            ri = (flatnonzero(R))[ri]
            ii = inds[:,1]

            matchx  = rx[ri]
            matchy  = ry[ri]
            matchdx = ix[ii] - matchx
            matchdy = iy[ii] - matchy
            ok = (matchdx >= -dcell) * (matchdx <= dcell) * (matchdy >= -dcell) * (matchdy <= dcell)
            #matchx = matchx[ok]
            #matchy = matchy[ok]
            matchdx = matchdx[ok]
            matchdy = matchdy[ok]
            print 'Cut to %i within %g x %g square' % (sum(ok), dcell*2, dcell*2)

            # Subplot places plots left-to-right, TOP-to-BOTTOM.
            #subplot(nh, nw, thisbin+1)
            subplot(nh, nw, 1 + ((nh - i - 1)*nw + j))

            plot(matchdx, matchdy, 'ro', mec='r', mfc='r', ms=5, alpha=0.2)
            plot(matchdx, matchdy, 'ro', mec='r', mfc='none', ms=5, alpha=0.2)
            #plot([0], [0], 'bo', mec='b', mfc='none')
            axhline(0, color='k', alpha=0.5)
            axvline(0, color='k', alpha=0.5)
            xticks([],[])
            yticks([],[])
            axis('scaled')
            axis([-dcell, dcell, -dcell, dcell])


    fn = prefix + '-missed.png'
    print 'Saving', fn
    savefig(fn)

    
def plotDistortion(sip, W, H, ncells, prefix, title, exaggerate=1.):
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

