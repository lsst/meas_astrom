import matplotlib
matplotlib.use('Agg')
from matplotlib.font_manager import FontProperties

from pylab import *
from numpy import array

import lsst.afw.geom.geomLib as afwGeom
import lsst.afw.coord.coordLib as afwCoord

from astrometry.libkd import spherematch

def _getplotdata(format='png'):
    import cStringIO
    io = cStringIO.StringIO()
    savefig(io, format=format)
    val = io.getvalue()
    io.close()
    return val

def _output(fn, format, write):
    if write:
        savefig(fn)
    else:
        return {fn: _getplotdata(format)}

def plotMatches(imgsources, refsources, matches, wcs, W, H, prefix,
                saveplot=True, format='png'):
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
    #print p1, p2, p3, p4, p5

    figlegend((p1, p2, p3, p4), #, p5),
              ('Image sources (brightest 200)',
               'Image sources (rest)',
               'Reference sources',
               'Matches',),
              'center right',
              numpoints=1,
              prop=FontProperties(size='small'))

    fn = prefix + '-matches.' + format
    return _output(fn, format, saveplot)

def plotPhotometry(imgsources, refsources, matches, prefix, band=None,
                   zp=None,
                   saveplot=True, format='png'):
    print '%i ref sources' % len(refsources)
    print '%i image sources' % len(imgsources)
    print '%i matches' % len(matches)

    # In the "matches" list:
    #    m.first  is catalog
    #    m.second is image

    # In this function, the "m" prefix stands for "matched",
    # "u" stands for "unmatched".
    # "i" indicates indices into an array, so "mrefi" are matched
    # reference star indices.

    # *sigh*, turn these into Python lists, so we have the "index" function.
    refsources = [s for s in refsources]
    imgsources = [s for s in imgsources]
    
    # These are the indices in the "ref/imgsources" arrays of the matched stars.
    mrefi = array([refsources.index(m.first ) for m in matches])
    mimgi = array([imgsources.index(m.second) for m in matches])

    def flux2mag(f):
        return -2.5*log10(f)

    refmag = array([flux2mag(s.getPsfFlux()) for s in refsources])
    imgflux = array([s.getPsfFlux() for s in imgsources])

    # Cut to fluxes that aren't silly and get mags of matched sources.
    mimgflux = imgflux[mimgi]
    okflux = (mimgflux > 1)
    mimgmag = flux2mag(mimgflux[okflux])
    mrefmag  = (refmag[mrefi])[okflux]

    # The indices of unmatched stars...
    uimg = ones(len(imgflux), bool)
    uimg[mimgi] = False
    uref = ones(len(refmag), bool)
    uref[mrefi] = False

    # Get mags of unmatched sources.
    uimgflux = imgflux[uimg]
    okflux = (uimgflux > 1)
    uimgmag = flux2mag(uimgflux[okflux])
    urefmag = refmag[uref]

    clf()
    p1 = plot(mimgmag, mrefmag, 'b.', alpha=0.5)
    imag = append(mimgmag, uimgmag)
    axis([floor(min(imag))-0.5, ceil(max(imag)), floor(min(refmag))-0.5, ceil(max(refmag))])
    ax = axis()

    # Red tick marks show unmatched img sources
    dy = (ax[3]-ax[2]) * 0.05
    y1 = ones_like(uimgmag) * ax[3]
    p2 = plot(vstack((uimgmag, uimgmag)), vstack((y1, y1-dy)), 'r-', alpha=0.5)
    p2 = p2[0]
    # Blue tick marks show matched img sources
    y1 = ones_like(mimgmag) * ax[3]
    p3 = plot(vstack((mimgmag, mimgmag)), vstack((y1-(0.25*dy), y1-(1.25*dy))), 'b-', alpha=0.5)
    p3 = p3[0]

    # Red ticks for unmatched ref sources
    dx = (ax[1]-ax[0]) * 0.05
    x1 = ones_like(urefmag) * ax[1]
    p4 = plot(vstack((x1, x1-dx)), vstack((urefmag, urefmag)), 'r-', alpha=0.5)
    p4 = p4[0]
    # Blue ticks for matched ref sources
    x1 = ones_like(mrefmag) * ax[1]
    p5 = plot(vstack((x1-(0.25*dx), x1-(1.25*dx))), vstack((mrefmag, mrefmag)), 'b-', alpha=0.5)
    p5 = p5[0]

    # Legend entries:
    pp = [p1, p3, p2]
    pl = ['Matched sources', 'Matched sources', 'Unmatched sources']

    if zp is not None:
        X = array([ax[0], ax[1]])
        pzp = plot(X, X+zp, 'b--')
        pp.append(pzp)
        pl.append('Zeropoint')

    figlegend(pp, pl, 'center right', numpoints=1, prop=FontProperties(size='small'))

    # reverse axis directions.
    axis([ax[1],ax[0], ax[3], ax[2]])
    xlabel('Image instrumental mag')
    if band is not None:
        ylabel('Reference catalog: %s band (mag)' % band)
    else:
        ylabel('Reference catalog mag')

    fn = prefix + '-photom.' + format
    return _output(fn, format, saveplot)

def plotCorrespondences2(imgsources, refsources, matches, wcs, W, H, prefix,
                         saveplot=True, format='png'):
    from astrometry.util.plotshift import plotshift

    ix = array([s.getXAstrom() for s in imgsources])
    iy = array([s.getYAstrom() for s in imgsources])

    rx,ry = [],[]
    for r in refsources:
        xy = wcs.skyToPixel(r.getRa(), r.getDec())
        rx.append(xy[0])
        ry.append(xy[1])
    rx = array(rx)
    ry = array(ry)

    ixy = vstack((ix, iy)).T
    rxy = vstack((rx, ry)).T

    cell = 10
    plotshift(ixy, rxy, dcell=cell, ncells=9, W=W, H=H)
    fn = prefix + '-shift1.' + format
    P1 = _output(fn, format, saveplot)

    clf()
    hot()
    plotshift(ixy, rxy, dcell=cell, ncells=9, W=W, H=H, hist=True, nhistbins=2*cell+1)
    fn = prefix + '-shift2.' + format
    P2 = _output(fn, format, saveplot)

    cell = 2
    plotshift(ixy, rxy, dcell=cell, ncells=9, W=W, H=H)
    fn = prefix + '-shift3.' + format
    P3 = _output(fn, format, saveplot)

    clf()
    hot()
    plotshift(ixy, rxy, dcell=cell, ncells=9, W=W, H=H, hist=True, nhistbins=10*cell+1)
    fn = prefix + '-shift4.' + format
    P4 = _output(fn, format, saveplot)

    if not saveplot:
        P1.update(P2)
        P1.update(P3)
        P1.update(P4)
        return P1


def plotCorrespondences(imgsources, refsources, matches, wcs, W, H, prefix):
    ix = array([s.getXAstrom() for s in imgsources])
    iy = array([s.getYAstrom() for s in imgsources])

    rx,ry = [],[]
    for r in refsources:
        xy = wcs.skyToPixel(r.getRa(), r.getDec())
        rx.append(xy[0])
        ry.append(xy[1])
    rx = array(rx)
    ry = array(ry)

    # correspondences we could have hit...
    ixy = vstack((ix, iy)).T
    rxy = vstack((rx, ry)).T
    dcell = 50.
    radius = dcell * sqrt(2.)
    #print 'ixy', ixy.shape
    #print 'rxy', rxy.shape

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
    #print 'Grid cell size', cellsize
    #print 'N cells', nw, 'x', nh
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
            #print 'cell %i, %i' % (j, i)
            #print '%i ref sources' % sum(R)
            if sum(R) == 0:
                continue
            (inds,dists) = spherematch.match(rxy[R,:], ixy, radius)
            #print 'Found %i matches within %g pixels' % (len(dists), radius)
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
            #print 'Cut to %i within %g x %g square' % (sum(ok), dcell*2, dcell*2)

            # Subplot places plots left-to-right, TOP-to-BOTTOM.
            subplot(nh, nw, 1 + ((nh - i - 1)*nw + j))

            plot(matchdx, matchdy, 'ro', mec='r', mfc='r', ms=5, alpha=0.2)
            plot(matchdx, matchdy, 'ro', mec='r', mfc='none', ms=5, alpha=0.2)
            axhline(0, color='k', alpha=0.5)
            axvline(0, color='k', alpha=0.5)
            xticks([],[])
            yticks([],[])
            axis('scaled')
            axis([-dcell, dcell, -dcell, dcell])


    fn = prefix + '-missed.png'
    print 'Saving', fn
    savefig(fn)


def wcsPlots(wcs, imgsources, refsources, matches, W, H, prefix, titleprefix,
             plotdata=None, plotformat='png'):
    '''Create diagnostic plots for WCS determination.

    wcs -- an lsst.afw.image.Wcs
    imgsources -- an lsst.afw.detection.SourceSet, of sources found in the image.
    refsources -- an lsst.afw.detection.SourceSet, of sources in the reference catalog.
    matches -- an lsst.afw.detection.SourceMatchSet (vector of SourceMatch)
    W, H -- ints, the image width and height
    prefix -- the output filename prefix for the plots.
    '''

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
    print 'WCS plots'
    D = plotMatches(imgsources, refsources, matches, wcs, W, H, prefix,
                    saveplot=(plotdata is None), format=plotformat)
    if plotdata is not None:
        plotdata.update(D)
    D = plotPhotometry(imgsources, refsources, matches, prefix,
                       saveplot=(plotdata is None), format=plotformat)
    if plotdata is not None:
        plotdata.update(D)
    #plotCorrespondences(imgsources, refsources, matches, wcs, W, H, prefix)
    D = plotCorrespondences2(imgsources, refsources, matches, wcs, W, H, prefix,
                             saveplot=(plotdata is None), format=plotformat)
    if plotdata is not None:
        plotdata.update(D)

    
def plotDistortion(sip, W, H, ncells, prefix, titletxt, exaggerate=1.,
                   saveplot=True, format='png', suffix='-distort.'):
    '''
    Produces a plot showing the SIP distortion that was found, by drawing
    a grid and distorting it.  Allows exaggeration of the distortion for ease
    of visualization.

    sip -- an lsst.afw.image.TanWcs
    W, H -- the image size
    ncells -- the approximate number of grid cells to split the image into.
    prefix -- output plot filename prefix.
    exaggerate -- the factor by which to exaggerate the distortion.
    
    '''
    ncells = float(ncells)
    cellsize = sqrt(W * H / ncells)
    nw = int(floor(W / cellsize))
    nh = int(floor(H / cellsize))
    #print 'Grid cell size', cellsize
    #print 'N cells', nw, 'x', nh
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

    title(titletxt)

    fn = prefix + suffix + format
    return _output(fn, format, saveplot)

