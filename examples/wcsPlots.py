from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import range
import matplotlib
matplotlib.use('Agg')
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Ellipse
import numpy as np
import pylab as plt

import lsst.afw.geom.geomLib as afwGeom
import lsst.afw.coord.coordLib as afwCoord

# These only exist in recent Astrometry.net versions...
#from astrometry.libkd import spherematch
#from astrometry.util.plotshift import plotshift
# dstn copied a version into this dir until we uprev...
#from plotshift import plotshift


def _getplotdata(format='png'):
    import io
    io = io.StringIO()
    plt.savefig(io, format=format)
    val = io.getvalue()
    io.close()
    return val


def _output(fn, format, write):
    if write:
        plt.savefig(fn)
    else:
        return {fn: _getplotdata(format)}


def plotMatches(imgsources, refsources, matches, wcs, W, H, prefix,
                saveplot=True, format='png'):
    plt.clf()

    # Image sources
    ix = np.array([s.getXAstrom() for s in imgsources])
    iy = np.array([s.getYAstrom() for s in imgsources])
    iflux = np.array([s.getPsfFlux() for s in imgsources])
    I = np.argsort(-iflux)
    # First 200: red dots
    II = I[:200]
    p1 = plt.plot(ix[II], iy[II], 'r.', zorder=10)
    # Rest: tiny dots
    II = I[200:]
    p2 = plt.plot(ix[II], iy[II], 'r.', markersize=1, zorder=9)

    # Ref sources:
    # Only getRa() (not getRaAstrom(), getRaObject()) is non-zero.

    rx, ry = [], []
    for r in refsources:
        xy = wcs.skyToPixel(r.getRaDec())
        rx.append(xy[0])
        ry.append(xy[1])
    rx = np.array(rx)
    ry = np.array(ry)
    p3 = plt.plot(rx, ry, 'bo', mec='b', mfc='none', markersize=6, zorder=20)

    x, y = [], []
    dx, dy = [], []
    for m in matches:
        x0, x1 = m.first.getXAstrom(), m.second.getXAstrom()
        y0, y1 = m.first.getYAstrom(), m.second.getYAstrom()
        #plt.plot([x0, x1], [y0, y1], 'g.-')
        x.append(x0)
        y.append(y0)
        dx.append(x1-x0)
        dy.append(y1-y0)
    #plt.plot(x, y, 's', mec='g', mfc='none', markersize=5)
    p4 = plt.plot(x, y, 'o', mec='g', mfc='g', alpha=0.5, markersize=8, zorder=5)
    p5 = plt.quiver(x, y, dx, dy, angles='xy', scale=30., zorder=30)
    plt.axis('scaled')
    plt.axis([0, W, 0, H])
    # print p1, p2, p3, p4, p5

    plt.figlegend((p1, p2, p3, p4),  # , p5),
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
                   zp=None, delta=False, referrs=None, refstargal=None,
                   title=None,
                   saveplot=True, format='png'):
    print('%i ref sources' % len(refsources))
    print('%i image sources' % len(imgsources))
    print('%i matches' % len(matches))

    # In the "matches" list:
    #    m.first  is catalog
    #    m.second is image

    # In this function, the "m" prefix stands for "matched",
    # "u" stands for "unmatched".

    # *sigh*, turn these into Python lists, so we have the "index" function.
    refsources = [s for s in refsources]
    imgsources = [s for s in imgsources]

    # Now we build numpy int arrays for indexing into the "refsources" and
    # "imgsources" arrays.
    MR = []
    MI = []
    for m in matches:
        try:
            i = refsources.index(m.first)
        except ValueError:
            print('Match list reference source ID', m.first.getSourceId(), 'was not in the list of reference stars')
            continue
        try:
            j = imgsources.index(m.second)
        except ValueError:
            print('Match list source ID', m.second.getSourceId(), 'was not in the list of image sources')
            continue
        MR.append(i)
        MI.append(j)
    MR = np.array(MR)
    MI = np.array(MI)

    # Build numpy boolean arrays for indexing the unmatched stars.
    UR = np.ones(len(refsources), bool)
    UR[MR] = False
    UI = np.ones(len(imgsources), bool)
    UI[MI] = False

    def flux2mag(f):
        return -2.5*np.log10(f)

    refmag = np.array([flux2mag(s.getPsfFlux()) for s in refsources])
    imgflux = np.array([s.getPsfFlux() for s in imgsources])
    imgfluxerr = np.array([s.getPsfFluxErr() for s in imgsources])

    # Cut to fluxes that aren't silly and get mags of matched sources.
    okflux = (imgflux[MI] > 1)
    MI = MI[okflux]
    MR = MR[okflux]

    mimgflux = imgflux[MI]
    mimgmag = flux2mag(mimgflux)
    mimgmagerr = abs(2.5 / np.log(10.) * imgfluxerr[MI] / mimgflux)
    mrefmag = refmag[MR]

    # Get mags of unmatched sources.
    uimgflux = imgflux[UI]
    okflux = (uimgflux > 1)
    uimgmag = flux2mag(uimgflux[okflux])
    urefmag = refmag[UR]

    if False:
        unmatched = [imgsources[i] for i in np.flatnonzero(uimg)]
        uflux = np.array([s.getPsfFlux() for s in unmatched])
        I = np.argsort(-uflux)
        print('Unmatched image sources, by psf flux:')
        print('# FLUX, X, Y, RA, DEC')
        for i in I:
            u = unmatched[i]
            print(u.getPsfFlux(), u.getXAstrom(), u.getYAstrom(), u.getRa(), u.getDec())

        print('Matched image sources, by psf flux:')
        print('# FLUX, X, Y, RA, DEC')
        for i in mimgi:
            m = imgsources[i]
            print(m.getPsfFlux(), m.getXAstrom(), m.getYAstrom(), m.getRa(), m.getDec())

    # Legend entries:
    pp = []
    pl = []

    plt.clf()
    imag = np.append(mimgmag, uimgmag)

    mrefmagerr = None
    if referrs is not None:
        referrs = np.array(referrs)
        mrefmagerr = referrs[MR]

    if refstargal:
        assert(len(refstargal) == len(refsources))
        refstargal = np.array(refstargal).astype(bool)
        ptsets = [(np.logical_not(refstargal[MR]), 'g', 'Matched galaxies', 10),
                  (refstargal[MR], 'b', 'Matched stars', 12)]

    else:
        ptsets = [(np.ones_like(mrefmag).astype(bool), 'b', 'Matched sources', 10)]

    for I, c, leg, zo in ptsets:
        if delta:
            dm = mimgmag[I] - mrefmag[I] + zp
            xi = mrefmag[I]
            yi = dm
            dx = mrefmagerr
            dy = mimgmagerr
        else:
            xi = mimgmag[I]
            yi = mrefmag[I]
            dx = mimgmagerr
            dy = mrefmagerr

        p1 = plt.plot(xi, yi, '.', color=c, mfc=c, mec=c, alpha=0.5, zorder=zo)
        if dx is None or dy is None:
            # errorbars
            xerr, yerr = None, None
            if dx is not None:
                xerr = dx[I]
            if dy is not None:
                yerr = dy[I]
            plt.errorbar(xi, yi, xerr=xerr, yerr=yerr, ecolor=c, fmt=None, zorder=zo)
        else:
            # get the current axis
            ca = plt.gca()
            # add error ellipses
            for j, i in enumerate(np.flatnonzero(I)):
                a = Ellipse(xy=np.array([xi[j], yi[j]]),
                            width=dx[i]/2., height=dy[i]/2.,
                            alpha=0.5, fill=True, ec=c, fc=c, zorder=zo)
                ca.add_artist(a)
        pp.append(p1)
        pl.append(leg)

    if delta:
        m = max(abs(dm))
        plt.axis([np.floor(min(refmag))-0.5, np.ceil(max(refmag)),
                  -m, m])
    else:
        plt.axis([np.floor(min(imag))-0.5, np.ceil(max(imag)),
                  np.floor(min(refmag))-0.5, np.ceil(max(refmag))])
    ax = plt.axis()

    if not delta:
        # Red tick marks show unmatched img sources
        dy = (ax[3]-ax[2]) * 0.05
        y1 = np.ones_like(uimgmag) * ax[3]
        p2 = plt.plot(np.vstack((uimgmag, uimgmag)), np.vstack((y1, y1-dy)), 'r-', alpha=0.5)
        p2 = p2[0]
        # Blue tick marks show matched img sources
        y1 = np.ones_like(mimgmag) * ax[3]
        p3 = plt.plot(np.vstack((mimgmag, mimgmag)), np.vstack((y1-(0.25*dy), y1-(1.25*dy))), 'b-', alpha=0.5)
        p3 = p3[0]
        # Red ticks for unmatched ref sources
        dx = (ax[1]-ax[0]) * 0.05
        x1 = np.ones_like(urefmag) * ax[1]
        p4 = plt.plot(np.vstack((x1, x1-dx)), np.vstack((urefmag, urefmag)), 'r-', alpha=0.5)
        p4 = p4[0]
        # Blue ticks for matched ref sources
        x1 = np.ones_like(mrefmag) * ax[1]
        p5 = plt.plot(np.vstack((x1-(0.25*dx), x1-(1.25*dx))), np.vstack((mrefmag, mrefmag)), 'b-', alpha=0.5)
        p5 = p5[0]

    if zp is not None:
        if delta:
            pzp = plt.axhline(0, linestyle='--', color='b')
        else:
            X = np.array([ax[0], ax[1]])
            pzp = plt.plot(X, X+zp, 'b--')
        pp.append(pzp)
        pl.append('Zeropoint')

    # reverse axis directions.
    if delta:
        plt.axis([ax[1], ax[0], ax[2], ax[3]])
    else:
        plt.axis([ax[1], ax[0], ax[3], ax[2]])

    if band is not None:
        reflabel = 'Reference catalog: %s band (mag)' % band
    else:
        reflabel = 'Reference catalog mag'

    if delta:
        plt.xlabel(reflabel)
        plt.ylabel('Instrumental - Reference (mag)')
        fn = prefix + '-dphotom.' + format

        if zp is not None:
            # Make the plot area smaller to fit the twin axis
            pos = plt.gca().get_position()
            ll = pos.min
            sz = pos.size
            plt.gca().set_position(pos=[ll[0], ll[1], sz[0], sz[1]-0.05])
            # Put the title up top (otherwise it follows the axis)
            if title is not None:
                plt.figtext(0.5, 0.96, title, ha='center', va='top', fontsize='large')
                title = None

            ax2 = plt.twiny()

            # Red tick marks show unmatched img sources
            if zp is not None:
                dy = (ax[3]-ax[2]) * 0.05
                y1 = np.ones_like(uimgmag) * ax[3]
                p2 = plt.plot(np.vstack((uimgmag, uimgmag)) + zp, np.vstack((y1, y1-dy)), 'r-', alpha=0.5)
                p2 = p2[0]
                # Blue tick marks show matched img sources
                y1 = np.ones_like(mimgmag) * ax[3]
                p3 = plt.plot(np.vstack((mimgmag, mimgmag)) + zp,
                              np.vstack((y1-(0.25*dy), y1-(1.25*dy))), 'b-', alpha=0.5)
                p3 = p3[0]
            # Red ticks for unmatched ref sources
            y1 = np.ones_like(urefmag) * ax[2]
            p4 = plt.plot(np.vstack((urefmag, urefmag)), np.vstack((y1, y1+dy)), 'r-', alpha=0.5)
            p4 = p4[0]
            # Blue ticks for matched ref sources
            y1 = np.ones_like(mrefmag) * ax[2]
            p5 = plt.plot(np.vstack((mrefmag, mrefmag)), np.vstack(
                (y1+(0.25*dy), y1+(1.25*dy))), 'b-', alpha=0.5)
            p5 = p5[0]

            plt.xlim(ax[1]-zp, ax[0]-zp)
            plt.xlabel('Instrumental mag')

        legloc = 'lower right'

    else:
        plt.ylabel(reflabel)
        plt.xlabel('Image instrumental mag')
        fn = prefix + '-photom.' + format
        legloc = 'center right'

    if title is not None:
        plt.title(title)

    pp += [p3, p2]
    pl += ['Matched sources', 'Unmatched sources']
    plt.figlegend(pp, pl, legloc, numpoints=1, prop=FontProperties(size='small'))

    P1 = _output(fn, format, saveplot)

    if delta:
        plt.ylim(-0.5, 0.5)
        fn = prefix + '-dphotom2.' + format
        P2 = _output(fn, format, saveplot)
        if not saveplot:
            P1.update(P2)

    return P1


def plotCorrespondences2(imgsources, refsources, matches, wcs, W, H, prefix,
                         saveplot=True, format='png'):
    print('ix,iy')
    ix = np.array([s.getXAstrom() for s in imgsources])
    iy = np.array([s.getYAstrom() for s in imgsources])

    print('rx,ry')
    rx, ry = [], []
    for r in refsources:
        xy = wcs.skyToPixel(r.getRaDec())
        rx.append(xy[0])
        ry.append(xy[1])
    rx = np.array(rx)
    ry = np.array(ry)

    ixy = np.vstack((ix, iy)).T
    rxy = np.vstack((rx, ry)).T

    print('plotshift...')
    cell = 10
    plotshift(ixy, rxy, dcell=cell, ncells=9, W=W, H=H)
    fn = prefix + '-shift1.' + format
    P1 = _output(fn, format, saveplot)

    print('plotshift 2...')
    plt.clf()
    plt.hot()
    plotshift(ixy, rxy, dcell=cell, ncells=9, W=W, H=H, hist=True, nhistbins=2*cell+1)
    fn = prefix + '-shift2.' + format
    P2 = _output(fn, format, saveplot)

    print('plotshift 3...')
    cell = 2
    plotshift(ixy, rxy, dcell=cell, ncells=9, W=W, H=H)
    fn = prefix + '-shift3.' + format
    P3 = _output(fn, format, saveplot)

    print('plotshift 4...')
    plt.clf()
    plt.hot()
    plotshift(ixy, rxy, dcell=cell, ncells=9, W=W, H=H, hist=True, nhistbins=10*cell+1)
    fn = prefix + '-shift4.' + format
    P4 = _output(fn, format, saveplot)

    if not saveplot:
        P1.update(P2)
        P1.update(P3)
        P1.update(P4)
        return P1


def plotCorrespondences(imgsources, refsources, matches, wcs, W, H, prefix):
    ix = np.array([s.getXAstrom() for s in imgsources])
    iy = np.array([s.getYAstrom() for s in imgsources])

    rx, ry = [], []
    for r in refsources:
        xy = wcs.skyToPixel(r.getRaDec())
        rx.append(xy[0])
        ry.append(xy[1])
    rx = np.array(rx)
    ry = np.array(ry)

    # correspondences we could have hit...
    ixy = np.vstack((ix, iy)).T
    rxy = np.vstack((rx, ry)).T
    dcell = 50.
    radius = dcell * np.sqrt(2.)
    # print 'ixy', ixy.shape
    # print 'rxy', rxy.shape

    if False:
        (inds, dists) = spherematch.match(rxy, ixy, radius)
        mi = inds[:, 0]
        ii = inds[:, 1]
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
        print('Found %i matches within %g pixels' % (len(dists), radius))

    ncells = 18.
    cellsize = np.sqrt(W * H / ncells)
    nw = int(round(W / cellsize))
    nh = int(round(H / cellsize))
    # print 'Grid cell size', cellsize
    # print 'N cells', nw, 'x', nh
    edgesx = np.linspace(0, W, nw+1)
    edgesy = np.linspace(0, H, nh+1)

    binx = np.digitize(rx, edgesx)
    biny = np.digitize(ry, edgesy)
    binx = np.clip(binx - 1, 0, nw-1)
    biny = np.clip(biny - 1, 0, nh-1)

    bin = biny * nw + binx

    plt.clf()

    for i in range(nh):
        for j in range(nw):
            thisbin = i * nw + j
            R = (bin == thisbin)
            # print 'cell %i, %i' % (j, i)
            # print '%i ref sources' % sum(R)
            if sum(R) == 0:
                continue
            (inds, dists) = spherematch.match(rxy[R, :], ixy, radius)
            # print 'Found %i matches within %g pixels' % (len(dists), radius)
            ri = inds[:, 0]
            # un-cut ref inds...
            ri = (np.flatnonzero(R))[ri]
            ii = inds[:, 1]

            matchx = rx[ri]
            matchy = ry[ri]
            matchdx = ix[ii] - matchx
            matchdy = iy[ii] - matchy
            ok = (matchdx >= -dcell) * (matchdx <= dcell) * (matchdy >= -dcell) * (matchdy <= dcell)
            #matchx = matchx[ok]
            #matchy = matchy[ok]
            matchdx = matchdx[ok]
            matchdy = matchdy[ok]
            # print 'Cut to %i within %g x %g square' % (sum(ok), dcell*2, dcell*2)

            # Subplot places plots left-to-right, TOP-to-BOTTOM.
            plt.subplot(nh, nw, 1 + ((nh - i - 1)*nw + j))

            plt.plot(matchdx, matchdy, 'ro', mec='r', mfc='r', ms=5, alpha=0.2)
            plt.plot(matchdx, matchdy, 'ro', mec='r', mfc='none', ms=5, alpha=0.2)
            plt.axhline(0, color='k', alpha=0.5)
            plt.axvline(0, color='k', alpha=0.5)
            xticks([], [])
            yticks([], [])
            plt.axis('scaled')
            plt.axis([-dcell, dcell, -dcell, dcell])

    fn = prefix + '-missed.png'
    print('Saving', fn)
    plt.savefig(fn)


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
    print('WCS plots')
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
    cellsize = np.sqrt(W * H / ncells)
    nw = int(np.floor(W / cellsize))
    nh = int(np.floor(H / cellsize))
    # print 'Grid cell size', cellsize
    # print 'N cells', nw, 'x', nh
    cx = np.arange(nw+1) * cellsize + ((W - (nw*cellsize))/2.)
    cy = np.arange(nh+1) * cellsize + ((H - (nh*cellsize))/2.)

    # pixel step size for grid lines
    step = 50

    xx = np.arange(-step, W+2*step, step)
    yy = np.arange(-step, H+2*step, step)

    plt.clf()
    for y in cy:
        dx, dy = [], []
        for x in xx:
            pix = afwGeom.Point2D(x, y)
            distpix = sip.distortPixel(pix)
            dx.append(distpix[0])
            dy.append(distpix[1])
        plt.plot(xx, y*np.ones_like(xx), 'k-', zorder=10)
        dx = np.array(dx)
        dy = np.array(dy)
        if exaggerate != 1:
            dx += (exaggerate * (dx - xx))
            dy += (exaggerate * (dy - y))
        plt.plot(dx, dy, 'r-', zorder=20)

    for x in cx:
        dx, dy = [], []
        for y in yy:
            pix = afwGeom.Point2D(x, y)
            distpix = sip.distortPixel(pix)
            dx.append(distpix[0])
            dy.append(distpix[1])
        plt.plot(x*np.ones_like(yy), yy, 'k-', zorder=10)
        dx = np.array(dx)
        dy = np.array(dy)
        if exaggerate != 1:
            dx += (exaggerate * (dx - x))
            dy += (exaggerate * (dy - yy))
        plt.plot(dx, dy, 'r-', zorder=20)

    plt.axis('scaled')
    plt.axis([0, W, 0, H])

    plt.title(titletxt)

    fn = prefix + suffix + format
    return _output(fn, format, saveplot)
