import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

import os

import lsst.meas.astrom as measAstrom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.pex.logging as pexLog

import astropy.io.fits as astropy

'''
This file was produced by dstn trying to reproduce and diagnose the
error reported in ticket #2710, ie, that computing SIP polynomials
wasn't working for sources with positions very far from the origin.

The how-to-reproduce data included a "coaddSources.fits" table of
source positions produced by the standard afwTable tools, and a 1 GB
exposure containing a STG WCS.  I extracted the x,y, and flux columns
from the source table into xy2710.fits; I was also testing with stock
astrometry.net tools.  I pulled out the original WCS header into
t2710.wcs to avoid adding a 1 GB file here.  I also pulled out the
image offset x0,y0 size W,H, and filter.

One gotcha: you need an astrometry_net more recent than 0.30, because
I'm using the "spherematch" module there to match up sources.
'''

def showSipSolutions(srcs, wcs0, andDir, x0, y0, W, H, filterName,
                     plotPrefix):
                     
    '''
    srcs: afw Catalog of sources
    wcs0: original WCS
    andDir: astrometry_net_data directory
    '''
    imargs = dict(imageSize=(W,H), filterName=filterName, x0=x0, y0=y0)

    # Set up astrometry_net_data
    os.environ['ASTROMETRY_NET_DATA_DIR'] = andDir
    andConfig = measAstrom.AstrometryNetDataConfig()
    fn = os.path.join(andDir, 'andConfig.py')
    andConfig.load(fn)

    # Set up meas_astrom
    conf = measAstrom.ANetBasicAstrometryConfig(sipOrder=4)
    ast = measAstrom.ANetBasicAstrometryTask(conf, andConfig, logLevel=pexLog.Log.DEBUG)

    # What reference sources are in the original WCS
    refs = ast.getReferenceSourcesForWcs(wcs0, **imargs)
    print 'Got', len(refs), 'reference objects for initial WCS'

    # How does a straight TAN solution look?
    conf2 = measAstrom.ANetBasicAstrometryConfig(sipOrder=4, calculateSip=False)
    ast2 = measAstrom.ANetBasicAstrometryTask(conf2, andConfig, logLevel=pexLog.Log.DEBUG)
    solve = ast2.determineWcs2(srcs, **imargs)
    tanwcs = solve.tanWcs

    # How about if we fit a SIP WCS using the *original* WCS?
    wcs2 = ast.getSipWcsFromWcs(wcs0, (W,H), x0=x0, y0=y0)

    # (We determineWcs() for a SIP solution below...)
    
    # Make some plots in pixel space by pushing ref sources through WCSes
    rx0,ry0 = [],[]
    rx2,ry2 = [],[]
    rx3,ry3 = [],[]
    for src in refs:
        xy = wcs0.skyToPixel(src.getCoord())
        rx0.append(xy[0])
        ry0.append(xy[1])

        xy = tanwcs.skyToPixel(src.getCoord())
        rx2.append(xy[0])
        ry2.append(xy[1])

        xy = wcs2.skyToPixel(src.getCoord())
        rx3.append(xy[0])
        ry3.append(xy[1])
        
    rx0 = np.array(rx0)
    ry0 = np.array(ry0)
    rx2 = np.array(rx2)
    ry2 = np.array(ry2)
    rx3 = np.array(rx3)
    ry3 = np.array(ry3)

    x = np.array([src.getX() for src in srcs])
    y = np.array([src.getY() for src in srcs])
    
    from astrometry.libkd.spherematch import match
    from astrometry.util.plotutils import plothist,PlotSequence

    ps = PlotSequence(plotPrefix)

    # Match up various sources...
    R = 2.
    
    II,d = match(np.vstack((x,y)).T, np.vstack((rx0,ry0)).T, R)
    I = II[:,0]
    J = II[:,1]

    pa = dict(range=((-R,R),(-R,R)))
    
    plt.clf()
    plothist(x[I]-rx0[J], y[I]-ry0[J], 200, **pa)
    plt.title('Source positions - Reference positions (initial WCS)')
    plt.xlabel('delta-X (pixels)')
    plt.ylabel('delta-Y (pixels)')
    ps.savefig()

    II,d = match(np.vstack((x,y)).T, np.vstack((rx2,ry2)).T, R)
    I = II[:,0]
    J = II[:,1]
    
    plt.clf()
    plothist(x[I]-rx2[J], y[I]-ry2[J], 200, **pa)
    plt.title('Source positions - Reference positions (TAN WCS)')
    plt.xlabel('delta-X (pixels)')
    plt.ylabel('delta-Y (pixels)')
    ps.savefig()

    II,d = match(np.vstack((x,y)).T, np.vstack((rx3,ry3)).T, R)
    I = II[:,0]
    J = II[:,1]
    plt.clf()
    plothist(x[I]-rx3[J], y[I]-ry3[J], 200, **pa)
    plt.title('Source positions - Reference positions (SIP WCS #2)')
    plt.xlabel('delta-X (pixels)')
    plt.ylabel('delta-Y (pixels)')
    ps.savefig()

    II,d = match(np.vstack((rx0,ry0)).T, np.vstack((rx3,ry3)).T, R)
    I = II[:,0]
    J = II[:,1]
    plt.clf()
    plothist(rx0[I]-rx3[J], ry0[I]-ry3[J], 200, **pa)
    plt.title('Reference positions (Original WCS) - Reference positions (SIP WCS #2)')
    plt.xlabel('delta-X (pixels)')
    plt.ylabel('delta-Y (pixels)')
    ps.savefig()

    
    matches = solve.tanMatches
    msx,msy = [],[]
    mrx,mry = [],[]
    for m in matches:
        ref,src = m.first, m.second
        xy = tanwcs.skyToPixel(ref.getCoord())
        mrx.append(xy[0])
        mry.append(xy[1])
        msx.append(src.getX())
        msy.append(src.getY())
    
    plt.clf()
    #plt.plot(x, y, 'o', mec='r', mfc='none', ms=4)
    plt.plot(x, y, 'r.')
    plt.plot(msx, msy, 'o', mec='r')
    #plt.plot(rx0, ry0, 'o', mec='g', mfc='none', ms=4)
    plt.plot(rx0, ry0, 'g.')
    plt.plot(mrx, mry, 'gx')
    plt.title('TAN matches')
    ps.savefig()

    # Get SIP solution (4th order)
    
    solve = ast.determineWcs2(srcs, **imargs)
    wcs1 = solve.sipWcs
    #print 'wcs1:', wcs1.getFitsMetadata().toString()

    matches = solve.sipMatches
    msx,msy = [],[]
    mrx,mry = [],[]
    for m in matches:
        ref,src = m.first, m.second
        xy = tanwcs.skyToPixel(ref.getCoord())
        mrx.append(xy[0])
        mry.append(xy[1])
        msx.append(src.getX())
        msy.append(src.getY())
    
    plt.clf()
    plt.plot(x, y, 'r.')
    plt.plot(msx, msy, 'o', mec='r')
    plt.plot(rx0, ry0, 'g.')
    plt.plot(mrx, mry, 'gx')
    plt.title('SIP matches')
    ps.savefig()
    
    rx1,ry1 = [],[]
    for src in refs:
        xy = wcs1.skyToPixel(src.getCoord())
        rx1.append(xy[0])
        ry1.append(xy[1])
    rx1 = np.array(rx1)
    ry1 = np.array(ry1)

    plt.clf()
    plt.plot(x, y, 'o', mec='r', mfc='none')
    plt.plot(rx0, ry0, 'bx')
    plt.plot(rx1, ry1, 'g+')
    plt.plot(rx2, ry2, 'mx')
    plt.plot(rx3, ry3, 'r+')
    ps.savefig()

    plt.axis([x0, x0+500, y0, y0+500])
    ps.savefig()

    II,d = match(np.vstack((x,y)).T, np.vstack((rx1,ry1)).T, R)
    I = II[:,0]
    J = II[:,1]
    
    plt.clf()
    plothist(x[I]-rx1[J], y[I]-ry1[J], 200, **pa)
    plt.title('Source positions - Reference positions (SIP WCS)')
    plt.xlabel('delta-X (pixels)')
    plt.ylabel('delta-Y (pixels)')
    ps.savefig()


def readSourcesFromXyTable(xyfn):
    '''
    Read sources from a plain FITS table of x,y,flux and put
    into an afw Catalog.
    '''
    P = pyfits.open(xyfn)[1].data
    x = P['x']
    y = P['y']
    f = P['flux']
    srcSchema = afwTable.SourceTable.makeMinimalSchema()
    key = srcSchema.addField("centroid", type="PointD")
    srcSchema.addField("centroid.flags", type="Flag")
    srcSchema.addField("centroid.err", type="CovPointF")
    fkey = srcSchema.addField("flux", type=float)
    srcTable = afwTable.SourceTable.make(srcSchema)
    srcTable.defineCentroid("centroid")
    srcTable.definePsfFlux("flux")
    srcs = afwTable.SourceCatalog(srcTable)
    for xi,yi,fi in zip(x,y,f):
        src = srcs.addNew()
        if not (np.isfinite(xi) and np.isfinite(yi)):
            continue
        src.set(key.getX(), xi)
        src.set(key.getY(), yi)
        src.set(fkey, fi)
    return srcs
    
if __name__ == '__main__':
    mydir = os.path.dirname(__file__)
    # fitscopy coaddSources.fits"[col x=centroid_sdss[1]; y=centroid_sdss[2]; flux=flux_psf]" xy.fits
    xyfn = os.path.join(mydir, 'xy2710.fits')
    sources = readSourcesFromXyTable(xyfn)

    x0,y0 = 335750, 223750
    W,H = 8500, 12500
    filterName = 'i'    

    # Read original WCS
    fn = os.path.join(mydir, 't2710.wcs')
    hdr = afwImage.readMetadata(fn)
    wcs0 = afwImage.makeWcs(hdr)

    anddir = os.path.join(mydir, 'astrometry_net_data', 'ticket2710')

    plotPrefix = 't2710'

    showSipSolutions(sources, wcs0, anddir, x0, y0, W, H, filterName,
                     plotPrefix)
