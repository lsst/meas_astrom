import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

import os

import lsst.meas.astrom as measAstrom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.daf.base  as dafBase
import lsst.pex.logging as pexLog

import pyfits

def main():
    mydir = os.path.dirname(__file__)

    # Read sources
    # fitscopy coaddSources.fits"[col x=centroid_sdss[1]; y=centroid_sdss[2]; flux=flux_psf]" xy.fits
    fn = os.path.join(mydir, 'xy2710.fits')
    P = pyfits.open(fn)[1].data
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
        #src = srcTable.makeRecord()
        #src = srcs.makeRecord()
        src = srcs.addNew()
        if not (np.isfinite(xi) and np.isfinite(yi)):
            continue
        #print 'x,y', xi,yi
        src.set(key.getX(), xi)
        src.set(key.getY(), yi)
        src.set(fkey, fi)
    # for src in srcs:
    # print src.getX(), src.getY()
        
    x0,y0 = 335750, 223750
    W,H = 8500, 12500
    imargs = dict(imageSize=(W,H), filterName='i', x0=x0, y0=y0)
    
    # Read original WCS
    fn = os.path.join(mydir, 't2710.wcs')
    #    print help(afwImage.Wcs.readFits)
    hdr = afwImage.readMetadata(fn)
    wcs0 = afwImage.makeWcs(hdr)
    #wcs0 = afwImage.Wcs.readFits(fn, 0)
    #print 'WCS0:', wcs0
    #print 'wcs0:', wcs0.getFitsMetadata().toString()

    # They call measAstrom.determineWcs() to get a new WCS with SIP distortion
    conf = measAstrom.MeasAstromConfig(sipOrder=4)

    andConfig = measAstrom.AstrometryNetDataConfig()
    anddir = os.path.join(mydir, 'astrometry_net_data', 'ticket2710')
    os.environ['ASTROMETRY_NET_DATA_DIR'] = anddir
    fn = os.path.join(anddir, 'andConfig.py')
    andConfig.load(fn)

    ast = measAstrom.Astrometry(conf, andConfig, logLevel=pexLog.Log.DEBUG)

    refs = ast.getReferenceSourcesForWcs(wcs0, **imargs)
    print 'Got', len(refs), 'reference objects for initial WCS'


    conf2 = measAstrom.MeasAstromConfig(sipOrder=4, calculateSip=False)
    ast2 = measAstrom.Astrometry(conf2, andConfig, logLevel=pexLog.Log.DEBUG)
    solve = ast2.determineWcs2(srcs, **imargs)
    tanwcs = solve.getTanWcs()

    wcs2 = ast.getSipWcsFromWcs(wcs0, tanwcs, (W,H), x0=x0, y0=y0)


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


    from astrometry.libkd.spherematch import match
    from astrometry.util.plotutils import plothist,PlotSequence

    ps = PlotSequence('t2710')

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

    
    matches = solve.getTanMatches()
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



    solve = ast.determineWcs2(srcs, **imargs)
    wcs1 = solve.getSipWcs()
    #print 'wcs1:', wcs1.getFitsMetadata().toString()


    matches = solve.getSipMatches()
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
    
    
    

if __name__ == '__main__':
    main()
    
