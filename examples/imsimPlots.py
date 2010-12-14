from optparse import OptionParser
from math import hypot

import lsst.pex.policy as policy
import lsst.meas.astrom as measAstrom
import lsst.afw.image as afwImage
from lsst.pex.logging import Log
from lsst.afw.coord import DEGREES

import wcsPlots
import imsimUtils

def main():
    parser = OptionParser()
    imsimUtils.addOptions(parser)
    parser.add_option('--fixup', dest='fixup', action='store_true', default=False, help='Fix up problems with PT1.1-current outputs')
    parser.add_option('--photometry', '-P', dest='dophotometry', action='store_true', default=False, help='Make photometry plot?')
    parser.add_option('--corr', '-C', dest='docorr', action='store_true', default=False, help='Make correspondences plot?')
    parser.add_option('--distortion', '-D', dest='dodistortion', action='store_true', default=False, help='Make distortion plot?')
    (opt, args) = parser.parse_args()

    if (opt.dophotometry is False and
        opt.docorr       is False and
        opt.dodistortion is False):
        plots = None
    else:
        plots = []
        if opt.dophotometry:
            plots.append('photom')
        if opt.docorr:
            plots.append('corr')
        if opt.doistortion:
            plots.append('distortion')
        
    inButler = imsimUtils.getInputButler(opt)

    allkeys = imsimUtils.getAllKeys(opt, inButler)

    for keys in allkeys:
        plotsForField(inButler, keys, opt.fixup, plots)


def plotsForField(inButler, keys, fixup, plots=None):
    if plots is None:
        plots = ['photom','corr','distortion']
        
    filters = inButler.queryMetadata('raw', 'filter', **keys)
    print 'Filters:', filters
    filterName = filters[0]

    psources = inButler.get('icSrc', **keys)
    pmatches = inButler.get('icMatch', **keys)
    print 'Got sources', psources
    print 'Got matches', pmatches
    matchmeta = pmatches.getSourceMatchMetadata()
    matches = pmatches.getSourceMatches()
    print 'Match metadata:', matchmeta
    sources = psources.getSources()
    
    calexp = inButler.get('calexp', **keys)
    print 'Got calexp', calexp
    wcs = calexp.getWcs()
    print 'Got wcs', wcs
    print wcs.getFitsMetadata().toString()
    wcs = afwImage.cast_TanWcs(wcs)
    print 'After cast:', wcs

    photocal = calexp.getCalib()
    zp = photocal.getMagnitude(1.)
    print 'Zeropoint is', zp

    # ref sources
    W,H = calexp.getWidth(), calexp.getHeight()
    xc,yc = W/2., H/2.
    radec = wcs.pixelToSky(xc, yc)
    ra = radec.getLongitude(DEGREES)
    dec = radec.getLatitude(DEGREES)
    radius = wcs.pixelScale() * hypot(xc, yc) * 1.1
    print 'Image W,H', W,H
    print 'Image center RA,Dec', ra, dec
    print 'Searching radius', radius, 'arcsec'
    log = Log.getDefaultLog()
    log.setThreshold(Log.DEBUG);
    pol = policy.Policy()
    pol.set('matchThreshold', 30)
    solver = measAstrom.createSolver(pol, log)
    idName = 'id'
    # could get this from matchlist meta...
    anid = matchmeta.getInt('ANINDID')
    print 'Searching index with ID', anid
    print 'Using ID column name', idName
    print 'Using filter column name', filterName
    ref = solver.getCatalogue(ra, dec, radius, filterName, idName, anid)
    print 'Got', len(ref), 'reference catalog sources'

    if False:
        extras = solver.getCatalogueExtra(ra, dec, radius, ['id', 'starnotgal', filterName + '_err'], anid)
        #print 'Got extras:', extras
        #pass
        #inds = solver.getIndexList()
        #print 'Got inds:', inds
        #for ind in inds:
        #    print 'ind', ind
        #    print dir(ind)
        #X = starkd_search_stars_in_field(skdt, tanwcs, margin)
        

    keepref = []
    for i in xrange(len(ref)):
        #print ref[i].getXAstrom(), ref[i].getYAstrom(), ref[i].getRa(), ref[i].getDec()
        x,y = wcs.skyToPixel(ref[i].getRa(), ref[i].getDec())
        if x < 0 or y < 0 or x > W or y > H:
            continue
        ref[i].setXAstrom(x)
        ref[i].setYAstrom(y)
        keepref.append(ref[i])
    print 'Kept', len(keepref), 'reference sources'
    ref = keepref

    if False:
        m0 = matches[0]
        f,s = m0.first, m0.second
        print 'match 0: ref %i, source %i' % (f.getSourceId(), s.getSourceId())
        print '  ref x,y,flux = (%.1f, %.1f, %.1f)' % (f.getXAstrom(), f.getYAstrom(), f.getPsfFlux())
        print '  src x,y,flux = (%.1f, %.1f, %.1f)' % (s.getXAstrom(), s.getYAstrom(), s.getPsfFlux())

    measAstrom.joinMatchList(matches, ref, first=True, log=log)
    args = {}
    if fixup:
        # ugh, mask and offset req'd because source ids are assigned at write-time
        # and match list code made a deep copy before that.
        # (see svn+ssh://svn.lsstcorp.org/DMS/meas/astrom/tickets/1491-b r18027)
        args['mask'] = 0xffff
        args['offset'] = -1
    measAstrom.joinMatchList(matches, sources, first=False, log=log, **args)

    if False:
        for m in matches:
            x0,x1 = m.first.getXAstrom(), m.second.getXAstrom()
            y0,y1 = m.first.getYAstrom(), m.second.getYAstrom()
            print 'x,y, dx,dy', x0, y0, x1-x0, y1-y0

    if False:
        m0 = matches[0]
        f,s = m0.first, m0.second
        print 'match 0: ref %i, source %i' % (f.getSourceId(), s.getSourceId())
        print '  ref x,y,flux = (%.1f, %.1f, %.1f)' % (f.getXAstrom(), f.getYAstrom(), f.getPsfFlux())
        print '  src x,y,flux = (%.1f, %.1f, %.1f)' % (s.getXAstrom(), s.getYAstrom(), s.getPsfFlux())
        r,d = 2.31262000000000, 3.16386000000000
        x,y = wcs.skyToPixel(r,d)
        print 'x,y', x,y
        r2d2 = wcs.pixelToSky(x,y)
        r2 = r2d2.getLongitude(DEGREES)
        d2 = r2d2.getLatitude(DEGREES)
        print r,d
        print r2,d2


    visit = keys['visit']
    raft = keys['raft']
    sensor = keys['sensor']
    prefix = 'imsim-v%i-r%s-s%s' % (visit, raft.replace(',',''), sensor.replace(',',''))

    if 'photom' in plots:
        wcsPlots.plotPhotometry(sources, ref, matches, prefix, band=filterName, zp=zp)

    if 'matches' in plots:
        wcsPlots.plotMatches(sources, ref, matches, wcs, W, H, prefix)

    if 'corr' in plots:
        wcsPlots.plotCorrespondences2(sources, ref, matches, wcs, W, H, prefix)
        wcsPlots.plotCorrespondences(sources, ref, matches, wcs, W, H, prefix)

    if 'distortion' in plots:
        wcsPlots.plotDistortion(wcs, W, H, 400, prefix,
                                'SIP Distortion (exaggerated x 10)', exaggerate=10.)

        wcsPlots.plotDistortion(wcs, W, H, 400, prefix,
                                'SIP Distortion (exaggerated x 100)', exaggerate=100.,
                                suffix='-distort2.')


if __name__ == '__main__':
    main()

