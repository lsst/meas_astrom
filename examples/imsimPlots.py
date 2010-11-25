from optparse import OptionParser
from math import hypot

import lsst.daf.persistence as dafPersist
import lsst.pex.policy as policy
import lsst.meas.astrom as measAstrom
import lsst.afw.image as afwImage
from lsst.pex.logging import Log
from lsst.afw.coord import DEGREES
from lsst.obs.lsstSim import LsstSimMapper

import wcsPlots


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-i', '--input', dest='inRoot', default='.', help='input root')
    parser.add_option('-R', '--registry', help='registry', dest='registry')
    parser.add_option('-v', '--visit', type='int', dest='visit')
    parser.add_option('-r', '--raft', dest='raft')
    parser.add_option('-s', '--sensor', dest='sensor')
    (opt, args) = parser.parse_args()

    mapper = LsstSimMapper(root=opt.inRoot, registry=opt.registry)
    bf = dafPersist.ButlerFactory(mapper=mapper)
    inButler = bf.create()

    #bf = dafPersist.ButlerFactory(mapper=LsstSimMapper(
    #    root=opt.outRoot, registry=opt.registry))
    #outButler = bf.create() 

    keys = { 'visit': opt.visit,
             'raft': opt.raft,
             'sensor': opt.sensor,
             }

    filters = inButler.queryMetadata('raw', 'filter', **keys)
    print 'Filters:', filters
    filterName = filters[0]

    psources = inButler.get('icSrc', **keys)
    pmatches = inButler.get('icMatch', **keys)
    print 'Got sources', psources
    print 'Got matches', pmatches
    #print '  ', type(matches)
    #print matches.__subject__
    matchmeta = pmatches.getSourceMatchMetadata()
    matches = pmatches.getSourceMatches()

    sources = psources.getSources()
    
    calexp = inButler.get('calexp', **keys)
    print 'Got calexp', calexp
    wcs = calexp.getWcs()
    print 'Got wcs', wcs
    wcs = afwImage.cast_TanWcs(wcs)
    print 'After cast:', wcs
    
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
    ref = solver.getCatalogue(ra, dec, radius, filterName, idName, anid)
    print 'Got', len(ref), 'reference catalog sources'

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

    measAstrom.joinMatchList(matches, ref, first=True, log=log)
    measAstrom.joinMatchList(matches, sources, first=False, log=log)

    #for m in matches:
    #    x0,x1 = m.first.getXAstrom(), m.second.getXAstrom()
    #    y0,y1 = m.first.getYAstrom(), m.second.getYAstrom()
    #    print 'x,y, dx,dy', x0, y0, x1-x0, y1-y0

    prefix = 'imsim-v%i-r%s-s%s' % (opt.visit, opt.raft.replace(',',''), opt.sensor.replace(',',''))

    wcsPlots.plotMatches(sources, ref, matches, wcs, W, H, prefix)
    wcsPlots.plotPhotometry(sources, ref, matches, prefix, band=filterName)
    
    wcsPlots.plotCorrespondences2(sources, ref, matches, wcs, W, H, prefix)
    wcsPlots.plotCorrespondences(sources, ref, matches, wcs, W, H, prefix)

    wcsPlots.plotDistortion(wcs, W, H, 400, prefix,
                            'SIP Distortion (exaggerated x 10)', exaggerate=10.)
