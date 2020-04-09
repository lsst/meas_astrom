from optparse import OptionParser

import lsst.meas.astrom as measAstrom
from lsst.log import Log
import lsst.meas.algorithms.utils as maUtils

import wcsPlots
import imsimUtils

import numpy as np


def main():
    parser = OptionParser()
    imsimUtils.addOptions(parser)
    parser.add_option('--fixup', dest='fixup', action='store_true', default=False,
                      help='Fix up problems with PT1.1-current outputs')
    parser.add_option('--photometry', '-P', dest='dophotometry', action='store_true',
                      default=False, help='Make photometry plot?')
    parser.add_option('--corr', '-C', dest='docorr', action='store_true',
                      default=False, help='Make correspondences plot?')
    parser.add_option('--distortion', '-D', dest='dodistortion', action='store_true',
                      default=False, help='Make distortion plot?')
    parser.add_option('--matches', '-M', dest='domatches', action='store_true',
                      default=False, help='Make matches plot?')
    parser.add_option('--prefix', '-p', dest='prefix', default='', help='Plot output filename prefix')
    (opt, args) = parser.parse_args()

    if (opt.dophotometry is False
        and opt.docorr is False
        and opt.dodistortion is False
            and opt.domatches is False):
        plots = None
    else:
        plots = []
        if opt.dophotometry:
            plots.append('photom')
        if opt.docorr:
            plots.append('corr')
        if opt.dodistortion:
            plots.append('distortion')
        if opt.domatches:
            plots.append('matches')

    inButler = imsimUtils.getInputButler(opt)

    allkeys = imsimUtils.getAllKeys(opt, inButler)

    for keys in allkeys:
        plotsForField(inButler, keys, opt.fixup, plots, opt.prefix)


def plotsForField(inButler, keys, fixup, plots=None, prefix=''):
    if plots is None:
        plots = ['photom', 'matches', 'corr', 'distortion']

    filters = inButler.queryMetadata('raw', 'filter', **keys)
    print('Filters:', filters)
    filterName = filters[0]

    try:
        psources = inButler.get('icSrc', **keys)
        print('Got sources', psources)
    except Exception:
        print('"icSrc" not found.  Trying "src" instead.')
        psources = inButler.get('src', **keys)
        print('Got sources', psources)

    pmatches = inButler.get('icMatch', **keys)
    sources = psources.getSources()

    calexp = inButler.get('calexp', **keys)
    wcs = calexp.getWcs()

    photocal = calexp.getCalib()
    zp = photocal.getMagnitude(1.)
    print('Zeropoint is', zp)

    # ref sources
    W, H = calexp.getWidth(), calexp.getHeight()

    log = Log.getDefaultLogger()
    log.setLevel(Log.DEBUG)

    kwargs = {}
    if fixup:
        # ugh, mask and offset req'd because source ids are assigned at write-time
        # and match list code made a deep copy before that.
        # (see svn+ssh://svn.lsstcorp.org/DMS/meas/astrom/tickets/1491-b r18027)
        kwargs['sourceIdMask'] = 0xffff
        kwargs['sourceIdOffset'] = -1

    (matches, ref) = measAstrom.generateMatchesFromMatchList(
        pmatches, sources, wcs, W, H, returnRefs=True, log=log, **kwargs)
    print('Got', len(ref), 'reference catalog sources')

    # pull 'stargal' and 'referrs' arrays out of the reference sources
    fdict = maUtils.getDetectionFlags()
    starflag = int(fdict["STAR"])
    stargal = [bool((r.getFlagForDetection() & starflag) > 0)
               for r in ref]
    referrs = [float(r.getPsfInstFluxErr() / r.getPsfInstFlux() * 2.5 / -np.log(10))
               for r in ref]
    nstars = sum([1 for s in stargal if s])
    print('Number of sources with STAR set:', nstars)

    visit = keys['visit']
    raft = keys['raft']
    sensor = keys['sensor']
    prefix += 'imsim-v%i-r%s-s%s' % (visit, raft.replace(',', ''), sensor.replace(',', ''))

    if 'photom' in plots:
        print('photometry plots...')
        tt = 'LSST ImSim v%i r%s s%s' % (visit, raft.replace(',', ''), sensor.replace(',', ''))

        wcsPlots.plotPhotometry(sources, ref, matches, prefix, band=filterName,
                                zp=zp, referrs=referrs, refstargal=stargal, title=tt)
        wcsPlots.plotPhotometry(sources, ref, matches, prefix, band=filterName, zp=zp,
                                delta=True, referrs=referrs, refstargal=stargal, title=tt)

        # test w/ and w/o referrs and stargal.
        if False:
            wcsPlots.plotPhotometry(sources, ref, matches, prefix + 'A', band=filterName, zp=zp, title=tt)
            wcsPlots.plotPhotometry(sources, ref, matches, prefix + 'B',
                                    band=filterName, zp=zp, referrs=referrs, title=tt)
            wcsPlots.plotPhotometry(sources, ref, matches, prefix + 'C',
                                    band=filterName, zp=zp, refstargal=stargal, title=tt)

            wcsPlots.plotPhotometry(sources, ref, matches, prefix + 'A',
                                    band=filterName, zp=zp, delta=True, title=tt)
            wcsPlots.plotPhotometry(sources, ref, matches, prefix + 'B', band=filterName,
                                    zp=zp, delta=True, referrs=referrs, title=tt)
            wcsPlots.plotPhotometry(sources, ref, matches, prefix + 'C', band=filterName,
                                    zp=zp, delta=True, refstargal=stargal, title=tt)

    if 'matches' in plots:
        print('matches...')
        wcsPlots.plotMatches(sources, ref, matches, wcs, W, H, prefix)

    if 'distortion' in plots:
        print('distortion...')
        wcsPlots.plotDistortion(wcs, W, H, 400, prefix,
                                'SIP Distortion (exaggerated x 10)', exaggerate=10.)
        print('distortion...')
        wcsPlots.plotDistortion(wcs, W, H, 400, prefix,
                                'SIP Distortion (exaggerated x 100)', exaggerate=100.,
                                suffix='-distort2.')


if __name__ == '__main__':
    main()
