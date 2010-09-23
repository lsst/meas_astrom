import os
import sys
import types
from optparse import OptionParser

import lsst.meas.astrom as measAstrom
import lsst.pex.policy  as pexPolicy
import lsst.pex.logging as pexLog
import lsst.daf.persistence              as dafPersist
import lsst.daf.base                     as dafBase
import lsst.afw.image                    as afwImage
import lsst.afw.detection                as afwDet

import lsst.meas.astrom.net as astromNet
import lsst.meas.photocal as photocal

from astrometry.util.pyfits_utils import *
from numpy import array,unique

from sourceset_boost_to_fits import *

class ducky(object):
    def __str__(self):
        return 'ducky: ' + '; '.join([''+str(k)+'='+str(v) for k,v in self.__dict__.items()])

class fakeExposure(ducky):
    def __init__(self):
        self.wcs = None
        self.width, self.height = None,None
        self.xy0 = None
        self.filterName = None

    def getMaskedImage(self):
        fakey = ducky()
        fakey.xy0 = self.xy0
        fakey.getXY0 = lambda x: x.xy0
        return fakey

    def getFilter(self):
        fakey = ducky()
        fakey.filterName = self.filterName
        # advanced duck-typing: here we give the duck a member function "getName()"
        fakey.getName = types.MethodType(lambda x: x.filterName, fakey, fakey.__class__)
        return fakey

    def getWcs(self):
        return self.wcs
    def setWcs(self, wcs):
        self.wcs = wcs
    def hasWcs(self):
        return self.wcs is not None

    def getWidth(self):
        return self.width
    def setWidth(self, width):
        self.width = width

    def getHeight(self):
        return self.height
    def setHeight(self, height):
        self.height = height
        


def rerun(sourceset, policy=None, exposure=None, wcs=None,
          W=None, H=None, xy0=None, filtername=None, log=None, fieldname=None,
          plotprefix='', outwcsfn='out.wcs', doPlots=False):

    if exposure is None:
        # Create a duck of the appropriate quackiness.
        exposure = fakeExposure()
        exposure.xy0 = xy0
        exposure.filterName = filtername

    if wcs is not None:
        exposure.setWcs(wcs)
    if W is not None:
        exposure.setWidth(W)
    if H is not None:
        exposure.setHeight(H)

    if policy is None:
        policy = pexPolicy.Policy.createPolicy(pexPolicy.PolicyString(
            '''#<?cfg paf policy?>
            matchThreshold: 30
            numBrightStars: 1000
            blindSolve: true
            distanceForCatalogueMatchinArcsec: 5
            cleaningParameter: 3
            calculateSip: True
            sipOrder: 2'''))

    if exposure.getWidth() is None or exposure.getHeight() is None:
        print 'Warning: exposure image width or height is None'

    doTrim = False

    print 'Exposure:', exposure
    print 'Filter:', exposure.getFilter()
    print 'Filter name:', exposure.getFilter().getName()
    print
    print 'determineWcs()...'
    print

    # *sigh*.
    path=os.path.join(os.environ['ASTROMETRY_NET_DATA_DIR'], "metadata.paf")
    solver = astromNet.GlobalAstrometrySolution(path, log)
    matchThreshold = policy.get('matchThreshold')
    solver.setMatchThreshold(matchThreshold)

    (matchList,wcs,refstars,plots) = measAstrom.determineWcs(policy, exposure, sourceset, log=log,
                                                              doTrim=doTrim, solver=solver,
                                                              returnRefStars=True,
                                                              returnPlotData=True, plotFormat='png')

    print
    print 'determineWcs() finished.  Got:'
    print
    if wcs is None:
        print 'WCS determination failed.'
    else:
        print '%i matches' % len(matchList)

        for sm in matchList:
            print '  id %i (%.1f, %.1f) flux %.3g   ---   id %i (%.1f, %.1f) flux %.3g' % (
                sm.first.getId(),  sm.first.getXAstrom(), sm.first.getYAstrom(), sm.first.getPsfFlux(),
                sm.second.getId(), sm.second.getXAstrom(), sm.second.getYAstrom(), sm.second.getPsfFlux())

        srcids = array([sm.second.getId() for sm in matchList])
        print '%i matches; %i unique image sources' % (len(matchList), len(unique(srcids)))
        assert(len(srcids) == len(unique(srcids)))

        fitshdr = wcs.getFitsMetadata()
        print 'Found WCS:'
        print fitshdr.toString()

        # No dice: daf_persistence can't write PropertySets to FITS.
        if False:
            loc = dafPersist.LogicalLocation(outfn)
            storageList = dafPersist.StorageList()
            additionalData = dafBase.PropertySet()
            persistence = dafPersist.Persistence.getPersistence(pexPolicy.Policy())
            storageList.append(persistence.getPersistStorage("FitsStorage", loc))
            persistence.persist(fitshdr, storageList, additionalData)
            
        # Add IMAGEW, IMAGEH header cards.
        fitshdr.add('IMAGEW', W)
        fitshdr.add('IMAGEH', H)

        # Create a fake Image in order to abuse its writeFits() method.
        im = afwImage.ImageF()
        im.writeFits(outwcsfn, fitshdr)

    # Write matchList.

    # Do photocal too.
    print 'Doing photocal...'
    magObj = photocal.calcPhotoCal(matchList, log=log, goodFlagValue=0)
    print 'got:', magObj

    if doPlots:
        import lsst.meas.astrom.wcsPlots as wcsPlots
        wcsPlots.plotPhotometry(sourceset, refstars, matchList, prefix='photocal',
                                saveplot=True)
        from pylab import axis,plot,savefig,title
        ax = axis()
        zp = magObj.getMag(1.)
        print 'Zero-point:', zp
        plot([ax[0], ax[1]], [ax[0]+zp, ax[1]+zp], 'b-')
        axis(ax)
        if fieldname is not None:
            title(fieldname)
        savefig(plotprefix + 'photocal-zp.png')
        
        print 'plots:'
        for k,v in plots.items():
            print '  ',k,'len', len(v)
            f = open(plotprefix + k, 'wb')
            f.write(v)
            f.close()


def rerun_main(sysargs):
    parser = OptionParser(usage='%prog [options] <*.boost or *.fits SourceSets>')
    parser.add_option('-W', '--width', dest='width', type='int', help='Image width (pixels)')
    parser.add_option('-H', '--height', dest='height', type='int', help='Image height (pixels)')
    parser.add_option('-f', '--filter', dest='filter', help='Filter name')
    parser.add_option('-p', '--prefix', dest='plotprefix', help='Plot filename prefix')
    parser.add_option('-P', '--plots', dest='doplots', help='Make plots?', action='store_true')
    parser.add_option('-v', '--verbose', dest='verb', help='+verbose', action='store_true')
    parser.set_defaults(width=None, height=None, filter=None, verb=False, plotprefix='', doplots=False)
    opt,args = parser.parse_args(sysargs)
    if len(args) == 0:
        parser.print_help()
        return -1

    level = pexLog.Log.DEBUG if opt.verb else pexLog.Log.INFO
    log = pexLog.Log(pexLog.Log.getDefaultLog(), "rerun-wcs", level);

    for fn in args:
        if fn.endswith('.boost'):
            print 'Reading boost-format', fn
            ss = sourceset_read_boost(fn)
        else:
            print 'Reading FITS-format', fn
            T = fits_table(fn)
            ss = afwDet.SourceSet()
            C = T.columns()
            # FITS -> Source.setXXX()
            columnmap = { #'SourceId':'id',
                'XAstrom':'x',
                'YAstrom':'y',
                'XAstromErr': 'xerr',
                'YAstromErr': 'yerr',
                'Ra': 'ra',
                'Dec': 'dec',
                'Ixx': 'ixx',
                'Iyy': 'iyy',
                'PsfFlux': 'f_psf',
                'ApFlux': 'f_ap',
                'Id':'id',
                }
            conversions = {'id':int}

            for k,v in columnmap.items():
                if not v in C:
                    print 'Warning, column', v, 'is not in the FITS table -- won\'t set Source\'s', k
            for i in range(len(T)):
                src = afwDet.Source()
                ss.append(src)
            for k,v in columnmap.items():
                if not v in C:
                    continue
                totype = conversions.get(v, float)
                print 'Setting source attribute',k,'from FITS column', v
                for src,val in zip(ss, T.getcolumn(v)):
                    #print 'val is', type(val), val
                    getattr(src, 'set'+k)(totype(val))

        print 'Read %i sources' % (len(ss))
        srcids = array([s.getId() for s in ss])
        print '%i sources; %i unique IDs' % (len(ss), len(unique(srcids)))

        #if opt.verb:
        #    for i in range(100):
        #        print '  (%.1f, %.1f) psf flux %.1f' % (ss[i].getXAstrom(), ss[i].getYAstrom(), ss[i].getPsfFlux())
                                
        rerun(ss, W=opt.width, H=opt.height, filtername=opt.filter, log=log, fieldname=fn,
              plotprefix=opt.plotprefix, doPlots=opt.doplots)
        return 0
    
if __name__ == '__main__':
    sys.exit(rerun_main(sys.argv[1:]))

