import sys
import types
from optparse import OptionParser

import lsst.meas.astrom as measAstrom
import lsst.pex.policy  as pexPolicy
import lsst.pex.logging as pexLog
import lsst.daf.persistence              as dafPersist
import lsst.daf.base                     as dafBase
import lsst.afw.image                    as afwImage

from astrometry.util.pyfits_utils import *
from numpy import array

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
          W=None, H=None, xy0=None, filtername=None, log=None):

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
            numBrightStars: 50
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

    (matchList,wcs) = measAstrom.determineWcs(policy, exposure, sourceset, log=log,
                                              doTrim=doTrim)

    print
    print 'determineWcs() finished.  Got:'
    print
    print '%i matches' % len(matchList)
    #print 'WCS: ', wcs

    fitshdr = wcs.getFitsMetadata()
    #print 'FITS:', fitshdr
    print 'Found WCS:'
    print fitshdr.toString()

    # No dice: daf_persistence can't write PropertySets to FITS.
    if False:
        outfn = 'out.wcs'
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
    im.writeFits('out.wcs', fitshdr)



if __name__ == '__main__':
    parser = OptionParser(usage='%prog [options] <*.boost SourceSets>')
    # or *.fits of sourcesets>')
    parser.add_option('-W', '--width', dest='width', type='int', help='Image width (pixels)')
    parser.add_option('-H', '--height', dest='height', type='int', help='Image height (pixels)')
    parser.add_option('-f', '--filter', dest='filter', help='Filter name')
    parser.add_option('-v', '--verbose', dest='verb', help='+verbose', action='store_true')
    #parser.add_option('-x', '--x-column', dest='xcol', help='X column name (for FITS inputs)')
    parser.set_defaults(width=None, height=None, filter=None, verb=False)
    opt,args = parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        sys.exit(0)

    level = pexLog.Log.DEBUG if opt.verb else pexLog.Log.INFO
    log = pexLog.Log(pexLog.Log.getDefaultLog(), "rerun-wcs", level);

    for fn in args:
        print 'Reading', fn
        ss = sourceset_read_boost(fn)
        print 'Read %i sources' % (len(ss))

        rerun(ss, W=opt.width, H=opt.height, filtername=opt.filter, log=log)

