import sys
from optparse import OptionParser

import lsst.meas.astrom as measAstrom
import lsst.pex.policy  as pexPolicy

from astrometry.util.pyfits_utils import *
from numpy import array

from sourceset_boost_to_fits import *

class ducky(object):
    pass

class fakeExposure(ducky):
    def __init__(self):
        self.wcs = None
        self.width, self.height = None,None
        self.xy0 = None
        self.filterName = None

    def getMaskedImage(self):
        fakey = ducky()
        fakey.xy0 = self.xy0
        fakey.getXY0 = lambda(x): x.xy0
        return fakey

    def getFilter(self):
        fakey = ducky()
        fakey.filterName = self.filterName
        fakey.getName = lambda(x): x.filterName
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
          W=None, H=None, xy0=None, filtername=None):

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
        policy = pexPolicy.Policy.createPolicy('') #'''#<?cfg paf policy?>''')

    (matchList,wcs) = measAstrom.determineWcs(policy, exposure, sourceset)


if __name__ == '__main__':
    parser = OptionParser(usage='%prog [options] <*.boost SourceSets>')
    # or *.fits of sourcesets>')
    parser.add_option('-W', '--width', dest='width', type='float', help='Image width (pixels)')
    parser.add_option('-H', '--height', dest='height', type='float', help='Image height (pixels)')
    parser.add_option('-f', '--filter', dest='filter', help='Filter name')
    #parser.add_option('-x', '--x-column', dest='xcol', help='X column name (for FITS inputs)')
    parser.set_defaults(width=None, height=None, filter=None)
    opt,args = parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        sys.exit(0)

    for fn in args:
        #if fn.endswith('.boost'):
        print 'Reading', fn
        ss = sourceset_read_boost(fn)
        #else:
        #    ss = 
        rerun(ss, W=opt.width, H=opt.height, filtername=opt.filter)

