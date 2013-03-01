import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

import os

import lsst.meas.astrom as measAstrom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.daf.base  as dafBase

import pyfits

def main():
    mydir = os.path.dirname(__file__)

    # Read sources
    fn = os.path.join(mydir, 'xy2710.fits')
    P = pyfits.open(fn)[1].data
    x = P['x']
    y = P['y']
    srcSchema = afwTable.SourceTable.makeMinimalSchema()
    key = srcSchema.addField("centroid", type="PointD")
    srcTable = afwTable.SourceTable.make(srcSchema)
    srcTable.defineCentroid("centroid")
    for xi,yi in zip(x,y):
        src = srcTable.makeRecord()
        src.set(key.getX(), xi)
        src.set(key.getY(), yi)

    x0,y0 = 335750, 223750
    W,H = 8500, 12500
        
    # Read original WCS
    fn = os.path.join(mydir, 't2710.wcs')
    #    print help(afwImage.Wcs.readFits)
    wcs0 = afwImage.Wcs.readFits(fn, 0)
    #print dir(afwImage.Wcs)
    #wcs0 = afwImage.Wcs(fn)
    print 'WCS0:', wcs0
        


if __name__ == '__main__':
    main()
    
