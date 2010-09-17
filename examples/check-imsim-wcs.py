from rerun_wcs import rerun
from sourceset_boost_to_fits import *

import lsst.pex.logging as pexLog

if __name__ == '__main__':
    W,H = 4072,4000

    runset = ('v85470977-fr/R22/S00.boost',
              'v85471033-fg/R22/S00.boost',
              'v85501865-fi/R22/S00.boost',
              'v85505309-fy/R22/S00.boost',
              'v85518116-fu/R22/S00.boost',
              'v85563915-fz/R22/S00.boost')

    verb=False
    
    level = pexLog.Log.DEBUG if verb else pexLog.Log.INFO
    log = pexLog.Log(pexLog.Log.getDefaultLog(), "rerun-wcs", level);

    for x in runset:
        prefix = x.replace('/','-').replace('.boost','')
        filt = (x.split('-')[1])[1]

        fn = 'icSrc/' + x
        print 'Reading boost-format', fn
        ss = sourceset_read_boost(fn)

        rerun(ss, W=W, H=H, filtername=filt, plotprefix=prefix+'-',
              outwcsfn=prefix+'.wcs', log=log)

