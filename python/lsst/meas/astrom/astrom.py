import os

import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
#import lsst.meas.astrom as measAst

class Astrometry(object):
    import config
    ConfigClass = config.AstromConfig

    def __init__(self,
                 conf,
                 #andDir=None,
                 andConfig=None,
                 log=None,
                 logLevel=pexLog.Log.INFO):
        '''
        conf: an AstromConfig object.
        andConfig: an AstromNetDataConfig object
        log: a pexLogging.Log
        logLevel: if log is None, the log level to use
        '''
        self.config = conf
        if log is not None:
            self.log = log
        else:
            self.log = pexLog.Log(pexLog.Log.getDefaultLog(),
                                  'meas.astrom',
                                  logLevel)

        if andConfig is not None:
            self.andConfig = andConfig
        else:
            #dirnm = andDir
            #if dirnm is None:
            # ASSUME SETUP IN EUPS
            dirnm = os.environ.get('ASTROMETRY_NET_DATA_DIR')
            if dirnm is None:
                self.log.log(pexLog.Log.WARN, 'astrometry_net_data is not setup')
            else:
                fn = os.path.join(dirnm, 'metadata.paf')
                self.andConfig = pexConfig.Config.load(fn)

    def setAndConfig(self, andconfig):
        self.andConfig = andconfig

    def determineWcs(self,
                     sources,
                     exposure=None,
                     wcs=None,
                     imageSize=None,
                     radecCenter=None,
                     searchRadius=None,
                     pixelScale=None,
                     filterName=None,
                     doTrim=False,
                     forceImageSize=None):
        '''
        We dont really need an Exposure; we need:
          -an initial Wcs estimate;
          -the image size;
          -the filter
        (all of which are metadata of Exposure).

        We also need the estimated pixel scale, which again we can get
        from the initial Wcs, but it should be possible to specify it
        without needing a Wcs.

        Same with the estimated RA,Dec and search radius.

        filterName: string
        imageSize: (W,H) integer tuple/iterable
        '''
        if exposure is not None:
            if filterName is None:
                filterName = exposure.getFilter().getName()
            if imageSize is None:
                imageSize = (exposure.getWidth(), exposure.getHeight())
            
