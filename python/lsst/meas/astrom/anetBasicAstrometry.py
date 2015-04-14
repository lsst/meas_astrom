from __future__ import absolute_import, division

import os
import math
import sys

import numpy as np # for isfinite()

import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
from lsst.pex.config import Field, RangeField, ListField
import lsst.pex.exceptions as pexExceptions
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.algorithms.utils as maUtils
from .loadAstrometryNetObjects import LoadAstrometryNetObjectsTask

from .astrometryNetDataConfig import AstrometryNetDataConfig
from . import sip as astromSip

__all__ = ["InitialAstrometry", "ANetBasicAstrometryConfig", "ANetBasicAstrometryTask"]

class InitialAstrometry(object):
    '''
    Object returned by Astrometry.determineWcs

    getWcs(): sipWcs or tanWcs
    getMatches(): sipMatches or tanMatches

    Other fields are:
    solveQa (PropertyList)
    tanWcs (Wcs)
    tanMatches (MatchList)
    sipWcs (Wcs)
    sipMatches (MatchList)
    matchMetadata (PropertyList)
    '''
    def __init__(self):
        self.tanWcs = None
        self.tanMatches = None
        self.sipWcs = None
        self.sipMatches = None
        self.matchMetadata = dafBase.PropertyList()
        self.solveQa = None

    def getMatches(self):
        '''
        Get "sipMatches" -- MatchList using the SIP WCS solution, if it
        exists, or "tanMatches" -- MatchList using the TAN WCS solution
        otherwise.
        '''
        return self.sipMatches or self.tanMatches

    def getWcs(self):
        '''
        Returns the SIP WCS, if one was found, or a TAN WCS
        '''
        return self.sipWcs or self.tanWcs

    matches = property(getMatches)
    wcs = property(getWcs)
    
    ### "Not very pythonic!" complains Paul.
    # Consider these methods deprecated; if you want these elements, just
    # .grab them.
    def getSipWcs(self):
        return self.sipWcs
    def getTanWcs(self):
        return self.tanWcs
    def getSipMatches(self):
        return self.sipMatches
    def getTanMatches(self):
        return self.tanMatches
    def getMatchMetadata(self):
        return self.matchMetadata
    def getSolveQaMetadata(self):
        return self.solveQa


class ANetBasicAstrometryConfig(LoadAstrometryNetObjectsTask.ConfigClass):

    maxCpuTime = RangeField(
        '''Maximum CPU time to spend solving, in seconds''',
        float,
        default=0., min=0.)

    matchThreshold = RangeField(
        '''Matching threshold for Astrometry.net solver (log-odds)''',
        float,
        default=math.log(1e12), min=math.log(1e6))

    maxStars = RangeField(
        '''Maximum number of stars to use in Astrometry.net solving''',
        int,
        default=50, min=10)

    useWcsPixelScale = Field(
        '''Use the pixel scale from the input exposure\'s WCS headers?''',
        bool,
        default=True)

    useWcsRaDecCenter = Field(
        dtype=bool, default=True,
        doc='''Use the RA,Dec center information from the input exposure\'s WCS headers?''')

    useWcsParity = Field(
        dtype=bool, default=True,
        doc='''Use the parity (flip / handedness) of the image from the input exposure\'s WCS headers?''')

    raDecSearchRadius = RangeField(
        '''When useWcsRaDecCenter=True, this is the radius, in degrees, around the RA,Dec center specified in the input exposure\'s WCS to search for a solution.''',
        float,
        default=1., min=0.)

    pixelScaleUncertainty = RangeField(
        '''Range of pixel scales, around the value in the WCS header, to search.  If the value of this field is X and the nominal scale is S, the range searched will be  S/X to S*X''',
        float,
        default = 1.1, min=1.001)

    # forceParity?
    # forcePixelScale?
    # forceRaDecCenter?
    # forcePixelScaleRange?
    # doTrim ?

    # forceImageSize = Field(
    #     tuple,
    #     '''Ignore the size of the input exposure and assume this
    #     image size instead.''',
    #     optional=True,
    #     check=lambda x: (len(x) == 2 and
    #                      type(x[0]) is int and type(x[1]) is int))

    catalogMatchDist = RangeField(
        #afwGeom.Angle,
        '''When matching image to reference catalog stars, how big should
        the matching radius be?''',
        float,
        default=1.,#* afwGeom.arcseconds,
        min=0.)

    cleaningParameter = RangeField(
        '''Sigma-clipping parameter in sip/cleanBadPoints.py''',
        float,
        default=3., min=0.)

    calculateSip = Field(
        '''Compute polynomial SIP distortion terms?''',
        bool,
        default=True)

    sipOrder = RangeField(
        '''Polynomial order of SIP distortion terms''',
        int,
        default=4, min=2)

    badFlags = ListField(
        doc = "List of flags which cause a source to be rejected as bad",
        dtype = str,
        default = [
                   ]
        )                      

    allFluxes = Field(dtype=bool, default=True, doc="Retrieve all available fluxes (and errors) from catalog?")


class ANetBasicAstrometryTask(object):
    ConfigClass = ANetBasicAstrometryConfig

    '''
    About Astrometry.net index files (astrometry_net_data):

    There are three components of an index file: a list of stars
    (stored as a star kd-tree), a list of quadrangles of stars ("quad
    file") and a list of the shapes ("codes") of those quadrangles,
    stored as a code kd-tree.

    Each index covers a region of the sky, defined by healpix nside
    and number, and a range of angular scales.  In LSST, we share the
    list of stars in a part of the sky between multiple indexes.  That
    is, the star kd-tree is shared between multiple indices (quads and
    code kd-trees).  In the astrometry.net code, this is called a
    "multiindex".

    It is possible to "unload" and "reload" multiindex (and index)
    objects.  When "unloaded", they consume no FILE or mmap resources.

    The multiindex object holds the star kd-tree and gives each index
    object it holds a pointer to it, so it is necessary to
    multiindex_reload_starkd() before reloading the indices it holds.
    The multiindex_unload() method, on the other hand, unloads its
    starkd and unloads each index it holds.
    '''

    class _LoadedMIndexes(object):
        def __init__(self, multiInds):
            self.multiInds = multiInds
        def __enter__(self):
            for mi in self.multiInds:
                #print 'Loading', mi.name
                mi.reload()
            return self.multiInds
        def __exit__(self, typ, val, trace):
            for mi in self.multiInds:
                #print 'Unloading', mi.name
                mi.unload()

    def __init__(self,
                 config,
                 andConfig=None,
                 log=None,                 
                 logLevel=pexLog.Log.INFO):
        '''
        conf: an ANetBasicAstrometryConfig object
        andConfig: an AstromNetDataConfig object
        log: a pexLogging.Log
        logLevel: if log is None, the log level to use
        '''
        self.config = config
        self.refObjLoader = LoadAstrometryNetObjectsTask(self.config, andConfig=andConfig)
        if log is not None:
            self.log = log
        else:
            self.log = pexLog.Log(pexLog.Log.getDefaultLog(),
                                  'meas.astrom',
                                  logLevel)
        if andConfig is None:
            # ASSUME SETUP IN EUPS
            dirnm = os.environ.get('ASTROMETRY_NET_DATA_DIR')
            if dirnm is None:
                raise RuntimeError("astrometry_net_data is not setup")

            andConfig = AstrometryNetDataConfig()
            fn = os.path.join(dirnm, 'andConfig.py')
            if not os.path.exists(fn):
                raise RuntimeError('astrometry_net_data config file \"%s\" required but not found' % fn)
            andConfig.load(fn)

        self.andConfig = andConfig
        self._readIndexFiles()

    def _readIndexFiles(self):
        from . import astrometry_net as an
        # .multiInds: multi-index objects
        self.multiInds = []

        # merge indexFiles and multiIndexFiles; we'll treat both as
        # multiindex for simplicity.
        mifiles = ([(True,[fn,fn]) for fn  in self.andConfig.indexFiles] +
                     [(False,fns)    for fns in self.andConfig.multiIndexFiles])

        nMissing = 0
        for single,fns in mifiles:
            # First filename in "fns" is star kdtree, the rest are index files.
            fn = fns[0]
            if single:
                self.log.log(self.log.DEBUG, 'Adding index file %s' % fns[0])
            else:
                self.log.log(self.log.DEBUG, 'Adding multiindex files %s' % str(fns))
            fn2 = self._getIndexPath(fn)
            if fn2 is None:
                if single:
                    self.log.logdebug('Unable to find index file %s' % fn)
                else:
                    self.log.logdebug('Unable to find star part of multiindex file %s' % fn)
                nMissing += 1
                continue
            fn = fn2
            self.log.log(self.log.DEBUG, 'Path: %s' % fn)

            mi = an.multiindex_new(fn)
            if mi is None:
                raise RuntimeError('Failed to read stars from multiindex filename "%s"' % fn)
            for i,fn in enumerate(fns[1:]):
                self.log.log(self.log.DEBUG, 'Reading index from multiindex file "%s"' % fn)
                fn2 = self._getIndexPath(fn)
                if fn2 is None:
                    self.log.logdebug('Unable to find index part of multiindex file %s' % fn)
                    nMissing += 1
                    continue
                fn = fn2
                self.log.log(self.log.DEBUG, 'Path: %s' % fn)
                if an.multiindex_add_index(mi, fn, an.INDEX_ONLY_LOAD_METADATA):
                    raise RuntimeError('Failed to read index from multiindex filename "%s"' % fn)
                ind = mi[i]
                self.log.log(self.log.DEBUG, '  index %i, hp %i (nside %i), nstars %i, nquads %i' %
                                (ind.indexid, ind.healpix, ind.hpnside, ind.nstars, ind.nquads))
            an.multiindex_unload_starkd(mi)
            self.multiInds.append(mi)

        if len(self.multiInds) == 0:
            self.log.warn('Unable to find any index files')
        elif nMissing > 0:
            self.log.warn('Unable to find %d index files' % nMissing)

    def _debug(self, s):
        self.log.log(self.log.DEBUG, s)
    def _warn(self, s):
        self.log.log(self.log.WARN, s)

    def memusage(self, prefix=''):
        # Not logging at DEBUG: do nothing
        if self.log.getThreshold() > pexLog.Log.DEBUG:
            return
        from astrometry.util.ttime import get_memusage
        mu = get_memusage()
        ss = []
        for key in ['VmPeak', 'VmSize', 'VmRSS', 'VmData']:
            if key in mu:
                ss.append(key + ': ' + ' '.join(mu[key]))
        if 'mmaps' in mu:
            ss.append('Mmaps: %i' % len(mu['mmaps']))
        if 'mmaps_total' in mu:
            ss.append('Mmaps: %i kB' % (mu['mmaps_total'] / 1024))
        self.log.logdebug(prefix + 'Memory: ' + ', '.join(ss))

    def setAndConfig(self, andconfig):
        self.andConfig = andconfig

    def _getImageParams(self, wcs, exposure, filterName=None, imageSize=None,
                        x0=None, y0=None):
        if exposure is not None:
            ex0,ey0 = exposure.getX0(), exposure.getY0()
            if x0 is None:
                x0 = ex0
            if y0 is None:
                y0 = ey0
            self._debug('Got exposure x0,y0 = %i,%i' % (ex0,ey0))
            if filterName is None:
                filterName = exposure.getFilter().getName()
                self._debug('Setting filterName = "%s" from exposure metadata' % str(filterName))
            if imageSize is None:
                imageSize = (exposure.getWidth(), exposure.getHeight())
                self._debug('Setting image size = (%i, %i) from exposure metadata' % (imageSize))
        if x0 is None:
            x0 = 0
        if y0 is None:
            y0 = 0
        imageSize = afwGeom.Extent2I(*imageSize)
        return filterName, imageSize, x0, y0

    def useKnownWcs(self, sourceCat, wcs=None, exposure=None, filterName=None, imageSize=None,
                    x0=None, y0=None):
        """
        Returns an InitialAstrometry object, just like determineWcs,
        but assuming the given input WCS is correct.

        This is enabled by the ANetBasicAstrometryConfig
        'forceKnownWcs' option.  If you are using that option, you
        probably also want to turn OFF 'calculateSip'.

        This involves searching for reference sources within the WCS
        area, and matching them to the given 'sourceCat'.  If
        'calculateSip' is set, we will try to compute a TAN-SIP
        distortion correction.

        sourceCat: list of detected sources in this image.
        wcs: your known WCS
        exposure: the exposure holding metadata for this image.
        filterName: string, filter name, eg "i"
        x0,y0: image origin / offset; these coordinates along with the
           "imageSize" determine the bounding-box in pixel coordinates of
           the image in question; this is used for finding reference sources
           in the image, among other things.
        
        You MUST provide a WCS, either by providing the 'wcs' kwarg
        (an lsst.image.Wcs object), or by providing the 'exposure' on
        which we will call 'exposure.getWcs()'.

        You MUST provide a filter name, either by providing the
        'filterName' kwarg (a string), or by setting the 'exposure';
        we will call 'exposure.getFilter().getName()'.

        You MUST provide the image size, either by providing the
        'imageSize' kwargs, an (W,H) tuple of ints, or by providing
        the 'exposure'; we will call 'exposure.getWidth()' and
        'exposure.getHeight()'.

        Note, when modifying this function, that it is also called by
        'determineWcs' (via 'determineWcs2'), since the steps are all
        the same.
        """
        # return value:
        astrom = InitialAstrometry()

        if wcs is None:
            if exposure is None:
                raise RuntimeError('useKnownWcs: need either "wcs=" or "exposure=" kwarg.')
            wcs = exposure.getWcs()
            if wcs is None:
                raise RuntimeError('useKnownWcs: wcs==None and exposure.getWcs()==None.')
                
        filterName,imageSize,x0,y0 = self._getImageParams(exposure=exposure, wcs=wcs,
                                                          imageSize=imageSize,
                                                          filterName=filterName,
                                                          x0=x0, y0=y0)
        bbox = afwGeom.Box2I(afwGeom.Point2I(x0, y0), imageSize)
        refCat = self.refObjLoader.loadPixelBox(
            bbox = bbox,
            wcs = wcs,
            filterName = filterName,
            calib = None,
        ).refCat
        catids = [src.getId() for src in refCat]
        uids = set(catids)
        self.log.logdebug('%i reference sources; %i unique IDs' % (len(catids), len(uids)))
        matches = self._getMatchList(sourceCat, refCat, wcs)
        uniq = set([sm.second.getId() for sm in matches])
        if len(matches) != len(uniq):
            self._warn(('The list of matched stars contains duplicate reference source IDs ' +
                        '(%i sources, %i unique ids)') % (len(matches), len(uniq)))
        if len(matches) == 0:
            self._warn('No matches found between input sources and reference catalogue.')
            return astrom

        self._debug('%i reference objects match input sources using input WCS' % (len(matches)))
        astrom.tanMatches = matches
        astrom.tanWcs = wcs
        
        srcids = [s.getId() for s in sourceCat]
        for m in matches:
            assert(m.second.getId() in srcids)
            assert(m.second in sourceCat)

        if self.config.calculateSip:
            sipwcs,matches = self._calculateSipTerms(wcs, refCat, sourceCat, matches, imageSize, x0=x0, y0=y0)
            if sipwcs == wcs:
                self._debug('Failed to find a SIP WCS better than the initial one.')
            else:
                self._debug('%i reference objects match input sources using SIP WCS' % (len(matches)))
                astrom.sipWcs = sipwcs
                astrom.sipMatches = matches
                
        W,H = imageSize
        wcs = astrom.getWcs()
        # _getMatchList() modifies the source list RA,Dec coordinates.
        # Here, we make them consistent with the WCS we are returning.
        for src in sourceCat:
            src.updateCoord(wcs)
        astrom.matchMetadata = _createMetadata(W, H, x0, y0, wcs, filterName)
        return astrom

    def determineWcs(self,
                     sourceCat,
                     exposure,
                     **kwargs):
        """
        Finds a WCS solution for the given 'sources' in the given
        'exposure', getting other parameters from config.

        Valid kwargs include:

        'radecCenter', an afw.coord.Coord giving the RA,Dec position
           of the center of the field.  This is used to limit the
           search done by Astrometry.net (to make it faster and load
           fewer index files, thereby using less memory).  Defaults to
           the RA,Dec center from the exposure's WCS; turn that off
           with the boolean kwarg 'useRaDecCenter' or config option
           'useWcsRaDecCenter'

        'useRaDecCenter', a boolean.  Don't use the RA,Dec center from
           the exposure's initial WCS.

        'searchRadius', in degrees, to search for a solution around
           the given 'radecCenter'; default from config option
           'raDecSearchRadius'.

        'useParity': parity is the 'flip' of the image.  Knowing it
           reduces the search space (hence time) for Astrometry.net.
           The parity can be computed from the exposure's WCS (the
           sign of the determinant of the CD matrix); this option
           controls whether we do that or force Astrometry.net to
           search both parities.  Default from config.useWcsParity.

        'pixelScale': afwGeom.Angle, estimate of the angle-per-pixel
           (ie, arcseconds per pixel).  Defaults to a value derived
           from the exposure's WCS.  If enabled, this value, plus or
           minus config.pixelScaleUncertainty, will be used to limit
           Astrometry.net's search.

        'usePixelScale': boolean.  Use the pixel scale to limit
           Astrometry.net's search?  Defaults to config.useWcsPixelScale.

        'filterName', a string, the filter name of this image.  Will
           be mapped through the 'filterMap' config dictionary to a
           column name in the astrometry_net_data index FITS files.
           Defaults to the exposure.getFilter() value.

        'imageSize', a tuple (W,H) of integers, the image size.
           Defaults to the exposure.get{Width,Height}() values.

        """
        assert(exposure is not None)

        margs = kwargs.copy()
        if not 'searchRadius' in margs:
            margs.update(searchRadius = self.config.raDecSearchRadius * afwGeom.degrees)
        if not 'usePixelScale' in margs:
            margs.update(usePixelScale = self.config.useWcsPixelScale)
        if not 'useRaDecCenter' in margs:
            margs.update(useRaDecCenter = self.config.useWcsRaDecCenter)
        if not 'useParity' in margs:
            margs.update(useParity = self.config.useWcsParity)
        margs.update(exposure=exposure)
        return self.determineWcs2(sourceCat, **margs)

    def determineWcs2(self, sourceCat, **kwargs):
        '''
        Get a blind astrometric solution for the given catalog of sources.

        We need:
          -the image size;
          -the filter

        And if available, we can use:
          -an initial Wcs estimate;
             --> RA,Dec center
             --> pixel scale
             --> "parity"
             
        (all of which are metadata of Exposure).

        filterName: string
        imageSize: (W,H) integer tuple/iterable
        pixelScale: afwGeom::Angle per pixel.
        radecCenter: afwCoord::Coord
        '''
        wcs,qa = self.getBlindWcsSolution(sourceCat, **kwargs)
        kw = {}
        # Keys passed to useKnownWcs
        for key in ['exposure', 'filterName', 'imageSize', 'x0', 'y0']:
            if key in kwargs:
                kw[key] = kwargs[key]
        astrom = self.useKnownWcs(sourceCat, wcs=wcs, **kw)
        astrom.solveQa = qa
        astrom.tanWcs = wcs
        return astrom

    def getBlindWcsSolution(self, sourceCat, 
                            exposure=None,
                            wcs=None,
                            imageSize=None,
                            x0=None, y0=None,
                            radecCenter=None,
                            searchRadius=None,
                            pixelScale=None,
                            filterName=None,
                            doTrim=False,
                            usePixelScale=True,
                            useRaDecCenter=True,
                            useParity=True,
                            searchRadiusScale=2.):
        if not useRaDecCenter and radecCenter is not None:
            raise RuntimeError('radecCenter is set, but useRaDecCenter is False.  Make up your mind!')
        if not usePixelScale and pixelScale is not None:
            raise RuntimeError('pixelScale is set, but usePixelScale is False.  Make up your mind!')

        filterName,imageSize,x0,y0 = self._getImageParams(exposure=exposure, wcs=wcs,
                                                          imageSize=imageSize,
                                                          filterName=filterName,
                                                          x0=x0, y0=y0)

        if exposure is not None:
            if wcs is None:
                wcs = exposure.getWcs()
                self._debug('Setting initial WCS estimate from exposure metadata')

        if imageSize is None:
            # Could guess from the extent of the Sources...
            raise RuntimeError('Image size must be specified by passing "exposure" or "imageSize"')

        W,H = imageSize
        xc, yc = W/2. + 0.5 + x0, H/2. + 0.5 + y0
        parity = None

        if wcs is not None:
            if pixelScale is None:
                if usePixelScale:
                    pixelScale = wcs.pixelScale()
                    self._debug('Setting pixel scale estimate = %.3f from given WCS estimate' %
                                (pixelScale.asArcseconds()))

            if radecCenter is None:
                if useRaDecCenter:
                    radecCenter = wcs.pixelToSky(xc, yc)
                    self._debug(('Setting RA,Dec center estimate = (%.3f, %.3f) from given WCS '
                                 + 'estimate, using pixel center = (%.1f, %.1f)') %
                                (radecCenter.getLongitude().asDegrees(),
                                 radecCenter.getLatitude().asDegrees(), xc, yc))

            if searchRadius is None:
                if useRaDecCenter:
                    assert(pixelScale is not None)
                    searchRadius = (pixelScale * math.hypot(W,H)/2. *
                                    searchRadiusScale)
                    self._debug(('Using RA,Dec search radius = %.3f deg, from pixel scale, '
                                 + 'image size, and searchRadiusScale = %g') %
                                (searchRadius, searchRadiusScale))
            if useParity:
                parity = wcs.isFlipped()
                self._debug('Using parity = %s' % (parity and 'True' or 'False'))

        if doTrim:
            n = len(sourceCat)
            if exposure is not None:
                bbox = afwGeom.Box2D(exposure.getMaskedImage().getBBox())
            else:
                # CHECK -- half-pixel issues here?
                bbox = afwGeom.Box2D(afwGeom.Point2D(0.,0.), afwGeom.Point2D(W, H))
            sourceCat = self._trimBadPoints(sourceCat, bbox)
            self._debug("Trimming: kept %i of %i sources" % (n, len(sourceCat)))

        wcs,qa = self._solve(sourceCat, wcs, imageSize, pixelScale, radecCenter, searchRadius, parity,
                             filterName, xy0=(x0,y0))
        if wcs is None:
            raise RuntimeError("Unable to match sources with catalog.")
        self.log.info('Got astrometric solution from Astrometry.net')

        rdc = wcs.pixelToSky(xc, yc)
        self._debug('New WCS says image center pixel (%.1f, %.1f) -> RA,Dec (%.3f, %.3f)' %
                    (xc, yc, rdc.getLongitude().asDegrees(), rdc.getLatitude().asDegrees()))
        return wcs,qa

    def getSipWcsFromWcs(self, wcs, imageSize, x0=0, y0=0, ngrid=20,
                         linearizeAtCenter=True):
        '''
        This function allows one to get a TAN-SIP WCS, starting from
        an existing WCS.  It uses your WCS to compute a fake grid of
        corresponding "stars" in pixel and sky coords, and feeds that
        to the regular SIP code.

        linearizeCenter: if True, get a linear approximation of the input
          WCS at the image center and use that as the TAN initialization for
          the TAN-SIP solution.  You probably want this if your WCS has its
          CRPIX outside the image bounding box.
          
        '''
        # Ugh, build src and ref tables
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        key = srcSchema.addField("centroid", type="PointD")
        srcTable = afwTable.SourceTable.make(srcSchema)
        srcTable.defineCentroid("centroid")
        srcs = srcTable
        refs = afwTable.SimpleTable.make(afwTable.SimpleTable.makeMinimalSchema())
        cref = []
        csrc = []
        (W,H) = imageSize
        for xx in np.linspace(0., W, ngrid):
            for yy in np.linspace(0, H, ngrid):
                src = srcs.makeRecord()
                src.set(key.getX(), x0 + xx)
                src.set(key.getY(), y0 + yy)
                csrc.append(src)
                rd = wcs.pixelToSky(afwGeom.Point2D(xx + x0, yy + y0))
                ref = refs.makeRecord()
                ref.setCoord(rd)
                cref.append(ref)

        if linearizeAtCenter:
            # Linearize the original WCS around the image center to create a
            # TAN WCS.
            # Reference pixel in LSST coords
            crpix = afwGeom.Point2D(x0 + W/2. - 0.5, y0 + H/2. - 0.5)
            crval = wcs.pixelToSky(crpix)
            crval = crval.getPosition(afwGeom.degrees)
            # Linearize *AT* crval to get effective CD at crval.
            # (we use the default skyUnit of degrees as per WCS standard)
            aff = wcs.linearizePixelToSky(crval)
            cd = aff.getLinear().getMatrix()
            wcs = afwImage.Wcs(crval, crpix, cd)
                
        return self.getSipWcsFromCorrespondences(wcs, cref, csrc, (W,H),
                                                 x0=x0, y0=y0)

    
    def getSipWcsFromCorrespondences(self, origWcs, refCat, sourceCat, imageSize,
                                     x0=0, y0=0):
        '''
        Produces a SIP solution given a list of known correspondences.
        Unlike _calculateSipTerms, this does not iterate the solution;
        it assumes you have given it a good sets of corresponding stars.

        NOTE that "refCat" and "sourceCat" are assumed to be the same length;
        entries "refCat[i]" and "sourceCat[i]" are assumed to be correspondences.

        origWcs: the WCS to linearize in order to get the TAN part of the
           TAN-SIP WCS.

        refCat: reference source catalog

        sources: image sources

        imageSize, x0, y0: these describe the bounding-box of the image,
            which is used when computing reverse SIP polynomials.

        '''
        sipOrder = self.config.sipOrder
        bbox = afwGeom.Box2I(afwGeom.Point2I(x0,y0),
                             afwGeom.Extent2I(imageSize[0], imageSize[1]))
        matches = []
        for ci,si in zip(refCat, sourceCat):
            matches.append(afwTable.ReferenceMatch(ci, si, 0.))

        sipObject = astromSip.makeCreateWcsWithSip(matches, origWcs, sipOrder, bbox)
        return sipObject.getNewWcs()
    
    def _calculateSipTerms(self, origWcs, refCat, sourceCat, matches, imageSize,
                           x0=0, y0=0):
        '''
        Iteratively calculate SIP distortions and regenerate matches based on improved WCS.

        origWcs: original WCS object, probably (but not necessarily) a TAN WCS;
           this is used to set the baseline when determining whether a SIP
           solution is any better; it will be returned if no better SIP solution
           can be found.

        matches: list of supposedly matched sources, using the "origWcs".

        refCat: reference source catalog

        sourceCat: sources in the image to be solved

        imageSize, x0, y0: these determine the bounding-box of the image,
           which is used when finding reverse SIP coefficients.
        '''
        sipOrder = self.config.sipOrder
        wcs = origWcs
        bbox = afwGeom.Box2I(afwGeom.Point2I(x0,y0),
                             afwGeom.Extent2I(imageSize[0], imageSize[1]))

        i=0
        lastScatPix = None
        while True:
            try:
                sipObject = astromSip.makeCreateWcsWithSip(matches, wcs, sipOrder, bbox)
                if lastScatPix is None:
                    lastScatPix = sipObject.getLinearScatterInPixels()
                proposedWcs = sipObject.getNewWcs()
                scatPix = sipObject.getScatterInPixels()
                self.plotSolution(matches, proposedWcs, imageSize)
            except pexExceptions.Exception as e:
                self._warn('Failed to calculate distortion terms. Error: ' + str(e))
                break

            matchSize = len(matches)
            # use new WCS to get new matchlist.
            proposedMatchlist = self._getMatchList(sourceCat, refCat, proposedWcs)

            self._debug('SIP iteration %i: %i objects match.  Median scatter is %g arcsec = %g pixels (vs previous: %i matches, %g pixels)' %
                        (i, len(proposedMatchlist), sipObject.getScatterOnSky().asArcseconds(), scatPix, matchSize, lastScatPix))
            #self._debug('Proposed WCS: ' + proposedWcs.getFitsMetadata().toString())
            # Hack convergence tests
            if len(proposedMatchlist) < matchSize:
                break
            if len(proposedMatchlist) == matchSize and scatPix >= lastScatPix:
                break

            wcs = proposedWcs
            matches = proposedMatchlist
            lastScatPix = scatPix
            matchSize = len(matches)
            i += 1

        return wcs, matches

    def plotSolution(self, matches, wcs, imageSize):
        """Plot the solution, when debugging is turned on.

        @param matches   The list of matches
        @param wcs         The Wcs
        @param imageSize   2-tuple with the image size (W,H)
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display 
        if not display:
            return

        try:
            import matplotlib.pyplot as plt
            import numpy
        except ImportError:
            print >> sys.stderr, "Unable to import matplotlib"
            return

        fig = plt.figure(1)
        fig.clf()
        try:
            fig.canvas._tkcanvas._root().lift() # == Tk's raise, but raise is a python reserved word
        except:                                 # protect against API changes
            pass

        num = len(matches)
        x = numpy.zeros(num)
        y = numpy.zeros(num)
        dx = numpy.zeros(num)
        dy = numpy.zeros(num)
        for i, m in enumerate(matches):
            x[i] = m.second.getX()
            y[i] = m.second.getY()
            pixel = wcs.skyToPixel(m.first.getCoord())
            dx[i] = x[i] - pixel.getX()
            dy[i] = y[i] - pixel.getY()

        subplots = maUtils.makeSubplots(fig, 2, 2, xgutter=0.1, ygutter=0.1, pygutter=0.04)

        def plotNext(x, y, xLabel, yLabel, xMax):
            ax = subplots.next()
            ax.set_autoscalex_on(False)
            ax.set_xbound(lower=0, upper=xMax)
            ax.scatter(x, y)
            ax.set_xlabel(xLabel)
            ax.set_ylabel(yLabel)
            ax.axhline(0.0)

        plotNext(x, dx, "x", "dx", imageSize[0])
        plotNext(x, dy, "x", "dy", imageSize[0])
        plotNext(y, dx, "y", "dx", imageSize[1])
        plotNext(y, dy, "y", "dy", imageSize[1])

        fig.show()

        while True:
            try:
                reply = raw_input("Pausing for inspection, enter to continue... [hpQ] ").strip()
            except EOFError:
                reply = "n"

            reply = reply.split()
            if reply:
                reply = reply[0]
            else:
                reply = ""

            if reply in ("", "h", "p", "Q"):
                if reply == "h":
                    print "h[elp] p[db] Q[uit]"
                    continue
                elif reply == "p":
                    import pdb; pdb.set_trace() 
                elif reply == "Q":
                    sys.exit(1)
                break

    def _getMatchList(self, sourceCat, refCat, wcs):
        dist = self.config.catalogMatchDist * afwGeom.arcseconds
        clean = self.config.cleaningParameter
        matcher = astromSip.MatchSrcToCatalogue(refCat, sourceCat, wcs, dist)
        matches = matcher.getMatches()
        if matches is None:
            # Produce debugging stats...
            X = [src.getX() for src in sourceCat]
            Y = [src.getY() for src in sourceCat]
            R1 = [src.getRa().asDegrees() for src in sourceCat]
            D1 = [src.getDec().asDegrees() for src in sourceCat]
            R2 = [src.getRa().asDegrees() for src in refCat]
            D2 = [src.getDec().asDegrees() for src in refCat]
            # for src in sourceCat:
            #self._debug("source: x,y (%.1f, %.1f), RA,Dec (%.3f, %.3f)" %
            #(src.getX(), src.getY(), src.getRa().asDegrees(), src.getDec().asDegrees()))
            #for src in refCat:
            #self._debug("ref: RA,Dec (%.3f, %.3f)" %
            #(src.getRa().asDegrees(), src.getDec().asDegrees()))
            self.loginfo('_getMatchList: %i sources, %i reference sources' % (len(sourceCat), len(refCat)))
            if len(sourceCat):
                self.loginfo('Source range: x [%.1f, %.1f], y [%.1f, %.1f], RA [%.3f, %.3f], Dec [%.3f, %.3f]' %
                             (min(X), max(X), min(Y), max(Y), min(R1), max(R1), min(D1), max(D1)))
            if len(refCat):
                self.loginfo('Reference range: RA [%.3f, %.3f], Dec [%.3f, %.3f]' %
                             (min(R2), max(R2), min(D2), max(D2)))
            raise RuntimeError('No matches found between image and catalogue')
        matches = astromSip.cleanBadPoints.clean(matches, wcs, nsigma=clean)
        return matches

    def getColumnName(self, filterName, columnMap, default=None):
        '''
        Returns the column name in the astrometry_net_data index file that will be used
        for the given filter name.

        @param filterName   Name of filter used in exposure
        @param columnMap    Dict that maps filter names to column names
        @param default      Default column name
        '''
        filterName = self.config.filterMap.get(filterName, filterName) # Exposure filter --> desired filter
        try:
            return columnMap[filterName] # Desired filter --> a_n_d column name
        except KeyError:
            self.log.warn("No column in configuration for filter '%s'; using default '%s'" %
                          (filterName, default))
            return default

    def _solve(self, sourceCat, wcs, imageSize, pixelScale, radecCenter,
               searchRadius, parity, filterName=None, xy0=None):
        solver = self._getSolver()

        x0,y0 = 0,0
        if xy0 is not None:
            x0,y0 = xy0

        # select sources with valid x,y, flux
        xybb = afwGeom.Box2D()
        goodsources = afwTable.SourceCatalog(sourceCat.table)
        badkeys = [goodsources.getSchema().find(name).key for name in self.config.badFlags]

        for s in sourceCat:
            if np.isfinite(s.getX()) and np.isfinite(s.getY()) and np.isfinite(s.getPsfFlux()) and self._isGoodSource(s, badkeys) :
                goodsources.append(s)
                xybb.include(afwGeom.Point2D(s.getX() - x0, s.getY() - y0))
        self.log.info("Number of selected sources for astrometry : %d" %(len(goodsources)))
        if len(goodsources) < len(sourceCat):
            self.log.logdebug('Keeping %i of %i sources with finite X,Y positions and PSF flux' %
                              (len(goodsources), len(sourceCat)))
        self._debug(('Feeding sources in range x=[%.1f, %.1f], y=[%.1f, %.1f] ' +
                     '(after subtracting x0,y0 = %.1f,%.1f) to Astrometry.net') %
                    (xybb.getMinX(), xybb.getMaxX(), xybb.getMinY(), xybb.getMaxY(), x0, y0))
        # setStars sorts them by PSF flux.
        solver.setStars(goodsources, x0, y0)
        solver.setMaxStars(self.config.maxStars)
        solver.setImageSize(*imageSize)
        solver.setMatchThreshold(self.config.matchThreshold)
        raDecRadius = None
        if radecCenter is not None:
            raDecRadius = (radecCenter.getLongitude().asDegrees(),
                        radecCenter.getLatitude().asDegrees(),
                        searchRadius.asDegrees())
            solver.setRaDecRadius(*raDecRadius)
            self.log.logdebug('Searching for match around RA,Dec = (%g, %g) with radius %g deg' %
                              raDecRadius)

        if pixelScale is not None:
            dscale = self.config.pixelScaleUncertainty
            scale = pixelScale.asArcseconds()
            lo = scale / dscale
            hi = scale * dscale
            solver.setPixelScaleRange(lo, hi)
            self.log.logdebug('Searching for matches with pixel scale = %g +- %g %% -> range [%g, %g] arcsec/pix' %
                              (scale, 100.*(dscale-1.), lo, hi))

        if parity is not None:
            solver.setParity(parity)
            self.log.logdebug('Searching for match with parity = ' + str(parity))

        # Find and load index files within RA,Dec range and scale range.
        if radecCenter is not None:
            multiInds = self._getMIndexesWithinRange(radecCenter, searchRadius)
        else:
            multiInds = self.multiInds
        qlo,qhi = solver.getQuadSizeLow(), solver.getQuadSizeHigh()

        toload_multiInds = set()
        toload_inds = []
        for mi in multiInds:
            for i in range(len(mi)):
                ind = mi[i]
                if not ind.overlapsScaleRange(qlo, qhi):
                    continue
                toload_multiInds.add(mi)
                toload_inds.append(ind)

        with ANetBasicAstrometryTask._LoadedMIndexes(toload_multiInds):
            solver.addIndices(toload_inds)
            self.memusage('Index files loaded: ')

            cpulimit = self.config.maxCpuTime
            solver.run(cpulimit)

            self.memusage('Solving finished: ')

        self.memusage('Index files unloaded: ')

        if solver.didSolve():
            self.log.logdebug('Solved!')
            wcs = solver.getWcs()
            self.log.logdebug('WCS: %s' % wcs.getFitsMetadata().toString())

            if x0 != 0 or y0 != 0:
                wcs.shiftReferencePixel(x0, y0)
                self.log.logdebug('After shifting reference pixel by x0,y0 = (%i,%i), WCS is: %s' %
                                  (x0, y0, wcs.getFitsMetadata().toString()))

        else:
            self.log.warn('Did not get an astrometric solution from Astrometry.net')
            wcs = None
            # Gather debugging info...

            # -are there any reference stars in the proposed search area?
            # log the number found and discard the results
            if radecCenter is not None:
                self.refObjLoader.loadSkyCircle(radecCenter, searchRadius, filterName)

        qa = solver.getSolveStats()
        self.log.logdebug('qa: %s' % qa.toString())
        return wcs, qa

    def _isGoodSource(self, candsource, keys):
        for k in keys:
            if candsource.get(k):
                return False
        return True

    def _getIndexPath(self, fn):
        if os.path.isabs(fn):
            return fn
        andir = os.getenv('ASTROMETRY_NET_DATA_DIR')
        if andir is not None:
            fn2 = os.path.join(andir, fn)
            if os.path.exists(fn2):
                return fn2

        if os.path.exists(fn):
            return os.path.abspath(fn)
        else:
            return None

    def _getMIndexesWithinRange(self, ctrCoord, radius):
        '''
        ra,dec,radius: [deg], spatial cut based on the healpix of the index

        Returns list of multiindex objects within range.
        '''
        good = []
        raDeg = ctrCoord.getLongitude().asDegrees()
        decDeg = ctrCoord.getLatitude().asDegrees()
        radiusDeg = radius.asDegrees()
        for mi in self.multiInds:
            if mi.isWithinRange(raDeg, decDeg, radiusDeg):
                good.append(mi)
        return good

    def _getSolver(self):
        from . import astrometry_net as an
        solver = an.solver_new()
        # HACK, set huge default pixel scale range.
        lo,hi = 0.01, 3600.
        solver.setPixelScaleRange(lo, hi)
        return solver

    @staticmethod
    def _trimBadPoints(sourceCat, bbox, wcs=None):
        '''Remove elements from catalog whose xy positions are not within the given bbox.

        sourceCat:  a Catalog of SimpleRecord or SourceRecord objects
        bbox: an afwImage.Box2D
        wcs:  if not None, will be used to compute the xy positions on-the-fly;
              this is required when sources actually contains SimpleRecords.
        
        Returns:
        a list of Source objects with xAstrom,yAstrom within the bbox.
        '''
        keep = type(sourceCat)(sourceCat.table)
        for s in sourceCat:
            point = s.getCentroid() if wcs is None else wcs.skyToPixel(s.getCoord())
            if bbox.contains(point):
                keep.append(s)
        return keep

    def joinMatchListWithCatalog(self, packedMatches, sourceCat):
        '''
        This function is required to reconstitute a ReferenceMatchVector after being
        unpersisted.  The persisted form of a ReferenceMatchVector is the 
        normalized Catalog of IDs produced by afw.table.packMatches(), with the result of 
        InitialAstrometry.getMatchMetadata() in the associated tables\' metadata.

        The "live" form of a matchlist has links to
        the real record objects that are matched; it is "denormalized".
        This function takes a normalized match catalog, along with the catalog of
        sources to which the match catalog refers.  It fetches the reference
        sources that are within range, and then denormalizes the matches
        -- sets the "matches[*].first" and "matches[*].second" entries
        to point to the sources in the "sourceCat" argument, and to the
        reference sources fetched from the astrometry_net_data files.
    
        @param[in] packedMatches  Unpersisted match list (an lsst.afw.table.BaseCatalog).
                                  packedMatches.table.getMetadata() must contain the
                                  values from InitialAstrometry.getMatchMetadata()
        @param[in,out] sourceCat  Source catalog used for the 'second' side of the matches
                                  (an lsst.afw.table.SourceCatalog).  As a side effect,
                                  the catalog will be sorted by ID.
        
        @return An lsst.afw.table.ReferenceMatchVector of denormalized matches.
        '''
        matchmeta = packedMatches.table.getMetadata()
        version = matchmeta.getInt('SMATCHV')
        if version != 1:
            raise ValueError('SourceMatchVector version number is %i, not 1.' % version)
        filterName = matchmeta.getString('FILTER').strip()
        ctrCoord = afwCoord.IcrsCoord(
            matchmeta.getDouble('RA') * afwGeom.degrees,
            matchmeta.getDouble('DEC') * afwGeom.degrees,
        )        
        rad = matchmeta.getDouble('RADIUS') * afwGeom.degrees
        refCat = self.refObjLoader.loadSkyCircle(ctrCoord, rad, filterName).refCat
        refCat.sort()
        sourceCat.sort()
        return afwTable.unpackMatches(packedMatches, refCat, sourceCat)


def _createMetadata(width, height, x0, y0, wcs, filterName):
    """
    Create match metadata entries required for regenerating the catalog

    @param width Width of the image (pixels)
    @param height Height of the image (pixels)
    @param x0 x offset of image origin from parent (pixels)
    @param y0 y offset of image origin from parent (pixels)
    @param filterName Name of filter, used for magnitudes
    @return Metadata
    """
    meta = dafBase.PropertyList()

    # cache: field center and size.
    cx,cy = x0 + 0.5 + width/2., y0 + 0.5 + height/2.
    radec = wcs.pixelToSky(cx, cy).toIcrs()
    meta.add('RA', radec.getRa().asDegrees(), 'field center in degrees')
    meta.add('DEC', radec.getDec().asDegrees(), 'field center in degrees')
    imgSize = wcs.pixelScale() * math.hypot(width, height)/2.
    meta.add('RADIUS', imgSize.asDegrees(),
             'field radius in degrees, approximate')
    meta.add('SMATCHV', 1, 'SourceMatchVector version number')
    if filterName is not None:
        meta.add('FILTER', filterName, 'LSST filter name for tagalong data')
    return meta

def readMatches(butler, dataId, sourcesName='icSrc', matchesName='icMatch', config=ANetBasicAstrometryConfig(),
                sourcesFlags=afwTable.SOURCE_IO_NO_FOOTPRINTS):
    """Read matches, sources and catalogue; combine.

    @param butler Data butler
    @param dataId Data identifier for butler
    @param sourcesName Name for sources from butler
    @param matchesName Name for matches from butler
    @param sourcesFlags Flags to pass for source retrieval
    @returns Matches
    """
    sourceCat = butler.get(sourcesName, dataId, flags=sourcesFlags)
    packedMatches = butler.get(matchesName, dataId)
    astrom = ANetBasicAstrometryTask(config)
    return astrom.joinMatchListWithCatalog(packedMatches, sourceCat)
