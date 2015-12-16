from __future__ import absolute_import, division

import math
import sys

import numpy as np # for isfinite()

import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
from lsst.pex.config import Field, RangeField, ListField
import lsst.pex.exceptions as pexExceptions
import lsst.pipe.base as pipeBase
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.algorithms.utils as maUtils
from .loadAstrometryNetObjects import LoadAstrometryNetObjectsTask, LoadMultiIndexes
from .display import displayAstrometry
from .astromLib import makeMatchStatisticsInRadians

from . import sip as astromSip

__all__ = ["InitialAstrometry", "ANetBasicAstrometryConfig", "ANetBasicAstrometryTask"]

class InitialAstrometry(object):
    """
    Object returned by Astrometry.determineWcs

    getWcs(): sipWcs or tanWcs
    getMatches(): sipMatches or tanMatches

    Other fields are:
    solveQa (PropertyList)
    tanWcs (Wcs)
    tanMatches (MatchList)
    sipWcs (Wcs)
    sipMatches (MatchList)
    refCat (lsst::afw::table::SimpleCatalog)
    matchMeta (PropertyList)
    """
    def __init__(self):
        self.tanWcs = None
        self.tanMatches = None
        self.sipWcs = None
        self.sipMatches = None
        self.refCat = None
        self.matchMeta = dafBase.PropertyList()
        self.solveQa = None

    def getMatches(self):
        """
        Get "sipMatches" -- MatchList using the SIP WCS solution, if it
        exists, or "tanMatches" -- MatchList using the TAN WCS solution
        otherwise.
        """
        return self.sipMatches or self.tanMatches

    def getWcs(self):
        """
        Returns the SIP WCS, if one was found, or a TAN WCS
        """
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
        return self.matchMeta
    def getSolveQaMetadata(self):
        return self.solveQa


class ANetBasicAstrometryConfig(LoadAstrometryNetObjectsTask.ConfigClass):

    maxCpuTime = RangeField(
        doc = "Maximum CPU time to spend solving, in seconds",
        dtype = float,
        default = 0.,
        min = 0.,
    )
    matchThreshold = RangeField(
        doc = "Matching threshold for Astrometry.net solver (log-odds)",
        dtype = float,
        default=math.log(1e12),
        min=math.log(1e6),
    )
    maxStars = RangeField(
        doc = "Maximum number of stars to use in Astrometry.net solving",
        dtype = int,
        default=50,
        min=10,
    )
    useWcsPixelScale = Field(
        doc = "Use the pixel scale from the input exposure\'s WCS headers?",
        dtype = bool,
        default = True,
    )
    useWcsRaDecCenter = Field(
        doc="Use the RA,Dec center information from the input exposure\'s WCS headers?",
        dtype=bool,
        default=True,
    )
    useWcsParity = Field(
        doc="Use the parity (flip / handedness) of the image from the input exposure\'s WCS headers?",
        dtype=bool,
        default=True,
    )
    raDecSearchRadius = RangeField(
        doc = "When useWcsRaDecCenter=True, this is the radius, in degrees, around the RA,Dec center " +
            "specified in the input exposure\'s WCS to search for a solution.",
        dtype = float,
        default=1.0,
        min=0.0,
    )
    pixelScaleUncertainty = RangeField(
        doc = "Range of pixel scales, around the value in the WCS header, to search. " +
            "If the value of this field is X and the nominal scale is S, " +
            "the range searched will be  S/X to S*X",
        dtype = float,
        default = 1.1,
        min=1.001,
    )
    catalogMatchDist = RangeField(
        doc = "Matching radius (arcsec) for matching sources to reference objects",
        dtype = float,
        default=1.0,
        min=0.0,
    )
    cleaningParameter = RangeField(
        doc = "Sigma-clipping parameter in sip/cleanBadPoints.py",
        dtype = float,
        default=3.0,
        min=0.0,
    )
    calculateSip = Field(
        doc = "Compute polynomial SIP distortion terms?",
        dtype = bool,
        default=True,
    )
    sipOrder = RangeField(
        doc = "Polynomial order of SIP distortion terms",
        dtype = int,
        default=4,
        min=2,
    )
    badFlags = ListField(
        doc = "List of flags which cause a source to be rejected as bad",
        dtype = str,
        default = [
            "slot_Centroid_flag", # bad centroids
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_saturated",
            "base_PixelFlags_flag_crCenter", # cosmic rays
        ],
    )
    allFluxes = Field(
        doc = "Retrieve all available fluxes (and errors) from catalog?",
        dtype = bool,
        default = True,
    )
    maxIter = RangeField(
        doc = "maximum number of iterations of match sources and fit WCS" +
            "ignored if not fitting a WCS",
        dtype = int,
        default = 5,
        min = 1,
    )
    matchDistanceSigma = RangeField(
        doc = "The match and fit loop stops when maxMatchDist minimized: "
            " maxMatchDist = meanMatchDist + matchDistanceSigma*stdDevMatchDistance " +
            " (where the mean and std dev are computed using outlier rejection);" +
            " ignored if not fitting a WCS",
        dtype = float,
        default = 2,
        min = 0,
    )


class ANetBasicAstrometryTask(pipeBase.Task):
    """!Basic implemeentation of the astrometry.net astrometrical fitter

    A higher-level class ANetAstrometryTask takes care of dealing with the fact
    that the initial WCS is probably only a pure TAN SIP, yet we may have
    significant distortion and a good estimate for that distortion.


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
    """
    ConfigClass = ANetBasicAstrometryConfig
    _DefaultName = "aNetBasicAstrometry"
    def __init__(self,
                 config,
                 andConfig=None,
                 **kwargs):
        """!Construct an ANetBasicAstrometryTask

        @param[in] config  configuration (an instance of self.ConfigClass)
        @param[in] andConfig  astrometry.net data config (an instance of AstromNetDataConfig, or None);
            if None then use andConfig.py in the astrometry_net_data product (which must be setup)
        @param[in] kwargs  additional keyword arguments for pipe_base Task.\_\_init\_\_

        @throw RuntimeError if andConfig is None and the configuration cannot be found,
            either because astrometry_net_data is not setup in eups
            or because the setup version does not include the file "andConfig.py"
        """
        pipeBase.Task.__init__(self, config=config, **kwargs)
        self.config = config
        # this is not a subtask because it cannot safely be retargeted
        self.refObjLoader = LoadAstrometryNetObjectsTask(
            config = self.config,
            andConfig = andConfig,
            log = self.log,
            name = "loadAN",
        )
        self.refObjLoader._readIndexFiles()

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

    def _getImageParams(self, exposure=None, bbox=None, wcs=None, filterName=None, wcsRequired=True):
        """Get image parameters

        @param[in] exposure  exposure (an afwImage.Exposure) or None
        @param[in] bbox  bounding box (an afwGeom.Box2I) or None; if None then bbox must be specified
        @param[in] wcs  WCS (an afwImage.Wcs) or None; if None then exposure must be specified
        @param[in] filterName  filter name, a string, or None; if None exposure must be specified
        @param[in] wcsRequired  if True then either wcs must be specified or exposure must contain a wcs;
            if False then the returned wcs may be None
        @return these items:
        - bbox  bounding box; guaranteed to be set
        - wcs  WCS if known, else None
        - filterName filter name if known, else None
        @throw RuntimeError if bbox cannot be determined, or wcs cannot be determined and wcsRequired True
        """
        if exposure is not None:
            if bbox is None:
                bbox = exposure.getBBox()
                self.log.logdebug("Setting bbox = %s from exposure metadata" % (bbox,))
            if wcs is None:
                self.log.logdebug("Setting wcs from exposure metadata")
                wcs = exposure.getWcs()
            if filterName is None:
                filterName = exposure.getFilter().getName()
                self.log.logdebug("Setting filterName = %r from exposure metadata" % (filterName,))
        if bbox is None:
            raise RuntimeError("bbox or exposure must be specified")
        if wcs is None and wcsRequired:
            raise RuntimeError("wcs or exposure (with a WCS) must be specified")
        return bbox, wcs, filterName

    def useKnownWcs(self, sourceCat, wcs=None, exposure=None, filterName=None, bbox=None, calculateSip=None):
        """!Return an InitialAstrometry object, just like determineWcs,
        but assuming the given input WCS is correct.

        This involves searching for reference sources within the WCS
        area, and matching them to the given 'sourceCat'.  If
        'calculateSip' is set, we will try to compute a TAN-SIP
        distortion correction.

        @param[in] sourceCat  list of detected sources in this image.
        @param[in] wcs  your known WCS, or None to get from exposure
        @param[in] exposure  the exposure holding metadata for this image;
            if None then you must specify wcs, filterName and bbox
        @param[in] filterName  string, filter name, eg "i", or None to get from exposure`
        @param[in] bbox  bounding box of image, or None to get from exposure
        @param[in] calculateSip  calculate SIP distortion terms for the WCS? If None
            then use self.config.calculateSip. To disable WCS fitting set calculateSip=False

        @note this function is also called by 'determineWcs' (via 'determineWcs2'), since the steps are all
        the same.
        """
        # return value:
        astrom = InitialAstrometry()

        if calculateSip is None:
            calculateSip = self.config.calculateSip

        bbox, wcs, filterName = self._getImageParams(
            exposure = exposure,
            bbox = bbox,
            wcs = wcs,
            filterName = filterName,
            wcsRequired = True,
        )
        refCat = self.refObjLoader.loadPixelBox(
            bbox = bbox,
            wcs = wcs,
            filterName = filterName,
            calib = None,
        ).refCat
        astrom.refCat = refCat
        catids = [src.getId() for src in refCat]
        uids = set(catids)
        self.log.logdebug('%i reference sources; %i unique IDs' % (len(catids), len(uids)))
        matches = self._getMatchList(sourceCat, refCat, wcs)
        uniq = set([sm.second.getId() for sm in matches])
        if len(matches) != len(uniq):
            self.log.warn(('The list of matched stars contains duplicate reference source IDs ' +
                        '(%i sources, %i unique ids)') % (len(matches), len(uniq)))
        if len(matches) == 0:
            self.log.warn('No matches found between input sources and reference catalogue.')
            return astrom

        self.log.logdebug('%i reference objects match input sources using input WCS' % (len(matches)))
        astrom.tanMatches = matches
        astrom.tanWcs = wcs

        srcids = [s.getId() for s in sourceCat]
        for m in matches:
            assert(m.second.getId() in srcids)
            assert(m.second in sourceCat)

        if calculateSip:
            sipwcs,matches = self._calculateSipTerms(wcs, refCat, sourceCat, matches, bbox=bbox)
            if sipwcs == wcs:
                self.log.logdebug('Failed to find a SIP WCS better than the initial one.')
            else:
                self.log.logdebug('%i reference objects match input sources using SIP WCS' %
                    (len(matches),))
                astrom.sipWcs = sipwcs
                astrom.sipMatches = matches

        wcs = astrom.getWcs()
        # _getMatchList() modifies the source list RA,Dec coordinates.
        # Here, we make them consistent with the WCS we are returning.
        for src in sourceCat:
            src.updateCoord(wcs)
        astrom.matchMeta = _createMetadata(bbox, wcs, filterName)
        return astrom

    def determineWcs(self,
                     sourceCat,
                     exposure,
                     **kwargs):
        """Find a WCS solution for the given 'sourceCat' in the given
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

        'bbox', bounding box of exposure; defaults to exposure.getBBox()

        """
        assert(exposure is not None)

        margs = kwargs.copy()
        if 'searchRadius' not in margs:
            margs.update(searchRadius = self.config.raDecSearchRadius * afwGeom.degrees)
        if 'usePixelScale' not in margs:
            margs.update(usePixelScale = self.config.useWcsPixelScale)
        if 'useRaDecCenter' not in margs:
            margs.update(useRaDecCenter = self.config.useWcsRaDecCenter)
        if 'useParity' not in margs:
            margs.update(useParity = self.config.useWcsParity)
        margs.update(exposure=exposure)
        return self.determineWcs2(sourceCat=sourceCat, **margs)

    def determineWcs2(self, sourceCat, **kwargs):
        """Get a blind astrometric solution for the given catalog of sources.

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
        """
        wcs,qa = self.getBlindWcsSolution(sourceCat, **kwargs)
        kw = {}
        # Keys passed to useKnownWcs
        for key in ['exposure', 'bbox', 'filterName']:
            if key in kwargs:
                kw[key] = kwargs[key]
        astrom = self.useKnownWcs(sourceCat, wcs=wcs, **kw)
        astrom.solveQa = qa
        astrom.tanWcs = wcs
        return astrom

    def getBlindWcsSolution(self, sourceCat,
                            exposure=None,
                            wcs=None,
                            bbox=None,
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

        bbox, wcs, filterName = self._getImageParams(
            exposure = exposure,
            bbox = bbox,
            wcs = wcs,
            filterName = filterName,
            wcsRequired = False,
        )

        bboxD = afwGeom.Box2D(bbox)
        xc, yc = bboxD.getCenter()
        parity = None

        if wcs is not None:
            if pixelScale is None:
                if usePixelScale:
                    pixelScale = wcs.pixelScale()
                    self.log.logdebug('Setting pixel scale estimate = %.3f from given WCS estimate' %
                                (pixelScale.asArcseconds()))

            if radecCenter is None:
                if useRaDecCenter:
                    radecCenter = wcs.pixelToSky(xc, yc)
                    self.log.logdebug(('Setting RA,Dec center estimate = (%.3f, %.3f) from given WCS '
                                 + 'estimate, using pixel center = (%.1f, %.1f)') %
                                (radecCenter.getLongitude().asDegrees(),
                                 radecCenter.getLatitude().asDegrees(), xc, yc))

            if searchRadius is None:
                if useRaDecCenter:
                    assert(pixelScale is not None)
                    pixRadius = math.hypot(*bboxD.getDimensions()) / 2
                    searchRadius = (pixelScale * pixRadius * searchRadiusScale)
                    self.log.logdebug(('Using RA,Dec search radius = %.3f deg, from pixel scale, '
                                 + 'image size, and searchRadiusScale = %g') %
                                (searchRadius, searchRadiusScale))
            if useParity:
                parity = wcs.isFlipped()
                self.log.logdebug('Using parity = %s' % (parity and 'True' or 'False'))

        if doTrim:
            n = len(sourceCat)
            if exposure is not None:
                exposureBBoxD = afwGeom.Box2D(exposure.getMaskedImage().getBBox())
            else:
                exposureBBoxD = bboxD
            sourceCat = self._trimBadPoints(sourceCat, exposureBBoxD)
            self.log.logdebug("Trimming: kept %i of %i sources" % (n, len(sourceCat)))

        wcs,qa = self._solve(
            sourceCat = sourceCat,
            wcs = wcs,
            bbox = bbox,
            pixelScale = pixelScale,
            radecCenter = radecCenter,
            searchRadius = searchRadius,
            parity = parity,
            filterName = filterName,
        )
        if wcs is None:
            raise RuntimeError("Unable to match sources with catalog.")
        self.log.info('Got astrometric solution from Astrometry.net')

        rdc = wcs.pixelToSky(xc, yc)
        self.log.logdebug('New WCS says image center pixel (%.1f, %.1f) -> RA,Dec (%.3f, %.3f)' %
                    (xc, yc, rdc.getLongitude().asDegrees(), rdc.getLatitude().asDegrees()))
        return wcs,qa

    def getSipWcsFromWcs(self, wcs, bbox, ngrid=20, linearizeAtCenter=True):
        """!Get a TAN-SIP WCS, starting from an existing WCS.

        It uses your WCS to compute a fake grid of corresponding "stars" in pixel and sky coords,
        and feeds that to the regular SIP code.

        @param[in] wcs  initial WCS
        @param[in] bbox  bounding box of image
        @param[in] ngrid  number of grid points along x and y for fitting (fit at ngrid^2 points)
        @param[in] linearizeAtCenter  if True, get a linear approximation of the input
          WCS at the image center and use that as the TAN initialization for
          the TAN-SIP solution.  You probably want this if your WCS has its
          CRPIX outside the image bounding box.
        """
        # Ugh, build src and ref tables
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        key = srcSchema.addField("centroid", type="PointD")
        srcTable = afwTable.SourceTable.make(srcSchema)
        srcTable.defineCentroid("centroid")
        srcs = srcTable
        refs = afwTable.SimpleTable.make(afwTable.SimpleTable.makeMinimalSchema())
        cref = []
        csrc = []
        (W,H) = bbox.getDimensions()
        x0, y0 = bbox.getMin()
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
            crpix = afwGeom.Box2D(bbox).getCenter()
            crval = wcs.pixelToSky(crpix)
            crval = crval.getPosition(afwGeom.degrees)
            # Linearize *AT* crval to get effective CD at crval.
            # (we use the default skyUnit of degrees as per WCS standard)
            aff = wcs.linearizePixelToSky(crval)
            cd = aff.getLinear().getMatrix()
            wcs = afwImage.Wcs(crval, crpix, cd)

        return self.getSipWcsFromCorrespondences(wcs, cref, csrc, (W,H),
                                                 x0=x0, y0=y0)


    def getSipWcsFromCorrespondences(self, origWcs, refCat, sourceCat, bbox):
        """Produce a SIP solution given a list of known correspondences.

        Unlike _calculateSipTerms, this does not iterate the solution;
        it assumes you have given it a good sets of corresponding stars.

        NOTE that "refCat" and "sourceCat" are assumed to be the same length;
        entries "refCat[i]" and "sourceCat[i]" are assumed to be correspondences.

        @param[in] origWcs  the WCS to linearize in order to get the TAN part of the TAN-SIP WCS.
        @param[in] refCat  reference source catalog
        @param[in] sourceCat  source catalog
        @param[in] bbox  bounding box of image
        """
        sipOrder = self.config.sipOrder
        matches = []
        for ci,si in zip(refCat, sourceCat):
            matches.append(afwTable.ReferenceMatch(ci, si, 0.))

        sipObject = astromSip.makeCreateWcsWithSip(matches, origWcs, sipOrder, bbox)
        return sipObject.getNewWcs()

    def _calculateSipTerms(self, origWcs, refCat, sourceCat, matches, bbox):
        """!Iteratively calculate SIP distortions and regenerate matches based on improved WCS.

        @param[in] origWcs  original WCS object, probably (but not necessarily) a TAN WCS;
           this is used to set the baseline when determining whether a SIP
           solution is any better; it will be returned if no better SIP solution
           can be found.
        @param[in] refCat  reference source catalog
        @param[in] sourceCat  sources in the image to be solved
        @param[in] matches  list of supposedly matched sources, using the "origWcs".
        @param[in] bbox  bounding box of image, which is used when finding reverse SIP coefficients.
        """
        sipOrder = self.config.sipOrder
        wcs = origWcs

        lastMatchSize = len(matches)
        lastMatchStats = self._computeMatchStatsOnSky(wcs=wcs, matchList=matches)
        for i in range(self.config.maxIter):
            # fit SIP terms
            try:
                sipObject = astromSip.makeCreateWcsWithSip(matches, wcs, sipOrder, bbox)
                proposedWcs = sipObject.getNewWcs()
                self.plotSolution(matches, proposedWcs, bbox.getDimensions())
            except pexExceptions.Exception as e:
                self.log.warn('Failed to calculate distortion terms. Error: ' + str(e))
                break

            # use new WCS to get new matchlist.
            proposedMatchlist = self._getMatchList(sourceCat, refCat, proposedWcs)
            proposedMatchSize = len(proposedMatchlist)
            proposedMatchStats = self._computeMatchStatsOnSky(wcs=proposedWcs, matchList=proposedMatchlist)

            self.log.logdebug(
                "SIP iteration %i: %i objects match, previous = %i;" %
                    (i, proposedMatchSize, lastMatchSize) +
                " clipped mean scatter = %s arcsec, previous = %s; " %
                    (proposedMatchStats.distMean.asArcseconds(), lastMatchStats.distMean.asArcseconds()) +
                " max match dist = %s arcsec, previous = %s" %
                    (proposedMatchStats.maxMatchDist.asArcseconds(),
                        lastMatchStats.maxMatchDist.asArcseconds())
            )

            if lastMatchStats.maxMatchDist <= proposedMatchStats.maxMatchDist:
                self.log.logdebug(
                    "Fit WCS: use iter %s because max match distance no better in next iter: " % (i-1,) +
                    " %g < %g arcsec" % (lastMatchStats.maxMatchDist.asArcseconds(),
                                        proposedMatchStats.maxMatchDist.asArcseconds()))
                break


            wcs = proposedWcs
            matches = proposedMatchlist
            lastMatchSize = proposedMatchSize
            lastMatchStats = proposedMatchStats

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
                    import pdb
                    pdb.set_trace()
                elif reply == "Q":
                    sys.exit(1)
                break

    def _computeMatchStatsOnSky(self, wcs, matchList):
        """Compute on-sky radial distance statistics for a match list

        @param[in] wcs  WCS for match list; an lsst.afw.image.Wcs
        @param[in] matchList  list of matches between reference object and sources;
            a list of lsst.afw.table.ReferenceMatch;
            the source centroid and reference object coord are read

        @return a pipe_base Struct containing these fields:
        - distMean  clipped mean of on-sky radial separation
        - distStdDev  clipped standard deviation of on-sky radial separation
        - maxMatchDist  distMean + self.config.matchDistanceSigma*distStdDev
        """
        distStatsInRadians = makeMatchStatisticsInRadians(wcs, matchList,
            afwMath.MEANCLIP | afwMath.STDEVCLIP)
        distMean = distStatsInRadians.getValue(afwMath.MEANCLIP)*afwGeom.radians
        distStdDev = distStatsInRadians.getValue(afwMath.STDEVCLIP)*afwGeom.radians
        return pipeBase.Struct(
            distMean = distMean,
            distStdDev = distStdDev,
            maxMatchDist = distMean + self.config.matchDistanceSigma*distStdDev,
        )

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
            #self.log.logdebug("source: x,y (%.1f, %.1f), RA,Dec (%.3f, %.3f)" %
            #(src.getX(), src.getY(), src.getRa().asDegrees(), src.getDec().asDegrees()))
            #for src in refCat:
            #self.log.logdebug("ref: RA,Dec (%.3f, %.3f)" %
            #(src.getRa().asDegrees(), src.getDec().asDegrees()))
            self.loginfo('_getMatchList: %i sources, %i reference sources' % (len(sourceCat), len(refCat)))
            if len(sourceCat):
                self.loginfo(
                    'Source range: x [%.1f, %.1f], y [%.1f, %.1f], RA [%.3f, %.3f], Dec [%.3f, %.3f]' %
                    (min(X), max(X), min(Y), max(Y), min(R1), max(R1), min(D1), max(D1)))
            if len(refCat):
                self.loginfo('Reference range: RA [%.3f, %.3f], Dec [%.3f, %.3f]' %
                             (min(R2), max(R2), min(D2), max(D2)))
            raise RuntimeError('No matches found between image and catalogue')
        matches = astromSip.cleanBadPoints.clean(matches, wcs, nsigma=clean)
        return matches

    def getColumnName(self, filterName, columnMap, default=None):
        """
        Returns the column name in the astrometry_net_data index file that will be used
        for the given filter name.

        @param filterName   Name of filter used in exposure
        @param columnMap    Dict that maps filter names to column names
        @param default      Default column name
        """
        filterName = self.config.filterMap.get(filterName, filterName) # Exposure filter --> desired filter
        try:
            return columnMap[filterName] # Desired filter --> a_n_d column name
        except KeyError:
            self.log.warn("No column in configuration for filter '%s'; using default '%s'" %
                          (filterName, default))
            return default

    def _solve(self, sourceCat, wcs, bbox, pixelScale, radecCenter,
               searchRadius, parity, filterName=None):
        solver = self.refObjLoader._getSolver()

        imageSize = bbox.getDimensions()
        x0, y0 = bbox.getMin()

        # select sources with valid x,y, flux
        xybb = afwGeom.Box2D()
        goodsources = afwTable.SourceCatalog(sourceCat.table)
        badkeys = [goodsources.getSchema().find(name).key for name in self.config.badFlags]

        for s in sourceCat:
            if np.isfinite(s.getX()) and np.isfinite(s.getY()) and np.isfinite(s.getPsfFlux()) \
                and self._isGoodSource(s, badkeys):
                goodsources.append(s)
                xybb.include(afwGeom.Point2D(s.getX() - x0, s.getY() - y0))
        self.log.info("Number of selected sources for astrometry : %d" %(len(goodsources)))
        if len(goodsources) < len(sourceCat):
            self.log.logdebug('Keeping %i of %i sources with finite X,Y positions and PSF flux' %
                              (len(goodsources), len(sourceCat)))
        self.log.logdebug(('Feeding sources in range x=[%.1f, %.1f], y=[%.1f, %.1f] ' +
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
            self.log.logdebug(
                'Searching for matches with pixel scale = %g +- %g %% -> range [%g, %g] arcsec/pix' %
                (scale, 100.*(dscale-1.), lo, hi))

        if parity is not None:
            solver.setParity(parity)
            self.log.logdebug('Searching for match with parity = ' + str(parity))

        # Find and load index files within RA,Dec range and scale range.
        if radecCenter is not None:
            multiInds = self.refObjLoader._getMIndexesWithinRange(radecCenter, searchRadius)
        else:
            multiInds = self.refObjLoader.multiInds
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

        import lsstDebug
        if lsstDebug.Info(__name__).display:
            # Use separate context for display, since astrometry.net can segfault if we don't...
            with LoadMultiIndexes(toload_multiInds):
                displayAstrometry(refCat=self.refObjLoader.loadPixelBox(bbox, wcs, filterName).refCat,
                                  frame=lsstDebug.Info(__name__).frame, pause=lsstDebug.Info(__name__).pause)

        with LoadMultiIndexes(toload_multiInds):
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

    @staticmethod
    def _trimBadPoints(sourceCat, bbox, wcs=None):
        """Remove elements from catalog whose xy positions are not within the given bbox.

        sourceCat:  a Catalog of SimpleRecord or SourceRecord objects
        bbox: an afwImage.Box2D
        wcs:  if not None, will be used to compute the xy positions on-the-fly;
              this is required when sources actually contains SimpleRecords.

        Returns:
        a list of Source objects with xAstrom,yAstrom within the bbox.
        """
        keep = type(sourceCat)(sourceCat.table)
        for s in sourceCat:
            point = s.getCentroid() if wcs is None else wcs.skyToPixel(s.getCoord())
            if bbox.contains(point):
                keep.append(s)
        return keep

def _createMetadata(bbox, wcs, filterName):
    """
    Create match metadata entries required for regenerating the catalog

    @param bbox  bounding box of image (pixels)
    @param filterName Name of filter, used for magnitudes
    @return Metadata
    """
    meta = dafBase.PropertyList()

    bboxD = afwGeom.Box2D(bbox)
    cx, cy = bboxD.getCenter()
    radec = wcs.pixelToSky(cx, cy).toIcrs()
    meta.add('RA', radec.getRa().asDegrees(), 'field center in degrees')
    meta.add('DEC', radec.getDec().asDegrees(), 'field center in degrees')
    pixelRadius = math.hypot(*bboxD.getDimensions())/2.0
    skyRadius = wcs.pixelScale() * pixelRadius
    meta.add('RADIUS', skyRadius.asDegrees(),
             'field radius in degrees, approximate')
    meta.add('SMATCHV', 1, 'SourceMatchVector version number')
    if filterName is not None:
        meta.add('FILTER', filterName, 'LSST filter name for tagalong data')
    return meta
