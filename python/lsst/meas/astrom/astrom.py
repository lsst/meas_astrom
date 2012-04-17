import os
import math

import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.algorithms.utils as maUtils

from .config import MeasAstromConfig, AstrometryNetDataConfig
import sip as astromSip

import numpy as np # for isfinite()

# Object returned by determineWcs.
class InitialAstrometry(object):
    '''
    Fields set by determineWcs():

    solveQa (PropertyList)
    tanWcs (Wcs)
    tanMatches (MatchList)
    if sip:
       sipWcs (Wcs)
       sipMatches (MatchList)
    astrom.matchMetadata (PropertyList)
    astrom.wcs (= sipWcs, if available, or tanWcs)
    astrom.matches (= sipMatches if available, else tanMatches)
    
    '''
    def __init__(self):
        self.matches = None
        self.wcs = None
    def getMatches(self):
        return self.matches
    def getWcs(self):
        return self.wcs
    def getMatchMetadata(self):
        return getattr(self, 'matchMetadata', None)

class Astrometry(object):
    ConfigClass = MeasAstromConfig

    def __init__(self,
                 config,
                 andConfig=None,
                 log=None,
                 logLevel=pexLog.Log.INFO):
        '''
        conf: an AstromConfig object.
        andConfig: an AstromNetDataConfig object
        log: a pexLogging.Log
        logLevel: if log is None, the log level to use
        '''
        self.config = config
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
            andConfig.load(fn)

        self.andConfig = andConfig
        self._readIndexFiles()

    def _readIndexFiles(self):
        import astrometry_net as an
        self.inds = []
        for fn in self.andConfig.indexFiles:
            self.log.log(self.log.DEBUG, 'Adding index file %s' % fn)
            fn = self._getIndexPath(fn)
            self.log.log(self.log.DEBUG, 'Path: %s' % fn)
            ind = an.index_load(fn, an.INDEX_ONLY_LOAD_METADATA, None);
            if ind:
                self.inds.append(ind)
                self.log.log(self.log.DEBUG, '  index %i, hp %i (nside %i), nstars %i, nquads %i' %
                             (ind.indexid, ind.healpix, ind.hpnside, ind.nstars, ind.nquads))
            else:
                raise RuntimeError('Failed to read index file: "%s"' % fn)

    def _debug(self, s):
        self.log.log(self.log.DEBUG, s)
    def _warn(self, s):
        self.log.log(self.log.WARN, s)

    def setAndConfig(self, andconfig):
        self.andConfig = andconfig

    def determineWcs(self,
                     sources,
                     exposure,
                     **kwargs):
        '''
        Version of determineWcs(), meant for pipeline use, that gets
        almost all its parameters from config or reasonable defaults.
        '''
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

        return self.determineWcs2(sources, exposure, **margs)
        

    def determineWcs2(self,
                      sources,
                      exposure=None,
                      wcs=None,
                      imageSize=None,
                      radecCenter=None,
                      searchRadius=None,
                      pixelScale=None,
                      filterName=None,
                      doTrim=False,
                      usePixelScale=True,
                      useRaDecCenter=True,
                      useParity=True,
                      searchRadiusScale=2.):
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
        pixelScale: afwGeom::Angle per pixel.
        radecCenter: afwCoord::Coord
        '''

        if not useRaDecCenter and radecCenter is not None:
            raise RuntimeError('radecCenter is set, but useRaDecCenter is False.  Make up your mind!')
        if not usePixelScale and pixelScale is not None:
            raise RuntimeError('pixelScale is set, but usePixelScale is False.  Make up your mind!')
        
        # return value:
        astrom = InitialAstrometry()
        
        if exposure is not None:
            if filterName is None:
                filterName = exposure.getFilter().getName()
            if imageSize is None:
                imageSize = (exposure.getWidth(), exposure.getHeight())
            if wcs is None:
                wcs = exposure.getWcs()

        if imageSize is None:
            # Could guess from the extent of the Sources...
            raise RuntimeError('Image size must be specified by passing "exposure" or "imageSize"')
        W,H = imageSize
        xc, yc = W/2. + 0.5, H/2. + 0.5

        parity = None
        
        if wcs is not None:
            if pixelScale is None:
                if usePixelScale:
                    pixelScale = wcs.pixelScale()

            if radecCenter is None:
                if useRaDecCenter:
                    radecCenter = wcs.pixelToSky(xc, yc)

            if searchRadius is None:
                if useRaDecCenter:
                    assert(pixelScale is not None)
                    searchRadius = (pixelScale * math.hypot(W,H)/2. *
                                    searchRadiusScale)
                
            if useParity:
                parity = wcs.isFlipped()

        if doTrim:
            n = len(sources)
            if exposure is not None:
                bbox = afwGeom.Box2D(exposure.getMaskedImage().getBBox(afwImage.PARENT))
            else:
                # CHECK -- half-pixel issues here?
                bbox = afwGeom.Box2D(afwGeom.Point2D(0.,0.), afwGeom.Point2D(W, H))
            sources = _trimBadPoints(sources, bbox)
            self._debug("Trimming: kept %i of %i sources" % (n, len(sources)))

        wcs,qa = self._solve(sources, wcs, imageSize, pixelScale, radecCenter, searchRadius, parity,
                             filterName)
        if wcs is None:
            raise RuntimeError("Unable to match sources with catalog.")

        pixelMargin = 50.
        cat = self.getReferenceSourcesForWcs(wcs, imageSize, filterName, pixelMargin)

        catids = [src.getId() for src in cat]
        uids = set(catids)
        self.log.logdebug('%i reference sources; %i unique IDs' % (len(catids), len(uids)))

        matchList = self._getMatchList(sources, cat, wcs)

        uniq = set([sm.second.getId() for sm in matchList])
        if len(matchList) != len(uniq):
            self._warn('The list of matched stars contains duplicate reference source IDs (%i sources, %i unique ids)'
                       % (len(matchList), len(uniq)))
        if len(matchList) == 0:
            self._warn('No matches found between input sources and reference catalogue.')
            return astrom

        self._debug('%i reference objects match input sources using linear WCS' % (len(matchList)))

        astrom.solveQa = qa
        astrom.tanWcs = wcs
        astrom.tanMatches = matchList

        srcids = [s.getId() for s in sources]
        for m in matchList:
            assert(m.second.getId() in srcids)
            assert(m.second in sources)

        if self.config.calculateSip:
            wcs,matchList = self._calculateSipTerms(wcs, cat, sources, matchList)
            astrom.sipWcs = wcs
            astrom.sipMatches = matchList

        meta = _createMetadata(W, H, wcs, filterName)
        #matchListMeta = solver.getMatchedIndexMetadata()
        #moreMeta.combine(matchListMeta)

        astrom.matchMetadata = meta
        astrom.wcs = wcs

        astrom.matches = afwTable.ReferenceMatchVector()
        for m in matchList:
            astrom.matches.push_back(m)

        return astrom

    #### FIXME!
    def _calculateSipTerms(self, origWcs, cat, sources, matchList):
        '''Iteratively calculate sip distortions and regenerate matchList based on improved wcs'''
        sipOrder = self.config.sipOrder
        wcs = origWcs

        i=0
        while True:
            try:
                sipObject = astromSip.CreateWcsWithSip(matchList, wcs, sipOrder)
                proposedWcs = sipObject.getNewWcs()
            except LsstCppException, e:
                self._warn('Failed to calculate distortion terms. Error: ' + str(e))
                break

            matchSize = len(matchList)
            self._debug('Sip Iteration %i: %i objects match. rms scatter is %g arcsec or %g pixels' %
                        (i, matchSize, sipObject.getScatterOnSky().asArcseconds(), sipObject.getScatterInPixels()))
            # use new WCS to get new matchlist.
            proposedMatchlist = self._getMatchList(sources, cat, proposedWcs)
            if len(proposedMatchlist) <= matchSize:
                # We're regressing, so stop
                break
            wcs = proposedWcs
            matchList = proposedMatchlist
            matchSize = len(matchList)
            i += 1

        return wcs, matchList

    def _getMatchList(self, sources, cat, wcs):
        dist = self.config.catalogMatchDist * afwGeom.arcseconds
        clean = self.config.cleaningParameter
        matcher = astromSip.MatchSrcToCatalogue(cat, sources, wcs, dist)
        matchList = matcher.getMatches()
        if matchList is None:
            raise RuntimeError('No matches found between image and catalogue')
        matchList = astromSip.cleanBadPoints.clean(matchList, wcs, nsigma=clean)
        return matchList

    def getCatalogFilterName(self, filterName):
        '''
        Returns the column name in the astrometry_net_data index file that will be used
        for the given filter name.
        '''
        return self._mapFilterName(filterName, self.andConfig.defaultMagColumn)

    def _mapFilterName(self, filterName, default=None):
        filterName = self.config.filterMap.get(filterName, filterName) # Exposure filter --> desired filter
        try:
            return self.andConfig.magColumnMap[filterName] # Desired filter --> a_n_d column name
        except KeyError:
            self.log.warn("No mag column in configuration for filter '%s'; using default '%s'" %
                          (filterName, default))
            return default


    def getReferenceSourcesForWcs(self, wcs, imageSize, filterName, pixelMargin,
                                  trim=True):
        W,H = imageSize
        xc, yc = W/2. + 0.5, H/2. + 0.5
        rdc = wcs.pixelToSky(xc, yc)
        ra,dec = rdc.getLongitude(), rdc.getLatitude()
        pixelScale = wcs.pixelScale()
        rad = pixelScale * (math.hypot(W,H)/2. + pixelMargin)
        cat = self.getReferenceSources(ra, dec, rad, filterName)
        # NOTE: reference objects don't have (x,y) anymore, so we can't apply WCS to set x,y positions
        if trim:
            # cut to image bounds + margin.
            bbox = afwGeom.Box2D(afwGeom.Point2D(0.,0.), afwGeom.Point2D(W, H))
            bbox.grow(pixelMargin)
            cat = self._trimBadPoints(cat, bbox, wcs=wcs) # passing wcs says to compute x,y on-the-fly
        return cat


    def getReferenceSources(self, ra, dec, radius, filterName, allFluxes=False):
        '''
        Searches for reference-catalog sources (in the
        astrometry_net_data files) in the requested RA,Dec region
        (afwGeom::Angle objects), with the requested radius (also an
        Angle).  The flux values will be set based on the requested
        filter (None => default filter).
        
        Returns: an lsst.afw.table.SimpleCatalog of reference objects
        '''
        solver = self._getSolver()

        sgCol = self.andConfig.starGalaxyColumn
        varCol = self.andConfig.variableColumn
        idcolumn = self.andConfig.idColumn

        magCol = self.getCatalogFilterName(filterName)
        magerrCol = self.andConfig.magErrorColumnMap.get(filterName, None)

        if allFluxes:
            magCol = [magCol] + [x for x in self.andConfig.magColumnMap.keys() if x != magCol]
            tmp = [self.andConfig.magErrorColumnMap.get(x, None) for x in \
                       self.andConfig.magErrorColumnMap.keys()]
            magerrCol = [magerrCol] + [x for x in tmp if x != magerrCol]

        '''
        Note about multiple astrometry_net index files and duplicate IDs:

        -as of astrometry_net 0.30, we take a reference catalog and build
         a set of astrometry_net index files from it, with each one covering a
         region of sky and a range of angular scales.  The index files covering
         the same region of sky at different scales use exactly the same stars.
         Therefore, if we search every index file, we will get multiple copies of
         each reference star (one from each index file).
         For now, we have the "unique_ids" option to solver.getCatalog().
         -recall that the index files to be used are specified in the
          AstrometryNetDataConfig.indexFiles flat list.

        -as of astrometry_net 0.40, we have the capability to share
         the reference stars between index files (called
         "multiindex"), so we will no longer have to repeat the
         reference stars in each index.  We will, however, have to
         change the way the index files are configured to take
         advantage of this functionality.  Once this is in place, we
         can eliminate the horrid ID checking and deduplication (in solver.getCatalog()).
         -multiindex files will be specified in the
          AstrometryNetDatConfig.multiIndexFiles list-of-lists; first
          element is the filename containing the stars, subsequent
          elements are filenames containing the index structures.
          We may be able to backwards-compatibly build this from the flat indexFiles
          list if we assume things about the filenames.
        '''
        cat = solver.getCatalog(self.inds, ra.asDegrees(), dec.asDegrees(), radius.asDegrees(),
                                idcolumn, magCol, magerrCol, sgCol, varCol)
        del solver
        return cat

    def _solve(self, sources, wcs, imageSize, pixelScale, radecCenter,
               searchRadius, parity, filterName=None):
        solver = self._getSolver()

        # select sources with valid x,y, flux
        goodsources = afwTable.SourceCatalog(sources.table)
        for s in sources:
            if np.isfinite(s.getX()) and np.isfinite(s.getY()) and np.isfinite(s.getPsfFlux()):
                goodsources.append(s)
        if len(goodsources) < len(sources):
            self.log.logdebug('Keeping %i of %i sources with finite X,Y positions and PSF flux' %
                              (len(goodsources), len(sources)))
        # setStars sorts them by PSF flux.
        solver.setStars(goodsources)
        solver.setMaxStars(self.config.maxStars)
        solver.setImageSize(*imageSize)
        solver.setMatchThreshold(self.config.matchThreshold)
        if radecCenter is not None:
            ra = radecCenter.getLongitude().asDegrees()
            dec = radecCenter.getLatitude().asDegrees()
            solver.setRaDecRadius(ra, dec, searchRadius.asDegrees())
            self.log.logdebug('Searching for match around RA,Dec = (%g, %g) with radius %g deg' %
                              (ra, dec, searchRadius.asDegrees()))

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

        solver.addIndices(self.inds)
        active = solver.getActiveIndexFiles()
        self.log.logdebug('Searching for match in %i of %i index files: [ ' % (len(active), len(self.inds)) +
                          ', '.join(ind.indexname for ind in active) + ' ]')

        cpulimit = self.config.maxCpuTime

        solver.run(cpulimit)
        if solver.didSolve():
            self.log.logdebug('Solved!')
            wcs = solver.getWcs()
            self.log.logdebug('WCS: %s' % wcs.getFitsMetadata().toString())
            
        else:
            self.log.warn('Did not get an astrometric solution from Astrometry.net')
            wcs = None
            # Gather debugging info...

            # -are there any reference stars in the proposed search area?
            if radecCenter is not None:
                ra = radecCenter.getLongitude()
                dec = radecCenter.getLatitude()
                refs = self.getReferenceSources(ra, dec, searchRadius, filterName)
                self.log.info('Searching around RA,Dec = (%g,%g) with radius %g deg yields %i reference-catalog sources' %
                              (ra.asDegrees(), dec.asDegrees(), searchRadius.asDegrees(), len(refs)))

        qa = solver.getSolveStats()
        self.log.logdebug('qa: %s' % qa.toString())
        return wcs, qa

    def _getIndexPath(self, fn):
        if os.path.isabs(fn):
            return fn
        andir = os.getenv('ASTROMETRY_NET_DATA_DIR')
        if andir is not None:
            fn2 = os.path.join(andir, fn)
            if os.path.exists(fn2):
                return fn2
        fn2 = os.path.abspath(fn)
        return fn2
                    

    def _getSolver(self):
        import astrometry_net as an
        solver = an.solver_new()
        # HACK, set huge default pixel scale range.
        lo,hi = 0.01, 3600.
        solver.setPixelScaleRange(lo, hi)
        return solver

    @staticmethod
    def _trimBadPoints(sources, bbox, wcs=None):
        '''Remove elements from catalog whose xy positions are not within the given bbox.

        sources:  a Catalog of SimpleRecord or SourceRecord objects
        bbox: an afwImage.Box2D
        wcs:  if not None, will be used to compute the xy positions on-the-fly;
              this is required when sources actually contains SimpleRecords.
        
        Returns:
        a list of Source objects with xAstrom,yAstrom within the bbox.
        '''
        keep = type(sources)(sources.table)
        for s in sources:
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
        -- sets the "matchList[*].first" and "matchList[*].second" entries
        to point to the sources in the "sources" argument, and to the
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
        # all in deg.
        ra = matchmeta.getDouble('RA') * afwGeom.degrees
        dec = matchmeta.getDouble('DEC') * afwGeom.degrees
        rad = matchmeta.getDouble('RADIUS') * afwGeom.degrees
        self.log.logdebug('Searching RA,Dec %.3f,%.3f, radius %.1f arcsec, filter "%s"' %
                          (ra.asDegrees(), dec.asDegrees(), rad.asArcseconds(), filterName))
        refCat = self.getReferenceSources(ra, dec, rad, filterName)
        self.log.logdebug('Found %i reference catalog sources in range' % len(refCat))
        refCat.sort()
        sourceCat.sort()
        return afwTable.unpackMatches(packedMatches, refCat, sourceCat)


def _createMetadata(width, height, wcs, filterName):
    """
    Create match metadata entries required for regenerating the catalog

    @param width Width of the image (pixels)
    @param height Height of the image (pixels)
    @param filterName Name of filter, used for magnitudes
    @return Metadata
    """
    meta = dafBase.PropertyList()

    #andata = os.environ.get('ASTROMETRY_NET_DATA_DIR')
    #if andata is None:
    #    meta.add('ANEUPS', 'none', 'ASTROMETRY_NET_DATA_DIR')
    #else:
    #    andata = os.path.basename(andata)
    #    meta.add('ANEUPS', andata, 'ASTROMETRY_NET_DATA_DIR')

    # cache: field center and size.  These may be off by 1/2 or 1 or 3/2 pixels.
    cx,cy = 0.5 + width/2., 0.5 + height/2.
    radec = wcs.pixelToSky(cx, cy).toIcrs()
    meta.add('RA', radec.getRa().asDegrees(), 'field center in degrees')
    meta.add('DEC', radec.getDec().asDegrees(), 'field center in degrees')
    imgSize = wcs.pixelScale() * math.hypot(width, height)/2.
    meta.add('RADIUS', imgSize.asDegrees(),
             'field radius in degrees, approximate')
    meta.add('SMATCHV', 1, 'SourceMatchVector version number')
    if filterName is not None:
        meta.add('FILTER', filterName, 'LSST filter name for tagalong data')
    #meta.add('STARGAL', stargalName, 'star/galaxy name for tagalong data')
    #meta.add('VARIABLE', variableName, 'variability name for tagalong data')
    #meta.add('MAGERR', magerrName, 'magnitude error name for tagalong data')
    return meta

def readMatches(butler, dataId, sourcesName='icSrc', matchesName='icMatch'):
    """Read matches, sources and catalogue; combine.

    @param butler Data butler
    @param dataId Data identifier for butler
    @param sourcesName Name for sources from butler
    @param matchesName Name for matches from butler
    @returns Matches
    """
    sources = butler.get(sourcesName, dataId)
    packedMatches = butler.get(matchesName, dataId)
    
    astrom = Astrometry(MeasAstromConfig())
    return astrom.joinMatchListWithCatalog(packedMatches, sources)

