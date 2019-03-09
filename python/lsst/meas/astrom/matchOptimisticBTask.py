
__all__ = ["MatchOptimisticBTask", "MatchOptimisticBConfig",
           "MatchTolerance"]

import math

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .setMatchDistance import setMatchDistance
from .matchOptimisticB import matchOptimisticB, MatchOptimisticBControl


class MatchTolerance:
    """Stores match tolerances for use in `lsst.meas.astrom.AstrometryTask` and
    later iterations of the matcher.

    MatchOptimsiticBTask relies on a maximum distance for matching
    set by either the default in MatchOptimisticBConfig or the 2 sigma
    scatter found after AstrometryTask has fit for a wcs.

    Parameters
    ----------
    maxMatchDist : `lsst.geom.Angle`
        Current maximum distance to consider a match.
    """

    def __init__(self, maxMatchDist=None):
        self.maxMatchDist = maxMatchDist


class MatchOptimisticBConfig(pexConfig.Config):
    """Configuration for MatchOptimisticBTask
    """
    maxMatchDistArcSec = pexConfig.RangeField(
        doc="Maximum separation between reference objects and sources "
        "beyond which they will not be considered a match (arcsec)",
        dtype=float,
        default=3,
        min=0,
    )
    numBrightStars = pexConfig.RangeField(
        doc="Number of bright stars to use",
        dtype=int,
        default=50,
        min=2,
    )
    minMatchedPairs = pexConfig.RangeField(
        doc="Minimum number of matched pairs; see also minFracMatchedPairs",
        dtype=int,
        default=30,
        min=2,
    )
    minFracMatchedPairs = pexConfig.RangeField(
        doc="Minimum number of matched pairs as a fraction of the smaller of "
        "the number of reference stars or the number of good sources; "
        "the actual minimum is the smaller of this value or minMatchedPairs",
        dtype=float,
        default=0.3,
        min=0,
        max=1,
    )
    maxOffsetPix = pexConfig.RangeField(
        doc="Maximum allowed shift of WCS, due to matching (pixel). "
            "When changing this value, the LoadReferenceObjectsConfig.pixelMargin should also be updated.",
        dtype=int,
        default=300,
        max=4000,
    )
    maxRotationDeg = pexConfig.RangeField(
        doc="Rotation angle allowed between sources and position reference objects (degrees)",
        dtype=float,
        default=1.0,
        max=6.0,
    )
    allowedNonperpDeg = pexConfig.RangeField(
        doc="Allowed non-perpendicularity of x and y (degree)",
        dtype=float,
        default=3.0,
        max=45.0,
    )
    numPointsForShape = pexConfig.Field(
        doc="number of points to define a shape for matching",
        dtype=int,
        default=6,
    )
    maxDeterminant = pexConfig.Field(
        doc="maximum determinant of linear transformation matrix for a usable solution",
        dtype=float,
        default=0.02,
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources for cross-matching",
        default="matcher"
    )

    def setDefaults(self):
        sourceSelector = self.sourceSelector["matcher"]
        sourceSelector.setDefaults()


# The following block adds links to this task from the Task Documentation page.
# \addtogroup LSST_task_documentation
# \{
# \page measAstrom_matchOptimisticBTask
# \ref MatchOptimisticBTask "MatchOptimisticBTask"
# Match sources to reference objects
# \}


class MatchOptimisticBTask(pipeBase.Task):
    """Match sources to reference objects using the Optimistic Pattern Matcher
    B algorithm of Tabur 2007.
    """
    ConfigClass = MatchOptimisticBConfig
    _DefaultName = "matchObjectsToSources"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("sourceSelector")

    def filterStars(self, refCat):
        """Extra filtering pass; subclass if desired.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            Catalog of reference objects.

        Returns
        -------
        trimmedRefCat : `lsst.afw.table.SimpleCatalog`
            Reference catalog with some filtering applied. Currently no
            filtering is applied.
        """
        return refCat

    @pipeBase.timeMethod
    def matchObjectsToSources(self, refCat, sourceCat, wcs, refFluxField,
                              match_tolerance=None):
        """Match sources to position reference stars.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            Reference catalog to match.
        sourceCat : `lsst.afw.table.SourceCatalog`
            Source catalog to match.
        wcs : `lsst.afw.geom.SkyWcs`
            Current WCS of the  exposure containing the sources.
        refFluxField : `str`
            Name of the reference catalog filter to use.
        match_tolerance : `lsst.meas.astrom.MatchTolerance`
            Object containing information from previous
            `lsst.meas.astrom.AstrometryTask` match/fit cycles for use in
            matching. If `None` is config defaults.

        Returns
        -------
        matchResult : `lsst.pipe.base.Struct`
            Result struct with components

            - ``matches`` : List of matches with distance below the maximum match
              distance (`list` of `lsst.afw.table.ReferenceMatch`).
            - ``useableSourceCat`` : Catalog of sources matched and suited for
              WCS fitting (`lsst.afw.table.SourceCatalog`).
            - ``match_tolerance`` : MatchTolerance object updated from this
              match iteration (`lsst.meas.astrom.MatchTolerance`).
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        preNumObj = len(refCat)
        refCat = self.filterStars(refCat)
        numRefObj = len(refCat)

        if self.log:
            self.log.info("filterStars purged %d reference stars, leaving %d stars" %
                          (preNumObj - numRefObj, numRefObj))

        if match_tolerance is None:
            match_tolerance = MatchTolerance()

        # usableSourceCat: sources that are good but may be saturated
        numSources = len(sourceCat)
        selectedSources = self.sourceSelector.run(sourceCat)
        usableSourceCat = selectedSources.sourceCat
        numUsableSources = len(usableSourceCat)
        self.log.info("Purged %d unusable sources, leaving %d usable sources" %
                      (numSources - numUsableSources, numUsableSources))

        if len(usableSourceCat) == 0:
            raise pipeBase.TaskError("No sources are usable")

        del sourceCat  # avoid accidentally using sourceCat; use usableSourceCat or goodSourceCat from now on

        minMatchedPairs = min(self.config.minMatchedPairs,
                              int(self.config.minFracMatchedPairs * min([len(refCat), len(usableSourceCat)])))

        # match usable (possibly saturated) sources and then purge saturated sources from the match list
        usableMatches = self._doMatch(
            refCat=refCat,
            sourceCat=usableSourceCat,
            wcs=wcs,
            refFluxField=refFluxField,
            numUsableSources=numUsableSources,
            minMatchedPairs=minMatchedPairs,
            maxMatchDist=match_tolerance.maxMatchDist,
            sourceFluxField=self.sourceSelector.fluxField,
            verbose=debug.verbose,
        )

        # cull non-good sources
        matches = []
        self._getIsGoodKeys(usableSourceCat.schema)
        for match in usableMatches:
            if self._isGoodTest(match.second):
                # Append the isGood match.
                matches.append(match)

        self.log.debug("Found %d usable matches, of which %d had good sources",
                       len(usableMatches), len(matches))

        if len(matches) == 0:
            raise RuntimeError("Unable to match sources")

        self.log.info("Matched %d sources" % len(matches))
        if len(matches) < minMatchedPairs:
            self.log.warn("Number of matches is smaller than request")

        return pipeBase.Struct(
            matches=matches,
            usableSourceCat=usableSourceCat,
            match_tolerance=match_tolerance,
        )

    def _getIsGoodKeys(self, schema):
        """Retrieve the keys needed for the isGoodTest from the source catalog
        schema.

        Parameters
        ----------
        schema : `lsst.afw.table.Schema`
            Source schema to retrieve `lsst.afw.table.Key` s from.
        """
        self.edgeKey = schema["base_PixelFlags_flag_edge"].asKey()
        self.interpolatedCenterKey = schema["base_PixelFlags_flag_interpolatedCenter"].asKey()
        self.saturatedKey = schema["base_PixelFlags_flag_saturated"].asKey()

    def _isGoodTest(self, source):
        """Test that an object is good for use in the WCS fitter.

        This is a hard coded version of the isGood flag from the old SourceInfo
        class that used to be part of this class.

        Parameters
        ----------
        source : `lsst.afw.table.SourceRecord`
            Source to test.

        Returns
        -------
        isGood : `bool`
            Source passes CCD edge and saturated tests.
        """
        return (not source.get(self.edgeKey) and
                not source.get(self.interpolatedCenterKey) and
                not source.get(self.saturatedKey))

    @pipeBase.timeMethod
    def _doMatch(self, refCat, sourceCat, wcs, refFluxField, numUsableSources, minMatchedPairs,
                 maxMatchDist, sourceFluxField, verbose):
        """Implementation of matching sources to position reference stars.

        Unlike matchObjectsToSources, this method does not check if the sources
        are suitable.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            Catalog of reference objects.
        sourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of detected sources.
        wcs : `lsst.afw.geom.SkyWcs`
            Current best WCS of the image.
        refFluxFioeld : `str`
            Name of flux field in refCat to use.
        numUsableSources : `int`
            Total number of source usable for matching.
        mintMatchPairs : `int`
            Minimum number of objects to match between the refCat and sourceCat
            to consider a valid match.
        maxMatchDist : `lsst.geom.Angle`
            Maximum separation to considering a reference and a source a match.
        sourceFluxField : `str`
            Name of source catalog flux field.
        verbose : `bool`
            Print diagnostic information std::cout

        Returns
        -------
        matches : `list` of `lsst.afw.table.ReferenceMatch`
        """
        numSources = len(sourceCat)
        posRefBegInd = numUsableSources - numSources
        if maxMatchDist is None:
            maxMatchDistArcSec = self.config.maxMatchDistArcSec
        else:
            maxMatchDistArcSec = min(maxMatchDist.asArcseconds(), self.config.maxMatchDistArcSec)
        configMatchDistPix = maxMatchDistArcSec/wcs.getPixelScale().asArcseconds()

        matchControl = MatchOptimisticBControl()
        matchControl.refFluxField = refFluxField
        matchControl.sourceFluxField = sourceFluxField
        matchControl.numBrightStars = self.config.numBrightStars
        matchControl.minMatchedPairs = self.config.minMatchedPairs
        matchControl.maxOffsetPix = self.config.maxOffsetPix
        matchControl.numPointsForShape = self.config.numPointsForShape
        matchControl.maxDeterminant = self.config.maxDeterminant

        for maxRotInd in range(4):
            matchControl.maxRotationDeg = self.config.maxRotationDeg * math.pow(2.0, 0.5*maxRotInd)
            for matchRadInd in range(3):
                matchControl.matchingAllowancePix = configMatchDistPix * math.pow(1.25, matchRadInd)

                for angleDiffInd in range(3):
                    matchControl.allowedNonperpDeg = self.config.allowedNonperpDeg*(angleDiffInd+1)
                    matches = matchOptimisticB(
                        refCat,
                        sourceCat,
                        matchControl,
                        wcs,
                        posRefBegInd,
                        verbose,
                    )
                    if matches is not None and len(matches) > 0:
                        setMatchDistance(matches)
                        return matches
        return matches
