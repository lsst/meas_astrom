
import numpy as np
from scipy.spatial import cKDTree

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwgeom
import lsst.afw.table as afwTable

from .matchOptimisticBTask import MatchTolerance

from .pessimistic_pattern_matcher_b_3D import PessimisticPatternMatcherB

__all__ = ["MatchPessimisticBTask", "MatchPessimisticBConfig",
           "MatchTolerancePessimistic"]


class MatchTolerancePessimistic(MatchTolerance):
    """Stores match tolerances for use in AstrometryTask and later
    iterations of the matcher.

    MatchPessimisticBTask relies on several state variables to be
    preserved over different iterations in the
    AstrometryTask.matchAndFitWcs loop of AstrometryTask.

    Parameters
    ----------
    maxMatchDist : `lsst.geom.Angle`
        Maximum distance to consider a match from the previous match/fit
        iteration.
    autoMaxMatchDist : `lsst.geom.Angle`
        Automated estimation of the maxMatchDist from the sky statistics of the
        source and reference catalogs.
    maxShift : `lsst.geom.Angle`
        Maximum shift found in the previous match/fit cycle.
    lastMatchedPattern : `int`
        Index of the last source pattern that was matched into the reference
        data.
    failedPatternList : `list` of `int`
        Previous matches were found to be false positives.
    PPMbObj : `lsst.meas.astrom.PessimisticPatternMatcherB`
        Initialized Pessimistic pattern matcher object. Storing this prevents
        the need for recalculation of the searchable distances in the PPMB.
    """

    def __init__(self, maxMatchDist=None, autoMaxMatchDist=None,
                 maxShift=None, lastMatchedPattern=None,
                 failedPatternList=None, PPMbObj=None):
        self.maxMatchDist = maxMatchDist
        self.autoMaxMatchDist = autoMaxMatchDist
        self.maxShift = maxShift
        self.lastMatchedPattern = lastMatchedPattern
        self.PPMbObj = PPMbObj
        if failedPatternList is None:
            self.failedPatternList = []
        else:
            self.failedPatternList = failedPatternList


class MatchPessimisticBConfig(pexConfig.Config):
    """Configuration for MatchPessimisticBTask
    """
    numBrightStars = pexConfig.RangeField(
        doc="Number of bright stars to use. Sets the max number of patterns "
            "that can be tested.",
        dtype=int,
        default=200,
        min=2,
    )
    minMatchedPairs = pexConfig.RangeField(
        doc="Minimum number of matched pairs; see also minFracMatchedPairs.",
        dtype=int,
        default=30,
        min=2,
    )
    minFracMatchedPairs = pexConfig.RangeField(
        doc="Minimum number of matched pairs as a fraction of the smaller of "
            "the number of reference stars or the number of good sources; "
            "the actual minimum is the smaller of this value or "
            "minMatchedPairs.",
        dtype=float,
        default=0.3,
        min=0,
        max=1,
    )
    matcherIterations = pexConfig.RangeField(
        doc="Number of softening iterations in matcher.",
        dtype=int,
        default=5,
        min=1,
    )
    maxOffsetPix = pexConfig.RangeField(
        doc="Maximum allowed shift of WCS, due to matching (pixel). "
            "When changing this value, the "
            "LoadReferenceObjectsConfig.pixelMargin should also be updated.",
        dtype=int,
        default=300,
        max=4000,
    )
    maxRotationDeg = pexConfig.RangeField(
        doc="Rotation angle allowed between sources and position reference "
            "objects (degrees).",
        dtype=float,
        default=1.0,
        max=6.0,
    )
    numPointsForShape = pexConfig.Field(
        doc="Number of points to define a shape for matching.",
        dtype=int,
        default=6,
    )
    numPointsForShapeAttempt = pexConfig.Field(
        doc="Number of points to try for creating a shape. This value should "
            "be greater than or equal to numPointsForShape. Besides "
            "loosening the signal to noise cut in the 'matcher' SourceSelector, "
            "increasing this number will solve CCDs where no match was found.",
        dtype=int,
        default=6,
    )
    minMatchDistPixels = pexConfig.RangeField(
        doc="Distance in units of pixels to always consider a source-"
            "reference pair a match. This prevents the astrometric fitter "
            "from over-fitting and removing stars that should be matched and "
            "allows for inclusion of new matches as the wcs improves.",
        dtype=float,
        default=1.0,
        min=0.0,
        max=6.0,
    )
    numPatternConsensus = pexConfig.Field(
        doc="Number of implied shift/rotations from patterns that must agree "
            "before it a given shift/rotation is accepted. This is only used "
            "after the first softening iteration fails and if both the "
            "number of reference and source objects is greater than "
            "numBrightStars.",
        dtype=int,
        default=3,
    )
    numRefRequireConsensus = pexConfig.Field(
        doc="If the available reference objects exceeds this number, "
            "consensus/pessimistic mode will enforced regardless of the "
            "number of available sources. Below this optimistic mode ("
            "exit at first match rather than requiring numPatternConsensus to "
            "be matched) can be used. If more sources are required to match, "
            "decrease the signal to noise cut in the sourceSelector.",
        dtype=int,
        default=1000,
    )
    maxRefObjects = pexConfig.RangeField(
        doc="Maximum number of reference objects to use for the matcher. The "
            "absolute maximum allowed for is 2 ** 16 for memory reasons.",
        dtype=int,
        default=2**16,
        min=0,
        max=2**16 + 1,
    )

    def validate(self):
        pexConfig.Config.validate(self)
        if self.numPointsForShapeAttempt < self.numPointsForShape:
            raise ValueError("numPointsForShapeAttempt must be greater than "
                             "or equal to numPointsForShape.")
        if self.numPointsForShape > self.numBrightStars:
            raise ValueError("numBrightStars must be greater than "
                             "numPointsForShape.")


# The following block adds links to this task from the Task Documentation page.
# \addtogroup LSST_task_documentation
# \{
# \page measAstrom_MatchPessimisticBTask
# \ref MatchPessimisticBTask "MatchPessimisticBTask"
# Match sources to reference objects
# \}


class MatchPessimisticBTask(pipeBase.Task):
    """Match sources to reference objects.
    """

    ConfigClass = MatchPessimisticBConfig
    _DefaultName = "matchObjectsToSources"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)

    @pipeBase.timeMethod
    def matchObjectsToSources(self, refCat, sourceCat, wcs, sourceFluxField, refFluxField,
                              match_tolerance=None):
        """Match sources to position reference stars

        refCat : `lsst.afw.table.SimpleCatalog`
            catalog of reference objects that overlap the exposure; reads
            fields for:

            - coord
            - the specified flux field

        sourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of sources found on an exposure.  This should already be
            down-selected to "good"/"usable" sources in the calling Task.
        wcs : `lsst.afw.geom.SkyWcs`
            estimated WCS
        sourceFluxField: `str`
            field of sourceCat to use for flux
        refFluxField : `str`
            field of refCat to use for flux
        match_tolerance : `lsst.meas.astrom.MatchTolerancePessimistic`
            is a MatchTolerance class object or `None`. This this class is used
            to communicate state between AstrometryTask and MatcherTask.
            AstrometryTask will also set the MatchTolerance class variable
            maxMatchDist based on the scatter AstrometryTask has found after
            fitting for the wcs.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``matches`` : source to reference matches found (`list` of
              `lsst.afw.table.ReferenceMatch`)
            - ``usableSourceCat`` : a catalog of sources potentially usable for
              matching and WCS fitting (`lsst.afw.table.SourceCatalog`).
            - ``match_tolerance`` : a MatchTolerance object containing the
              resulting state variables from the match
              (`lsst.meas.astrom.MatchTolerancePessimistic`).
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        # If we get an empty tolerance struct create the variables we need for
        # this matcher.
        if match_tolerance is None:
            match_tolerance = MatchTolerancePessimistic()

        # Make a name alias here for consistency with older code, and to make
        # it clear that this is a good/usable (cleaned) source catalog.
        goodSourceCat = sourceCat

        numUsableSources = len(goodSourceCat)

        if len(goodSourceCat) == 0:
            raise pipeBase.TaskError("No sources are good")

        minMatchedPairs = min(self.config.minMatchedPairs,
                              int(self.config.minFracMatchedPairs *
                                  min([len(refCat), len(goodSourceCat)])))

        if len(refCat) > self.config.maxRefObjects:
            self.log.warn(
                "WARNING: Reference catalog larger that maximum allowed. "
                "Trimming to %i" % self.config.maxRefObjects)
            trimmedRefCat = self._filterRefCat(refCat, refFluxField)
        else:
            trimmedRefCat = refCat

        doMatchReturn = self._doMatch(
            refCat=trimmedRefCat,
            sourceCat=goodSourceCat,
            wcs=wcs,
            refFluxField=refFluxField,
            numUsableSources=numUsableSources,
            minMatchedPairs=minMatchedPairs,
            match_tolerance=match_tolerance,
            sourceFluxField=sourceFluxField,
            verbose=debug.verbose,
        )
        matches = doMatchReturn.matches
        match_tolerance = doMatchReturn.match_tolerance

        if len(matches) == 0:
            raise RuntimeError("Unable to match sources")

        self.log.info("Matched %d sources" % len(matches))
        if len(matches) < minMatchedPairs:
            self.log.warn("Number of matches is smaller than request")

        return pipeBase.Struct(
            matches=matches,
            usableSourceCat=goodSourceCat,
            match_tolerance=match_tolerance,
        )

    def _filterRefCat(self, refCat, refFluxField):
        """Sub-select a number of reference objects starting from the brightest
        and maxing out at the number specified by maxRefObjects in the config.

        No trimming is done if len(refCat) > config.maxRefObjects.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            Catalog of reference objects to trim.
        refFluxField : `str`
            field of refCat to use for flux
        Returns
        -------
        outCat : `lsst.afw.table.SimpleCatalog`
            Catalog trimmed to the number set in the task config from the
            brightest flux down.
        """
        # Find the flux cut that gives us the desired number of objects.
        if len(refCat) <= self.config.maxRefObjects:
            return refCat
        fluxArray = refCat.get(refFluxField)
        sortedFluxArray = fluxArray[fluxArray.argsort()]
        minFlux = sortedFluxArray[-(self.config.maxRefObjects + 1)]

        selected = (refCat.get(refFluxField) > minFlux)

        outCat = afwTable.SimpleCatalog(refCat.schema)
        outCat.reserve(self.config.maxRefObjects)
        outCat.extend(refCat[selected])

        return outCat

    @pipeBase.timeMethod
    def _doMatch(self, refCat, sourceCat, wcs, refFluxField, numUsableSources,
                 minMatchedPairs, match_tolerance, sourceFluxField, verbose):
        """Implementation of matching sources to position reference objects

        Unlike matchObjectsToSources, this method does not check if the sources
        are suitable.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            catalog of position reference objects that overlap an exposure
        sourceCat : `lsst.afw.table.SourceCatalog`
            catalog of sources found on the exposure
        wcs : `lsst.afw.geom.SkyWcs`
            estimated WCS of exposure
        refFluxField : `str`
            field of refCat to use for flux
        numUsableSources : `int`
            number of usable sources (sources with known centroid that are not
            near the edge, but may be saturated)
        minMatchedPairs : `int`
            minimum number of matches
        match_tolerance : `lsst.meas.astrom.MatchTolerancePessimistic`
            a MatchTolerance object containing variables specifying matcher
            tolerances and state from possible previous runs.
        sourceFluxField : `str`
            Name of the flux field in the source catalog.
        verbose : `bool`
            Set true to print diagnostic information to std::cout

        Returns
        -------
        result :
            Results struct with components:

            - ``matches`` : a list the matches found
              (`list` of `lsst.afw.table.ReferenceMatch`).
            - ``match_tolerance`` : MatchTolerance containing updated values from
              this fit iteration (`lsst.meas.astrom.MatchTolerancePessimistic`)
        """

        # Load the source and reference catalog as spherical points
        # in numpy array. We do this rather than relying on internal
        # lsst C objects for simplicity and because we require
        # objects contiguous in memory. We need to do these slightly
        # differently for the reference and source cats as they are
        # different catalog objects with different fields.
        src_array = np.empty((len(sourceCat), 4), dtype=np.float64)
        for src_idx, srcObj in enumerate(sourceCat):
            coord = wcs.pixelToSky(srcObj.getCentroid())
            theta = np.pi / 2 - coord.getLatitude().asRadians()
            phi = coord.getLongitude().asRadians()
            flux = srcObj[sourceFluxField]
            src_array[src_idx, :] = \
                self._latlong_flux_to_xyz_mag(theta, phi, flux)

        if match_tolerance.PPMbObj is None or \
           match_tolerance.autoMaxMatchDist is None:
            # The reference catalog is fixed per AstrometryTask so we only
            # create the data needed if this is the first step in the match
            # fit cycle.
            ref_array = np.empty((len(refCat), 4), dtype=np.float64)
            for ref_idx, refObj in enumerate(refCat):
                theta = np.pi / 2 - refObj.getDec().asRadians()
                phi = refObj.getRa().asRadians()
                flux = refObj[refFluxField]
                ref_array[ref_idx, :] = \
                    self._latlong_flux_to_xyz_mag(theta, phi, flux)
            # Create our matcher object.
            match_tolerance.PPMbObj = PessimisticPatternMatcherB(
                ref_array[:, :3], self.log)
            self.log.debug("Computing source statistics...")
            maxMatchDistArcSecSrc = self._get_pair_pattern_statistics(
                src_array)
            self.log.debug("Computing reference statistics...")
            maxMatchDistArcSecRef = self._get_pair_pattern_statistics(
                ref_array)
            maxMatchDistArcSec = np.max((
                self.config.minMatchDistPixels *
                wcs.getPixelScale().asArcseconds(),
                np.min((maxMatchDistArcSecSrc,
                        maxMatchDistArcSecRef))))
            match_tolerance.autoMaxMatchDist = afwgeom.Angle(
                maxMatchDistArcSec, afwgeom.arcseconds)

        # Set configurable defaults when we encounter None type or set
        # state based on previous run of AstrometryTask._matchAndFitWcs.
        if match_tolerance.maxShift is None:
            maxShiftArcseconds = (self.config.maxOffsetPix *
                                  wcs.getPixelScale().asArcseconds())
        else:
            # We don't want to clamp down too hard on the allowed shift so
            # we test that the smallest we ever allow is the pixel scale.
            maxShiftArcseconds = np.max(
                (match_tolerance.maxShift.asArcseconds(),
                 self.config.minMatchDistPixels *
                 wcs.getPixelScale().asArcseconds()))

        # If our tolerances are not set from a previous run, estimate a
        # starting tolerance guess from the statistics of patterns we can
        # create on both the source and reference catalog. We use the smaller
        # of the two.
        if match_tolerance.maxMatchDist is None:
            match_tolerance.maxMatchDist = match_tolerance.autoMaxMatchDist
        else:
            maxMatchDistArcSec = np.max(
                (self.config.minMatchDistPixels *
                 wcs.getPixelScale().asArcseconds(),
                 np.min((match_tolerance.maxMatchDist.asArcseconds(),
                         match_tolerance.autoMaxMatchDist.asArcseconds()))))

        # Make sure the data we are considering is dense enough to require
        # the consensus mode of the matcher. If not default to Optimistic
        # pattern matcher behavior. We enforce pessimistic mode if the
        # reference cat is sufficiently large, avoiding false positives.
        numConsensus = self.config.numPatternConsensus
        if len(refCat) < self.config.numRefRequireConsensus:
            minObjectsForConsensus = \
                self.config.numBrightStars + \
                self.config.numPointsForShapeAttempt
            if len(refCat) < minObjectsForConsensus or \
               len(sourceCat) < minObjectsForConsensus:
                numConsensus = 1

        self.log.debug("Current tol maxDist: %.4f arcsec" %
                       maxMatchDistArcSec)
        self.log.debug("Current shift: %.4f arcsec" %
                       maxShiftArcseconds)

        match_found = False
        # Start the iteration over our tolerances.
        for soften_dist in range(self.config.matcherIterations):
            if soften_dist == 0 and \
               match_tolerance.lastMatchedPattern is not None:
                # If we are on the first, most stringent tolerance,
                # and have already found a match, the matcher should behave
                # like an optimistic pattern matcher. Exiting at the first
                # match.
                run_n_consent = 1
            else:
                # If we fail or this is the first match attempt, set the
                # pattern consensus to the specified config value.
                run_n_consent = numConsensus
            # We double the match dist tolerance each round and add 1 to the
            # to the number of candidate spokes to check.
            matcher_struct = match_tolerance.PPMbObj.match(
                source_array=src_array,
                n_check=self.config.numPointsForShapeAttempt,
                n_match=self.config.numPointsForShape,
                n_agree=run_n_consent,
                max_n_patterns=self.config.numBrightStars,
                max_shift=maxShiftArcseconds,
                max_rotation=self.config.maxRotationDeg,
                max_dist=maxMatchDistArcSec * 2. ** soften_dist,
                min_matches=minMatchedPairs,
                pattern_skip_array=np.array(
                    match_tolerance.failedPatternList)
            )

            if soften_dist == 0 and \
               len(matcher_struct.match_ids) == 0 and \
               match_tolerance.lastMatchedPattern is not None:
                # If we found a pattern on a previous match-fit iteration and
                # can't find an optimistic match on our first try with the
                # tolerances as found in the previous match-fit,
                # the match we found in the last iteration was likely bad. We
                # append the bad match's index to the a list of
                # patterns/matches to skip on subsequent iterations.
                match_tolerance.failedPatternList.append(
                    match_tolerance.lastMatchedPattern)
                match_tolerance.lastMatchedPattern = None
                maxShiftArcseconds = \
                    self.config.maxOffsetPix * wcs.getPixelScale().asArcseconds()
            elif len(matcher_struct.match_ids) > 0:
                # Match found, save a bit a state regarding this pattern
                # in the match tolerance class object and exit.
                match_tolerance.maxShift = \
                    matcher_struct.shift * afwgeom.arcseconds
                match_tolerance.lastMatchedPattern = \
                    matcher_struct.pattern_idx
                match_found = True
                break

        # If we didn't find a match, exit early.
        if not match_found:
            return pipeBase.Struct(
                matches=[],
                match_tolerance=match_tolerance,
            )

        # The matcher returns all the nearest neighbors that agree between
        # the reference and source catalog. For the current astrometric solver
        # we need to remove as many false positives as possible before sending
        # the matches off to the solver. The low value of 100 and high value of
        # 2 are the low number of sigma and high respectively. The exact values
        # were found after testing on data of various reference/source
        # densities and astrometric distortion quality, specifically the
        # visits:  HSC (3358), DECam (406285, 410827),
        # CFHT (793169, 896070, 980526).
        distances_arcsec = np.degrees(matcher_struct.distances_rad) * 3600
        dist_cut_arcsec = np.max(
            (np.degrees(matcher_struct.max_dist_rad) * 3600,
             self.config.minMatchDistPixels * wcs.getPixelScale().asArcseconds()))

        # A match has been found, return our list of matches and
        # return.
        matches = []
        for match_id_pair, dist_arcsec in zip(matcher_struct.match_ids,
                                              distances_arcsec):
            if dist_arcsec < dist_cut_arcsec:
                match = afwTable.ReferenceMatch()
                match.first = refCat[int(match_id_pair[1])]
                match.second = sourceCat[int(match_id_pair[0])]
                # We compute the true distance along and sphere. This isn't
                # used in the WCS fitter however it is used in the unittest
                # to confirm the matches computed.
                match.distance = match.first.getCoord().separation(
                    match.second.getCoord()).asArcseconds()
                matches.append(match)

        return pipeBase.Struct(
            matches=matches,
            match_tolerance=match_tolerance,
        )

    def _latlong_flux_to_xyz_mag(self, theta, phi, flux):
        """Convert angles theta and phi and a flux into unit sphere
        x, y, z, and a relative magnitude.

        Takes in a afw catalog object and converts the catalog object RA, DECs
        to points on the unit sphere. Also converts the flux into a simple,
        non-zero-pointed magnitude for relative sorting.

        Parameters
        ----------
        theta : `float`
            Angle from the north pole (z axis) of the sphere
        phi : `float`
            Rotation around the sphere

        Return
        ------
        output_array : `numpy.ndarray`, (N, 4)
            Spherical unit vector x, y, z  with flux.
        """
        output_array = np.empty(4, dtype=np.float64)
        output_array[0] = np.sin(theta)*np.cos(phi)
        output_array[1] = np.sin(theta)*np.sin(phi)
        output_array[2] = np.cos(theta)
        if flux > 0:
            output_array[3] = -2.5 * np.log10(flux)
        else:
            # Set flux to a very faint mag if its for some reason it
            # does not exist
            output_array[3] = 99.

        return output_array

    def _get_pair_pattern_statistics(self, cat_array):
        """ Compute the tolerances for the matcher automatically by comparing
        pinwheel patterns as we would in the matcher.

        We test how similar the patterns we can create from a given set of
        objects by computing the spoke lengths for each pattern and sorting
        them from smallest to largest. The match tolerance is the average
        distance per spoke between the closest two patterns in the sorted
        spoke length space.

        Parameters
        ----------
        cat_array : `numpy.ndarray`, (N, 3)
            array of 3 vectors representing the x, y, z position of catalog
            objects on the unit sphere.

        Returns
        -------
        dist_tol : `float`
            Suggested max match tolerance distance calculated from comparisons
            between pinwheel patterns used in optimistic/pessimistic pattern
            matcher.
        """

        self.log.debug("Starting automated tolerance calculation...")

        # Create an empty array of all the patterns we possibly make
        # sorting from brightest to faintest.
        pattern_array = np.empty(
            (cat_array.shape[0] - self.config.numPointsForShape,
             self.config.numPointsForShape - 1))
        flux_args_array = np.argsort(cat_array[:, -1])

        # Sort our input array.
        tmp_sort_array = cat_array[flux_args_array]

        # Start making patterns.
        for start_idx in range(cat_array.shape[0] -
                               self.config.numPointsForShape):
            pattern_points = tmp_sort_array[start_idx:start_idx +
                                            self.config.numPointsForShape, :-1]
            pattern_delta = pattern_points[1:, :] - pattern_points[0, :]
            pattern_array[start_idx, :] = np.sqrt(
                pattern_delta[:, 0] ** 2 +
                pattern_delta[:, 1] ** 2 +
                pattern_delta[:, 2] ** 2)

            # When we store the length of each spoke in our pattern we
            # sort from shortest to longest so we have a defined space
            # to compare them in.
            pattern_array[start_idx, :] = pattern_array[
                start_idx, np.argsort(pattern_array[start_idx, :])]

        # Create a searchable tree object of the patterns and find
        # for any given pattern the closest pattern in the sorted
        # spoke length space.
        dist_tree = cKDTree(
            pattern_array[:, :(self.config.numPointsForShape - 1)])
        dist_nearest_array, ids = dist_tree.query(
            pattern_array[:, :(self.config.numPointsForShape - 1)], k=2)
        dist_nearest_array = dist_nearest_array[:, 1]
        dist_nearest_array.sort()

        # We use the two closest patterns to set our tolerance.
        dist_idx = 0
        dist_tol = (np.degrees(dist_nearest_array[dist_idx]) * 3600. /
                    (self.config.numPointsForShape - 1.))

        self.log.debug("Automated tolerance")
        self.log.debug("\tdistance/match tol: %.4f [arcsec]" % dist_tol)

        return dist_tol
