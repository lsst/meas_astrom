from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import sigmaclip

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwgeom
import lsst.afw.table as afwTable
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .matchOptimisticB import MatchTolerance

from .pessimistic_pattern_matcher_b_3D import PessimisticPatternMatcherB

__all__ = ["MatchPessimisticBTask", "MatchPessimisticBConfig",
           "MatchTolerancePessimistic"]


class MatchTolerancePessimistic(MatchTolerance):
    """Stores match tolerances for use in AstrometryTask and later
    iterations of the matcher.

    Attributes
    ----------
    maxMatchDist : lsst.afw.geom.Angle
    autoMaxMatchDist : lsst.afw.geom.Angle
    maxShift : lsst.afw.geom.Angle
    lastMatchedPattern : int
    failedPatternList : list of ints
    """

    def __init__(self, maxMatchDist=None, autoMaxMatchDist=None,
                 maxShift=None, lastMatchedPattern=None,
                 failedPatternList=None):
        """Construct a MatchPessimisticTolerance

        MatchPessimisticBTask relies on several state variables to be
        preserved over different iterations in the
        AstrometryTask.matchAndFitWcs loop of AstrometryTask.

        Parameters
        ----------
        maxMatchDist : afw.geom.Angle
            Current 2 sigma scatter from the previous matched wcs (if it
            exists. It is None if this is the first iteration.)
        autoMatxMatchDist : afw.geom.Angle
            Result of the automated match tolerance generation.
        maxShift  : afw.geom.Angle
            None for the first iteration or is the magnitude of the previous
            iteration's wcs shift.
        lastMatchedPattern : int
            Reference to the position in the magnitude sorted source array
            where a successful pattern match was found.
        failedPatternList : list of ints
            List of ints specifying indicies in the magnitude sourced source
            array to skip. These are skipped are previous iterations that are
            likely false positives due to the code having to soften after a
            pattern is matched.
        """
        self.maxMatchDist = maxMatchDist
        self.autoMaxMatchDist = autoMaxMatchDist
        self.maxShift = maxShift
        self.lastMatchedPattern = lastMatchedPattern
        if failedPatternList is None:
            self.failedPatternList = []


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
            "loosening the signal to noise cut in the matcherSourceSelector, "
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
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources for cross-matching. The default "
            "matcherSourceSelector removes objects with low S/N, bad "
            "saturated objects, edge objects, and interpolated objects.",
        default="matcherPessimistic"
    )

    def setDefaults(self):
        sourceSelector = self.sourceSelector["matcherPessimistic"]
        sourceSelector.setDefaults()

    def validate(self):
        pexConfig.Config.validate(self)
        if self.numPointsForShapeAttempt < self.numPointsForShape:
            raise ValueError("numPointsForShapeAttempt must be greater than "
                             "or equal to numPointsForShape.")


# The following block adds links to this task from the Task Documentation page.
# \addtogroup LSST_task_documentation
# \{
# \page measAstrom_MatchPessimisticBTask
# \ref MatchPessimisticBTask "MatchPessimisticBTask"
# Match sources to reference objects
# \}


class MatchPessimisticBTask(pipeBase.Task):
    """!Match sources to reference objects

    @anchor MatchPessimisticBTask_

    @section meas_astrom_MatchPessimisticB_Contents Contents

     - @ref meas_astrom_MatchPessimisticB_Purpose
     - @ref meas_astrom_MatchPessimisticB_Initialize
     - @ref meas_astrom_MatchPessimisticB_IO
     - @ref meas_astrom_MatchPessimisticB_Config
     - @ref meas_astrom_MatchPessimisticB_Example
     - @ref meas_astrom_MatchPessimisticB_Debug

    @section meas_astrom_MatchPessimisticB_Purpose  Description

    Match sources to reference objects. This is often done as a preliminary
    step to fitting an astrometric or photometric solution. For details about
    the matching algorithm see pessimistic_pattern_matcher_b_3D.py

    @section meas_astrom_MatchPessimisticB_Initialize   Task initialization

    @copydoc \_\_init\_\_

    @section meas_astrom_MatchPessimisticB_IO       Invoking the Task

    @copydoc matchObjectsToSources

    @section meas_astrom_MatchPessimisticB_Config       Configuration
    parameters

    See @ref MatchPessimisticBConfig

    To modify the tests for good sources for matching, create a new
    sourceSelector class in meas_algorithms and use it in the config.

    @section meas_astrom_MatchPessimisticB_Example  A complete example of
    using MatchPessimisticBTask

    MatchPessimisticBTask is a subtask of AstrometryTask, which is called by
    PhotoCalTask. See \ref meas_photocal_photocal_Example.

    @section meas_astrom_MatchPessimisticB_Debug        Debug variables

    The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink
    interface supports a flag @c -d to import @b debug.py from your
    @c PYTHONPATH; see @ref baseDebug for more about @b debug.py files.

    The available variables in MatchPessimisticBTask are:
    <DL>
      <DT> @c verbose (bool)
      <DD> If True then the matcher prints debug messages to stdout
    </DL>

    To investigate the @ref meas_astrom_MatchPessimisticB_Debug, put something
    like
    @code{.py}
        import lsstDebug
        def DebugInfo(name):
            # N.b. lsstDebug.Info(name) would call us recursively
            debug = lsstDebug.getInfo(name)
            if name == "lsst.pipe.tasks.astrometry":
                debug.verbose = True

            return debug

        lsstDebug.Info = DebugInfo
    @endcode
    into your debug.py file and run this task with the @c --debug flag.
    """

    ConfigClass = MatchPessimisticBConfig
    _DefaultName = "matchObjectsToSources"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("sourceSelector")

    @pipeBase.timeMethod
    def matchObjectsToSources(self, refCat, sourceCat, wcs, refFluxField,
                              match_tolerance=None):
        """!Match sources to position reference stars

        @param[in] refCat  catalog of reference objects that overlap the
        exposure; reads fields for:
            - coord
            - the specified flux field
        @param[in] sourceCat  catalog of sources found on an exposure;
            Please check the required fields of your specified source selector
            that the correct flags are present.
        @param[in] wcs  estimated WCS
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] match_tolerance is a MatchTolerance class object or None.
            This this class is used to communicate state between AstrometryTask
            and MatcherTask. AstrometryTask will also set the MatchTolerance
            class variable maxMatchDist based on the scatter AstrometryTask has
            found after fitting for the wcs.
        @return an lsst.pipe.base.Struct with fields:
        - matches  a list of matches, each instance of
            lsst.afw.table.ReferenceMatch
        - usableSourcCat  a catalog of sources potentially usable for
            matching.
        - match_tolerance a MatchTolerance object containing the resulting
            state variables from the match.
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        # If we get an empty tolerance struct create the variables we need for
        # this matcher.
        if match_tolerance is None:
            match_tolerance = MatchTolerancePessimistic()

        # usableSourceCat: sources that are good but may be saturated
        numSources = len(sourceCat)
        selectedSources = self.sourceSelector.selectSources(sourceCat)
        goodSourceCat = selectedSources.sourceCat
        numUsableSources = len(goodSourceCat)
        self.log.info("Purged %d sources, leaving %d good sources" %
                      (numSources - numUsableSources, numUsableSources))

        if len(goodSourceCat) == 0:
            raise pipeBase.TaskError("No sources are good")

        # avoid accidentally using sourceCat; use goodSourceCat from now on
        del sourceCat

        minMatchedPairs = min(self.config.minMatchedPairs,
                              int(self.config.minFracMatchedPairs *
                                  min([len(refCat), len(goodSourceCat)])))

        doMatchReturn = self._doMatch(
            refCat=refCat,
            sourceCat=goodSourceCat,
            wcs=wcs,
            refFluxField=refFluxField,
            numUsableSources=numUsableSources,
            minMatchedPairs=minMatchedPairs,
            match_tolerance=match_tolerance,
            sourceFluxField=self.sourceSelector.fluxField,
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

    @pipeBase.timeMethod
    def _doMatch(self, refCat, sourceCat, wcs, refFluxField, numUsableSources,
                 minMatchedPairs, match_tolerance, sourceFluxField, verbose):
        """!Implementation of matching sources to position reference stars

        Unlike matchObjectsToSources, this method does not check if the sources
        are suitable.

        @param[in] refCat  catalog of position reference stars that overlap an
            exposure
        @param[in] sourceCat  catalog of sources found on the exposure
        @param[in] wcs  estimated WCS of exposure
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] numUsableSources  number of usable sources (sources with
            known centroid that are not near the edge, but may be saturated)
        @param[in] minMatchedPairs  minimum number of matches
        @param[in] match_tolerance a MatchTolerance object containing
            variables specifying matcher tolerances and state from possible
            previous runs.
        @param[in] sourceInfo  SourceInfo for the sourceCat
        @param[in] verbose  true to print diagnostic information to std::cout

        @return a list of matches, an instance of
            lsst.afw.table.ReferenceMatch, a MatchTolerance object
        """

        # Load the source and reference catalog as spherical points
        # in numpy array. We do this rather than relying on internal
        # lsst C objects for simplicity and because we require
        # objects contiguous in memory. We need to do these slightly
        # differently for the reference and source cats as they are
        # different catalog objects with different fields.
        ref_array = np.empty((len(refCat), 4), dtype=np.float64)
        for ref_idx, refObj in enumerate(refCat):
            theta = np.pi / 2 - refObj.getDec().asRadians()
            phi = refObj.getRa().asRadians()
            flux = refObj[refFluxField]
            ref_array[ref_idx, :] = \
                self._latlong_flux_to_xyz_mag(theta, phi, flux)

        src_array = np.empty((len(sourceCat), 4), dtype=np.float64)
        for src_idx, srcObj in enumerate(sourceCat):
            coord = wcs.pixelToSky(srcObj.getCentroid())
            theta = np.pi / 2 - coord.getLatitude().asRadians()
            phi = coord.getLongitude().asRadians()
            flux = srcObj.getPsfFlux()
            src_array[src_idx, :] = \
                self._latlong_flux_to_xyz_mag(theta, phi, flux)

        # Set configurable defaults when we encounter None type or set
        # state based on previous run of AstrometryTask._matchAndFitWcs.
        if match_tolerance.maxShift is None:
            maxShiftArcseconds = (self.config.maxOffsetPix *
                                  wcs.pixelScale().asArcseconds())
        else:
            # We don't want to clamp down too hard on the allowed shift so
            # we test that the smallest we ever allow is the pixel scale.
            maxShiftArcseconds = np.max(
                (match_tolerance.maxShift.asArcseconds(),
                 self.config.minMatchDistPixels *
                 wcs.pixelScale().asArcseconds()))

        # If our tolerances are not set from a previous run, estimate a
        # starting tolerance guess from the statistics of patterns we can
        # create on both the source and reference catalog. We use the smaller
        # of the two.
        if match_tolerance.maxMatchDist is None:
            self.log.debug("Computing source statistics...")
            maxMatchDistArcSecSrc = self._get_pair_pattern_statistics(
                src_array)
            self.log.debug("Computing reference statistics...")
            maxMatchDistArcSecRef = self._get_pair_pattern_statistics(
                ref_array)
            maxMatchDistArcSec = np.min((maxMatchDistArcSecSrc,
                                         maxMatchDistArcSecRef))
            match_tolerance.autoMaxDist = afwgeom.Angle(maxMatchDistArcSec,
                                                        afwgeom.arcseconds)
        else:
            maxMatchDistArcSec = np.max(
                (self.config.minMatchDistPixels *
                 wcs.pixelScale().asArcseconds(),
                 np.min((match_tolerance.maxMatchDist.asArcseconds(),
                         match_tolerance.autoMaxDist.asArcseconds()))))

        # Make sure the data we are considering is dense enough to require
        # the consensus mode of the matcher. If not default to Optimistic
        # pattern matcher behavior.
        numConsensus = self.config.numPatternConsensus
        minObjectsForConsensus = \
            self.config.numBrightStars + self.config.numPointsForShapeAttempt
        if ref_array.shape[0] < minObjectsForConsensus or \
           src_array.shape[0] < minObjectsForConsensus:
            numConsensus = 1

        self.log.debug("Current tol maxDist: %.4f arcsec" %
                       maxMatchDistArcSec)
        self.log.debug("Current shift: %.4f arcsec" %
                       maxShiftArcseconds)

        # Create our matcher object.
        pyPPMb = PessimisticPatternMatcherB(ref_array[:, :3], self.log)

        match_found = False
        # Start the iteration over our tolerances.
        for try_idx in range(self.config.matcherIterations):
            for soften_idx in range(self.config.matcherIterations):
                if try_idx == 0 and soften_idx == 0 and \
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
                matcher_struct = pyPPMb.match(
                    source_array=src_array,
                    n_check=self.config.numPointsForShapeAttempt + try_idx,
                    n_match=self.config.numPointsForShape,
                    n_agree=run_n_consent,
                    max_n_patterns=self.config.numBrightStars,
                    max_shift=maxShiftArcseconds,
                    max_rotation=self.config.maxRotationDeg,
                    max_dist=maxMatchDistArcSec * 2. ** soften_idx,
                    min_matches=minMatchedPairs,
                    pattern_skip_array=np.array(
                        match_tolerance.failedPatternList)
                )

                if try_idx == 0 and soften_idx == 0 and \
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
                        self.config.maxOffsetPix * wcs.pixelScale().asArcseconds()
                    continue
                elif len(matcher_struct.match_ids) > 0:
                    # Match found, save a bit a state regarding this pattern
                    # in the match tolerance class object and exit.
                    match_tolerance.maxShift = \
                        matcher_struct.shift * afwgeom.arcseconds
                    match_tolerance.lastMatchedPattern = \
                        matcher_struct.pattern_idx
                    match_found = True
                    break
            if match_found:
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
        # the matches off to the solver.
        distances_arcsec = np.degrees(matcher_struct.distances_rad) * 3600 
        clip_max_dist = np.max(
            (sigmaclip(distances_arcsec, low=100, high=2)[-1],
             self.config.minMatchDistPixels * wcs.pixelScale().asArcseconds())
        )
        # We pick the largest of the the unsoftened maxMatchDistArcSec or
        # a user specified distance in pixels. This prevents the
        # AstrometryTask._matchAndFitWCS from over-fitting to a small number of
        # objects and also allows the WCS fitter to bring in more matches as
        # the WCS fit improves.
        if not np.isfinite(clip_max_dist):
            clip_max_dist = maxMatchDistArcSec
        
        if clip_max_dist < maxMatchDistArcSec and \
           len(distances_arcsec[distances_arcsec < clip_max_dist]) < \
           minMatchedPairs:
            dist_cut_arcsec = maxMatchDistArcSec
        else:
            dist_cut_arcsec = np.min((clip_max_dist, maxMatchDistArcSec))

        # A match has been found, return our list of matches and
        # return.
        matches = []
        for match_id_pair, dist_arcsec in zip(matcher_struct.match_ids,
                                              distances_arcsec):
            if dist_arcsec < dist_cut_arcsec:
                match = afwTable.ReferenceMatch()
                match.first = refCat[match_id_pair[1]]
                match.second = sourceCat[match_id_pair[0]]
                # We compute the true distance along and sphere instead
                # and store it in units of arcseconds. The previous
                # distances we used were approximate.
                match.distance = match.first.getCoord().angularSeparation(
                    match.second.getCoord()).asArcseconds()
                matches.append(match)

        return pipeBase.Struct(
            matches=matches,
            match_tolerance=match_tolerance,
        )

    def _latlong_flux_to_xyz_mag(self, theta, phi, flux):
        r"""Convert angles theta and phi and a flux into unit sphere
        x, y, z, and a relative magnitude.

        Takes in a afw catalog object and converts the catalog object RA, DECs
        to points on the unit sphere. Also converts the flux into a simple,
        non-zero-pointed magnitude for relative sorting.

        Parameters
        ----------
        theta : float
            Angle from the north pole (z axis) of the sphere
        phi : float
            Rotation around the sphere

        Return
        ------
        float array
            Spherical unit vector x, y, z  on the unit-sphere.
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
        cat_array : float array
            array of 3 vectors representing the x, y, z position of catalog
            objects on the unit sphere.

        Returns
        -------
        float
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
