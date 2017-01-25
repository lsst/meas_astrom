from __future__ import absolute_import, division, print_function
from builtins import range
from builtins import object
import math

import numpy as np
from scipy.spatial import cKDTree

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .setMatchDistance import setMatchDistance
from .optimistic_pattern_matcher_b_3D import OptimisticPatternMatcherB

__all__ = ["MatchOptimisticBTask", "MatchOptimisticBConfig"]
__deg_to_rad__ = np.pi/180.


class MatchOptimisticBConfig(pexConfig.Config):
    """Configuration for MatchOptimisticBTask
    """
    maxMatchDistArcSec = pexConfig.RangeField(
        doc="Maximum separation between reference objects and sources "
        "beyond which they will not be considered a match (arcsec)",
        dtype=float,
        default=15,
        min=0,
    )
    maxOffsetPix = pexConfig.RangeField(
        doc="Max offset in pixels to work with HSC.",
        dtype=float,
        default=750,
        min=0,
    )
    maxAngTol = pexConfig.RangeField(
        doc="Maximum angle allowed for pattern in degress",
        dtype=float,
        default=0.5,
        min=0,
    )
    numPatterns = pexConfig.RangeField(
        doc="Number of patterns to attempt before exiting",
        dtype=int,
        default=100,
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
        default=0.15,
        min=0,
        max=1,
    )
    maxShift = pexConfig.RangeField(
        doc="Maximum allowed shift of WCS, due to matching (arcsec)",
        dtype=int,
        default=50,
        max=400,
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
    numPointsForShapeAttempt = pexConfig.Field(
        doc="number of points to try to match for a shape",
        dtype=int,
        default=9,
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
    """!Match sources to reference objects

    @anchor MatchOptimisticBTask_

    @section meas_astrom_matchOptimisticB_Contents Contents

     - @ref meas_astrom_matchOptimisticB_Purpose
     - @ref meas_astrom_matchOptimisticB_Initialize
     - @ref meas_astrom_matchOptimisticB_IO
     - @ref meas_astrom_matchOptimisticB_Config
     - @ref meas_astrom_matchOptimisticB_Example
     - @ref meas_astrom_matchOptimisticB_Debug

    @section meas_astrom_matchOptimisticB_Purpose  Description

    Match sources to reference objects. This is often done as a preliminary step to fitting an astrometric
    or photometric solution. For details about the matching algorithm see matchOptimisticB.h

    @section meas_astrom_matchOptimisticB_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_matchOptimisticB_IO       Invoking the Task

    @copydoc matchObjectsToSources

    @section meas_astrom_matchOptimisticB_Config       Configuration parameters

    See @ref MatchOptimisticBConfig

    To modify the tests for usable sources and good sources, subclass SourceInfo and
    set MatchOptimisticBTask.SourceInfoClass to your subclass.

    @section meas_astrom_matchOptimisticB_Example  A complete example of using MatchOptimisticBTask

    MatchOptimisticBTask is a subtask of AstrometryTask, which is called by PhotoCalTask.
    See \ref meas_photocal_photocal_Example.

    @section meas_astrom_matchOptimisticB_Debug        Debug variables

    The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink interface supports a
    flag @c -d to import @b debug.py from your @c PYTHONPATH; see @ref baseDebug for more about
    @b debug.py files.

    The available variables in MatchOptimisticBTask are:
    <DL>
      <DT> @c verbose (bool)
      <DD> If True then the matcher prints debug messages to stdout
    </DL>

    To investigate the @ref meas_astrom_matchOptimisticB_Debug, put something like
    @code{.py}
        import lsstDebug
        def DebugInfo(name):
            debug = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
            if name == "lsst.pipe.tasks.astrometry":
                debug.verbose = True

            return debug

        lsstDebug.Info = DebugInfo
    @endcode
    into your debug.py file and run this task with the @c --debug flag.
    """
    ConfigClass = MatchOptimisticBConfig
    _DefaultName = "matchObjectsToSources"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("sourceSelector")

    def filterStars(self, refCat):
        """Extra filtering pass; subclass if desired
        """
        return refCat

    @pipeBase.timeMethod
    def matchObjectsToSources(self, refCat, sourceCat, wcs, refFluxField, maxShift=None,
                              maxMatchDist=None):
        """!Match sources to position reference stars

        @param[in] refCat  catalog of reference objects that overlap the exposure; reads fields for:
            - coord
            - the specified flux field
        @param[in] sourceCat  catalog of sources found on an exposure; reads fields for:
            - centroid
            - centroid flag
            - edge flag
            - saturated flag
            - aperture flux, if found, else PSF flux
        @param[in] wcs  estimated WCS
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] maxMatchDist  maximum on-sky distance between reference objects and sources
            (an lsst.afw.geom.Angle); if specified then the smaller of config.maxMatchDistArcSec or
            maxMatchDist is used; if None then config.maxMatchDistArcSec is used
        @return an lsst.pipe.base.Struct with fields:
        - matches  a list of matches, each instance of lsst.afw.table.ReferenceMatch
        - usableSourcCat  a catalog of sources potentially usable for matching.
            For this fitter usable sources include unresolved sources not too near the edge.
            It includes saturated sources, even those these are removed from the final match list,
            because saturated sources may be used to determine the match list.
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        preNumObj = len(refCat)
        refCat = self.filterStars(refCat)
        numRefObj = len(refCat)

        if self.log:
            self.log.info("filterStars purged %d reference stars, leaving %d stars" %
                          (preNumObj - numRefObj, numRefObj))

        # usableSourceCat: sources that are good but may be saturated
        numSources = len(sourceCat)
        selectedSources = self.sourceSelector.selectSources(sourceCat)
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
        usableMatches, resShift = self._doMatch(
            refCat=refCat,
            sourceCat=usableSourceCat,
            wcs=wcs,
            refFluxField=refFluxField,
            numUsableSources=numUsableSources,
            minMatchedPairs=minMatchedPairs,
            maxShift=maxShift,
            maxMatchDist=maxMatchDist,
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
            resShift=resShift,
        )

    def _getIsGoodKeys(self, schema):
        self.edgeKey = schema["base_PixelFlags_flag_edge"].asKey()
        self.interpolatedCenterKey = schema["base_PixelFlags_flag_interpolatedCenter"].asKey()
        self.saturatedKey = schema["base_PixelFlags_flag_saturated"].asKey()

    def _isGoodTest(self, source):
        """
        This is a hard coded version of the isGood flag from the old SourceInfo class that used to be
        part of this class. This is done current as the API for sourceSelector does not currently
        support matchLists.
        """
        return (not source.get(self.edgeKey) and
                not source.get(self.interpolatedCenterKey) and
                not source.get(self.saturatedKey))
        
    @pipeBase.timeMethod
    def _doMatch(self, refCat, sourceCat, wcs, refFluxField, numUsableSources, minMatchedPairs,
                 maxMatchDist, maxShift, sourceInfo, verbose):
        """!Implementation of matching sources to position reference stars

        Unlike matchObjectsToSources, this method does not check if the sources are suitable.

        @param[in] refCat  catalog of position reference stars that overlap an exposure
        @param[in] sourceCat  catalog of sources found on the exposure
        @param[in] wcs  estimated WCS of exposure
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] numUsableSources  number of usable sources (sources with known centroid
            that are not near the edge, but may be saturated)
        @param[in] minMatchedPairs  minimum number of matches
        @param[in] maxMatchDist  maximum on-sky distance between reference objects and sources
            (an lsst.afw.geom.Angle); if specified then the smaller of config.maxMatchDistArcSec or
            maxMatchDist is used; if None then config.maxMatchDistArcSec is used
        @param[in] sourceInfo  SourceInfo for the sourceCat
        @param[in] verbose  true to print diagnostic information to std::cout

        @return a list of matches, an instance of lsst.afw.table.ReferenceMatch
        """
        if maxShift is None:
            maxShift = self.config.maxShift
        else:
            maxShift = max(0.2, min(maxShift, self.config.maxShift))
        # if maxMatchDist is None:
        #     maxMatchDistArcSec = self.config.maxMatchDistArcSec
        # else:
        #     maxMatchDistArcSec = min(maxMatchDist.asArcseconds(), self.config.maxMatchDistArcSec)
        # max_ang_tol = np.min((self.config.maxAngTol,
        #                       np.arctan(maxMatchDistArcSec/(0.2*2048*np.sqrt(2)))/__deg_to_rad__))
        max_rotation = self.config.maxRotationDeg

        ref_array = np.empty((len(refCat), 4))
        for ref_idx, refObj in enumerate(refCat):
            theta = np.pi/2 - refObj.getDec().asRadians()
            phi = refObj.getRa().asRadians()
            ref_array[ref_idx, 0] = np.sin(theta)*np.cos(phi)
            ref_array[ref_idx, 1] = np.sin(theta)*np.sin(phi)
            ref_array[ref_idx, 2] = np.cos(theta)
            ref_array[ref_idx, 3] = -2.5*np.log10(
                np.where(refObj[refFluxField] > 0,
                         refObj[refFluxField], 10**-32))

        src_array = np.empty((len(sourceCat), 4))
        for src_idx, srcObj in enumerate(sourceCat):
            coord = wcs.pixelToSky(srcObj.getCentroid())
            tmp_ra = coord.getLongitude().asRadians()
            tmp_dec = coord.getLatitude().asRadians()
            theta = np.pi/2 - tmp_dec
            phi = tmp_ra
            src_array[src_idx, 0] = np.sin(theta)*np.cos(phi)
            src_array[src_idx, 1] = np.sin(theta)*np.sin(phi)
            src_array[src_idx, 2] = np.cos(theta)
            src_array[src_idx, 3] = -2.5*np.log10(
                np.where(srcObj.getPsfFlux() > 0,
                         srcObj.getPsfFlux(), 10**-32))

        if maxMatchDist is None:
            if len(src_array) > len(ref_array):
                maxMatchDistArcSec, max_ang_tol = self._get_pair_pattern_statistics(
                    src_array)
            else:
                maxMatchDistArcSec, max_ang_tol = self._get_pair_pattern_statistics(
                    ref_array)
        else:
            maxMatchDistArcSec = min(maxMatchDist.asArcseconds(), self.config.maxMatchDistArcSec)
            max_ang_tol = np.min((self.config.maxAngTol,
                                 np.arctan(maxMatchDistArcSec/(0.2*1024*np.sqrt(2)))/__deg_to_rad__))

        pyOPMb = OptimisticPatternMatcherB(
            reference_catalog=ref_array, max_rotation_theta=maxShift/3600.,
            max_rotation_phi=max_rotation, dist_tol=maxMatchDistArcSec/3600.,
            max_dist_cand=100000, ang_tol=max_ang_tol,
            max_match_dist=np.min((self.config.maxMatchDistArcSec/3600.,
                                   2*maxMatchDistArcSec/3600.)),
            min_matches=minMatchedPairs, max_n_patterns=self.config.numPatterns)

        current_shift = None
        match_id_list = []
        dist_array = []
        for try_idx in xrange(4):
            match_id_list, dist_array = pyOPMb.match(src_array, self.config.numPointsForShapeAttempt + try_idx,
                                                     self.config.numPointsForShape)
            if len(match_id_list) > 0:
                current_shift = np.arccos(pyOPMb._cos_theta)*3600/__deg_to_rad__
                break
            else:
                maxShift *= 2
                maxShift = min((400., maxShift))
                maxMatchDistArcSec *= 2
                max_ang_tol *= 2
                max_rotation *=2
                pyOPMb._max_cos_theta = np.cos(maxShift/3600.*__deg_to_rad__)
                pyOPMb._max_cos_phi_sq = np.cos(max_rotation*__deg_to_rad__)**2
                pyOPMb._dist_tol = maxMatchDistArcSec/3600.*__deg_to_rad__
                pyOPMb._max_match_dist = maxMatchDistArcSec/3600.*__deg_to_rad__
                pyOPMb._ang_tol = max_ang_tol*__deg_to_rad__

        matches = afwTable.ReferenceMatchVector()
        for match_ids, dist in zip(match_id_list, dist_array):
            match = afwTable.ReferenceMatch()
            match.first = refCat.get(match_ids[1])
            match.second = sourceCat.get(match_ids[0])
            angular_sep = match.first.getCoord().angularSeparation(match.second.getCoord())
            match.distance = angular_sep.asRadians()
            matches.append(match)

        return matches, current_shift

    def _get_pair_pattern_statistics(self, cat_array):

        print("Starting automated tolerance calculation...")

        # Currently hard coded for a 6 point pattern.
        pattern_array = np.empty((cat_array.shape[0] - 6, 9))
        flux_args_array = np.argsort(cat_array[:, -1])

        tmp_sort_array = cat_array[flux_args_array]

        for start_idx in xrange(cat_array.shape[0] - 6):
            pattern_points = tmp_sort_array[start_idx:start_idx + 6, :-1]
            pattern_delta = pattern_points[1:, :] - pattern_points[0, :]
            pattern_array[start_idx, :5] = np.sqrt(pattern_delta[:, 0] ** 2 +
                                                  pattern_delta[:, 1] ** 2 +
                                                  pattern_delta[:, 2] ** 2)
            tmp_pattern_dot = (np.dot(pattern_delta[1:], pattern_delta[0]) /
                               (pattern_array[start_idx, 1:5] *
                                pattern_array[start_idx, 0]))
            tmp_pattern_cross = np.empty((4, 3))
            for idx in xrange(1, 5):
                tmp_pattern_cross[idx - 1, :] = (
                    np.cross(pattern_delta[idx], pattern_delta[0]) /
                    (pattern_array[start_idx, idx] *
                     pattern_array[start_idx, 0]))
            pattern_dot_cross = np.dot(tmp_pattern_cross,
                                       pattern_points[0])
            pattern_array[start_idx, 5:] = (np.sign(pattern_dot_cross) *
                                            np.arccos(tmp_pattern_dot))
            pattern_array[start_idx, :5] = pattern_array[
                start_idx, np.argsort(pattern_array[start_idx, :5])]
            pattern_array[start_idx, 5:] = pattern_array[
                start_idx, np.argsort(pattern_array[start_idx, 5:])]

        dist_tree = cKDTree(pattern_array[:, :5])
        theta_tree = cKDTree(pattern_array[:, 5:])

        dist_nearest_array, ids = dist_tree.query(pattern_array[:, :5], k=2)
        theta_nearest_array, ids = theta_tree.query(pattern_array[:, 5:], k=2)
        dist_nearest_array = dist_nearest_array[:, 1]
        theta_nearest_array = theta_nearest_array[:, 1]
        dist_nearest_array.sort()
        theta_nearest_array.sort()

        dist_idx = np.min((
            np.max((np.int64(np.ceil(dist_nearest_array.shape[0]*0.01)), 1)),
            dist_nearest_array.shape[0] - 1))
        theta_idx = np.min((
            np.max((np.int64(np.ceil(theta_nearest_array.shape[0]*0.01)), 1)),
            theta_nearest_array.shape[0] - 1))
        print("Index: %i %i" % (dist_idx, theta_idx))
        dist_tol = dist_nearest_array[dist_idx] * 3600. * (180. / np.pi) / np.sqrt(5.)
        theta_tol = theta_nearest_array[theta_idx] * (180. / np.pi) / np.sqrt(4.)

        print("New tolerances")
        print("\tdistance tol: %.4f [arcsec]" % dist_tol)
        print("\ttheta tol: %.4f [deg]" % theta_tol)

        return dist_tol, theta_tol
