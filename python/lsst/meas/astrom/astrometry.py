#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import absolute_import, division, print_function
from builtins import range

from lsst.afw.image.utils import getDistortedWcs
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .matchOptimisticB import MatchOptimisticBTask
from .fitTanSipWcs import FitTanSipWcsTask
from .display import displayAstrometry
from .astromLib import makeMatchStatistics
from .createMatchMetadata import createMatchMetadata


class AstrometryConfig(pexConfig.Config):
    matcher = pexConfig.ConfigurableField(
        target=MatchOptimisticBTask,
        doc="reference object/source matcher",
    )
    wcsFitter = pexConfig.ConfigurableField(
        target=FitTanSipWcsTask,
        doc="WCS fitter",
    )
    forceKnownWcs = pexConfig.Field(
        dtype=bool,
        doc="If True then load reference objects and match sources but do not fit a WCS; " +
        " this simply controls whether 'run' calls 'solve' or 'loadAndMatch'",
        default=False,
    )
    maxIter = pexConfig.RangeField(
        doc="maximum number of iterations of match sources and fit WCS" +
        "ignored if not fitting a WCS",
        dtype=int,
        default=3,
        min=1,
    )
    matchDistanceSigma = pexConfig.RangeField(
        doc="the maximum match distance is set to "
        " mean_match_distance + matchDistanceSigma*std_dev_match_distance; " +
        "ignored if not fitting a WCS",
        dtype=float,
        default=2,
        min=0,
    )
    minMatchDistanceArcSec = pexConfig.RangeField(
        doc="the match distance below which further iteration is pointless (arcsec); "
        "ignored if not fitting a WCS",
        dtype=float,
        default=0.001,
        min=0,
    )

# The following block adds links to this task from the Task Documentation page.
## \addtogroup LSST_task_documentation
## \{
## \page measAstrom_astrometryTask
## \ref AstrometryTask_ "AstrometryTask"
##      Match an input source catalog with objects from a reference catalog and solve for the WCS
## \}


class AstrometryTask(pipeBase.Task):
    """!Match an input source catalog with objects from a reference catalog and solve for the WCS

    @anchor AstrometryTask_

    @section meas_astrom_astrometry_Contents Contents

     - @ref meas_astrom_astrometry_Purpose
     - @ref meas_astrom_astrometry_Initialize
     - @ref meas_astrom_astrometry_IO
     - @ref meas_astrom_astrometry_Config
     - @ref meas_astrom_astrometry_Example
     - @ref meas_astrom_astrometry_Debug

    @section meas_astrom_astrometry_Purpose  Description

    Match input sourceCat with a reference catalog and solve for the Wcs

    There are three steps, each performed by different subtasks:
    - Find position reference stars that overlap the exposure
    - Match sourceCat to position reference stars
    - Fit a WCS based on the matches

    @section meas_astrom_astrometry_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_astrometry_IO       Invoking the Task

    @copydoc run

    @copydoc loadAndMatch

    @section meas_astrom_astrometry_Config       Configuration parameters

    See @ref AstrometryConfig

    @section meas_astrom_astrometry_Example  A complete example of using AstrometryTask

    See \ref meas_photocal_photocal_Example.

    @section meas_astrom_astrometry_Debug        Debug variables

    The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink interface supports a
    flag @c -d to import @b debug.py from your @c PYTHONPATH; see @ref baseDebug for more about
    @b debug.py files.

    The available variables in AstrometryTask are:
    <DL>
      <DT> @c display (bool)
      <DD> If True display information at three stages: after finding reference objects,
        after matching sources to reference objects, and after fitting the WCS; defaults to False
      <DT> @c frame (int)
      <DD> ds9 frame to use to display the reference objects; the next two frames are used
            to display the match list and the results of the final WCS; defaults to 0
    </DL>

    To investigate the @ref meas_astrom_astrometry_Debug, put something like
    @code{.py}
    import lsstDebug
    def DebugInfo(name):
        debug = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
        if name == "lsst.meas.astrom.astrometry":
            debug.display = True

        return debug

    lsstDebug.Info = DebugInfo
    @endcode
    into your debug.py file and run this task with the @c --debug flag.
    """
    ConfigClass = AstrometryConfig
    _DefaultName = "astrometricSolver"

    def __init__(self, refObjLoader, schema=None, **kwargs):
        """!Construct an AstrometryTask

        @param[in] refObjLoader A reference object loader object
        @param[in] schema  ignored; available for compatibility with an older astrometry task
        @param[in] kwargs  additional keyword arguments for pipe_base Task.\_\_init\_\_
        """
        pipeBase.Task.__init__(self, **kwargs)
        self.refObjLoader = refObjLoader
        self.makeSubtask("matcher")
        self.makeSubtask("wcsFitter")

    @pipeBase.timeMethod
    def run(self, exposure, sourceCat):
        """!Load reference objects, match sources and optionally fit a WCS

        This is a thin layer around solve or loadAndMatch, depending on config.forceKnownWcs

        @param[in,out] exposure  exposure whose WCS is to be fit
            The following are read only:
            - bbox
            - calib (may be absent)
            - filter (may be unset)
            - detector (if wcs is pure tangent; may be absent)
            The following are updated:
            - wcs (the initial value is used as an initial guess, and is required)
        @param[in] sourceCat  catalog of sources detected on the exposure (an lsst.afw.table.SourceCatalog)
        @return an lsst.pipe.base.Struct with these fields:
        - refCat  reference object catalog of objects that overlap the exposure (with some margin)
            (an lsst::afw::table::SimpleCatalog)
        - matches  list of reference object/source matches (an lsst.afw.table.ReferenceMatchVector)
        - scatterOnSky  median on-sky separation between reference objects and sources in "matches"
            (an lsst.afw.geom.Angle), or None if config.forceKnownWcs True
        - matchMeta  metadata needed to unpersist matches (an lsst.daf.base.PropertyList)
        """
        if self.config.forceKnownWcs:
            res = self.loadAndMatch(exposure=exposure, sourceCat=sourceCat)
            res.scatterOnSky = None
        else:
            res = self.solve(exposure=exposure, sourceCat=sourceCat)
        return res

    @pipeBase.timeMethod
    def loadAndMatch(self, exposure, sourceCat):
        """!Load reference objects overlapping an exposure and match to sources detected on that exposure

        @param[in] exposure  exposure that the sources overlap
        @param[in] sourceCat  catalog of sources detected on the exposure (an lsst.afw.table.SourceCatalog)

        @return an lsst.pipe.base.Struct with these fields:
        - refCat  reference object catalog of objects that overlap the exposure (with some margin)
            (an lsst::afw::table::SimpleCatalog)
        - matches  list of reference object/source matches (an lsst.afw.table.ReferenceMatchVector)
        - matchMeta  metadata needed to unpersist matches (an lsst.daf.base.PropertyList)

        @note ignores config.forceKnownWcs, config.maxIter, config.matchDistanceSigma
            and config.minMatchDistanceArcSec
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        matchMeta = createMatchMetadata(exposure, border=self.refObjLoader.config.pixelMargin)
        expMd = self._getExposureMetadata(exposure)

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            calib=expMd.calib,
        )

        matchRes = self.matcher.matchObjectsToSources(
            refCat=loadRes.refCat,
            sourceCat=sourceCat,
            wcs=expMd.wcs,
            refFluxField=loadRes.fluxField,
            maxMatchDist=None,
        )

        distStats = self._computeMatchStatsOnSky(matchRes.matches)
        self.log.info(
            "Found %d matches with scatter = %0.3f +- %0.3f arcsec; " %
            (len(matchRes.matches), distStats.distMean.asArcseconds(), distStats.distStdDev.asArcseconds())
        )

        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=loadRes.refCat,
                sourceCat=sourceCat,
                matches=matchRes.matches,
                exposure=exposure,
                bbox=expMd.bbox,
                frame=frame,
                title="Matches",
            )

        return pipeBase.Struct(
            refCat=loadRes.refCat,
            matches=matchRes.matches,
            matchMeta=matchMeta,
        )

    @pipeBase.timeMethod
    def solve(self, exposure, sourceCat):
        """!Load reference objects overlapping an exposure, match to sources and fit a WCS

        @return an lsst.pipe.base.Struct with these fields:
        - refCat  reference object catalog of objects that overlap the exposure (with some margin)
            (an lsst::afw::table::SimpleCatalog)
        - matches  list of reference object/source matches (an lsst.afw.table.ReferenceMatchVector)
        - scatterOnSky  median on-sky separation between reference objects and sources in "matches"
            (an lsst.afw.geom.Angle)
        - matchMeta  metadata needed to unpersist matches (an lsst.daf.base.PropertyList)

        @note ignores config.forceKnownWcs
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        matchMeta = createMatchMetadata(exposure, border=self.refObjLoader.config.pixelMargin)
        expMd = self._getExposureMetadata(exposure)

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            calib=expMd.calib,
        )
        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=loadRes.refCat,
                sourceCat=sourceCat,
                exposure=exposure,
                bbox=expMd.bbox,
                frame=frame,
                title="Reference catalog",
            )

        res = None
        wcs = expMd.wcs
        maxMatchDist = None
        for i in range(self.config.maxIter):
            iterNum = i + 1
            try:
                tryRes = self._matchAndFitWcs(  # refCat, sourceCat, refFluxField, bbox, wcs, exposure=None
                    refCat=loadRes.refCat,
                    sourceCat=sourceCat,
                    refFluxField=loadRes.fluxField,
                    bbox=expMd.bbox,
                    wcs=wcs,
                    exposure=exposure,
                    maxMatchDist=maxMatchDist,
                )
            except Exception as e:
                # if we have had a succeessful iteration then use that; otherwise fail
                if i > 0:
                    self.log.info("Fit WCS iter %d failed; using previous iteration: %s" % (iterNum, e))
                    iterNum -= 1
                    break
                else:
                    raise

            tryMatchDist = self._computeMatchStatsOnSky(tryRes.matches)
            self.log.debug(
                "Match and fit WCS iteration %d: found %d matches with scatter = %0.3f +- %0.3f arcsec; "
                "max match distance = %0.3f arcsec",
                iterNum, len(tryRes.matches), tryMatchDist.distMean.asArcseconds(),
                tryMatchDist.distStdDev.asArcseconds(), tryMatchDist.maxMatchDist.asArcseconds())
            if maxMatchDist is not None:
                if tryMatchDist.maxMatchDist >= maxMatchDist:
                    self.log.debug(
                        "Iteration %d had no better maxMatchDist; using previous iteration", iterNum)
                    iterNum -= 1
                    break

            maxMatchDist = tryMatchDist.maxMatchDist
            res = tryRes
            wcs = res.wcs
            if tryMatchDist.maxMatchDist.asArcseconds() < self.config.minMatchDistanceArcSec:
                self.log.debug(
                    "Max match distance = %0.3f arcsec < %0.3f = config.minMatchDistanceArcSec; "
                    "that's good enough",
                    tryMatchDist.maxMatchDist.asArcseconds(), self.config.minMatchDistanceArcSec)
                break

        self.log.info(
            "Matched and fit WCS in %d iterations; "
            "found %d matches with scatter = %0.3f +- %0.3f arcsec" %
            (iterNum, len(tryRes.matches), tryMatchDist.distMean.asArcseconds(),
                tryMatchDist.distStdDev.asArcseconds()))

        exposure.setWcs(res.wcs)

        return pipeBase.Struct(
            refCat=loadRes.refCat,
            matches=res.matches,
            scatterOnSky=res.scatterOnSky,
            matchMeta=matchMeta,
        )

    def _computeMatchStatsOnSky(self, matchList):
        """Compute on-sky radial distance statistics for a match list

        @param[in] matchList  list of matches between reference object and sources;
            the distance field is the only field read and it must be set to distance in radians

        @return a pipe_base Struct containing these fields:
        - distMean  clipped mean of on-sky radial separation
        - distStdDev  clipped standard deviation of on-sky radial separation
        - maxMatchDist  distMean + self.config.matchDistanceSigma*distStdDev
        """
        distStatsInRadians = makeMatchStatistics(matchList, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        distMean = distStatsInRadians.getValue(afwMath.MEANCLIP)*afwGeom.radians
        distStdDev = distStatsInRadians.getValue(afwMath.STDEVCLIP)*afwGeom.radians
        return pipeBase.Struct(
            distMean=distMean,
            distStdDev=distStdDev,
            maxMatchDist=distMean + self.config.matchDistanceSigma*distStdDev,
        )

    def _getExposureMetadata(self, exposure):
        """!Extract metadata from an exposure

        @return an lsst.pipe.base.Struct containing the following exposure metadata:
        - bbox: parent bounding box
        - wcs: WCS (an lsst.afw.image.Wcs)
        - calib calibration (an lsst.afw.image.Calib), or None if unknown
        - filterName: name of filter, or None if unknown
        """
        exposureInfo = exposure.getInfo()
        filterName = exposureInfo.getFilter().getName() or None
        if filterName == "_unknown_":
            filterName = None
        return pipeBase.Struct(
            bbox=exposure.getBBox(),
            wcs=getDistortedWcs(exposureInfo, log=self.log),
            calib=exposureInfo.getCalib() if exposureInfo.hasCalib() else None,
            filterName=filterName,
        )

    @pipeBase.timeMethod
    def _matchAndFitWcs(self, refCat, sourceCat, refFluxField, bbox, wcs, maxMatchDist=None,
                        exposure=None):
        """!Match sources to reference objects and fit a WCS

        @param[in] refCat  catalog of reference objects
        @param[in] sourceCat  catalog of sources detected on the exposure (an lsst.afw.table.SourceCatalog)
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] bbox  bounding box of exposure (an lsst.afw.geom.Box2I)
        @param[in] wcs  initial guess for WCS of exposure (an lsst.afw.image.Wcs)
        @param[in] maxMatchDist  maximum on-sky distance between reference objects and sources
            (an lsst.afw.geom.Angle); if None then use the matcher's default
        @param[in] exposure  exposure whose WCS is to be fit, or None; used only for the debug display

        @return an lsst.pipe.base.Struct with these fields:
        - matches  list of reference object/source matches (an lsst.afw.table.ReferenceMatchVector)
        - wcs  the fit WCS (an lsst.afw.image.Wcs)
        - scatterOnSky  median on-sky separation between reference objects and sources in "matches"
            (an lsst.afw.geom.Angle)
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)
        matchRes = self.matcher.matchObjectsToSources(
            refCat=refCat,
            sourceCat=sourceCat,
            wcs=wcs,
            refFluxField=refFluxField,
            maxMatchDist=maxMatchDist,
        )
        self.log.debug("Found %s matches", len(matchRes.matches))
        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=refCat,
                sourceCat=matchRes.usableSourceCat,
                matches=matchRes.matches,
                exposure=exposure,
                bbox=bbox,
                frame=frame + 1,
                title="Initial WCS",
            )

        self.log.debug("Fitting WCS")
        fitRes = self.wcsFitter.fitWcs(
            matches=matchRes.matches,
            initWcs=wcs,
            bbox=bbox,
            refCat=refCat,
            sourceCat=sourceCat,
            exposure=exposure,
        )
        fitWcs = fitRes.wcs
        scatterOnSky = fitRes.scatterOnSky
        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=refCat,
                sourceCat=matchRes.usableSourceCat,
                matches=matchRes.matches,
                exposure=exposure,
                bbox=bbox,
                frame=frame + 2,
                title="Fit TAN-SIP WCS",
            )

        return pipeBase.Struct(
            matches=matchRes.matches,
            wcs=fitWcs,
            scatterOnSky=scatterOnSky,
        )
