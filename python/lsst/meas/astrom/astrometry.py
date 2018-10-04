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

__all__ = ["AstrometryConfig", "AstrometryTask"]


import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .ref_match import RefMatchTask, RefMatchConfig
from .fitTanSipWcs import FitTanSipWcsTask
from .display import displayAstrometry


class AstrometryConfig(RefMatchConfig):
    wcsFitter = pexConfig.ConfigurableField(
        target=FitTanSipWcsTask,
        doc="WCS fitter",
    )
    forceKnownWcs = pexConfig.Field(
        dtype=bool,
        doc="If True then load reference objects and match sources but do not fit a WCS; "
        " this simply controls whether 'run' calls 'solve' or 'loadAndMatch'",
        default=False,
    )
    maxIter = pexConfig.RangeField(
        doc="maximum number of iterations of match sources and fit WCS"
        "ignored if not fitting a WCS",
        dtype=int,
        default=3,
        min=1,
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


class AstrometryTask(RefMatchTask):
    # Parameters
    # ----------
    # refObjLoader :
    #     A reference object loader object
    # schema :
    #     ignored; available for compatibility with an older astrometry task
    # kwargs :
    #     additional keyword arguments for pipe_base Task.\_\_init\_\_
    r"""!Match an input source catalog with objects from a reference catalog and solve for the WCS

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

    See \ref pipe_tasks_photocal_Example.

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
        RefMatchTask.__init__(self, refObjLoader, schema=schema, **kwargs)

        if schema is not None:
            self.usedKey = schema.addField("calib_astrometry_used", type="Flag",
                                           doc="set if source was used in astrometric calibration")
        else:
            self.usedKey = None

        self.makeSubtask("wcsFitter")

    @pipeBase.timeMethod
    def run(self, sourceCat, exposure):
        """Load reference objects, match sources and optionally fit a WCS

        This is a thin layer around solve or loadAndMatch, depending on config.forceKnownWcs

        Parameters
        ----------
        exposure :
            exposure whose WCS is to be fit
            The following are read only:

            - bbox
            - calib (may be absent)
            - filter (may be unset)
            - detector (if wcs is pure tangent; may be absent)
            The following are updated:
            - wcs (the initial value is used as an initial guess, and is required)

        sourceCat :lsst.afw.table.SourceCatalog
            catalog of sources detected on the exposure

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            with these fields:

        - ``refCat`` :  reference object catalog of objects that overlap the exposure (with some margin)
            (an lsst::afw::table::SimpleCatalog)
        - ``matches`` : astrometric matches, a list of lsst.afw.table.ReferenceMatch
        - ``scatterOnSky`` :  median on-sky separation between reference objects and sources in "matches"
            (an lsst.afw.geom.Angle), or None if config.forceKnownWcs True
        - ``matchMeta`` :  metadata needed to unpersist matches (an lsst.daf.base.PropertyList)
        """
        if self.config.forceKnownWcs:
            res = self.loadAndMatch(exposure=exposure, sourceCat=sourceCat)
            res.scatterOnSky = None
        else:
            res = self.solve(exposure=exposure, sourceCat=sourceCat)
        return res

    @pipeBase.timeMethod
    def solve(self, exposure, sourceCat):
        """Load reference objects overlapping an exposure, match to sources and fit a WCS

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            with these fields:

            - ``refCat`` :  reference object catalog of objects that overlap the exposure (with some margin)
                (an lsst::afw::table::SimpleCatalog)
            - ``matches`` :  astrometric matches, a list of lsst.afw.table.ReferenceMatch
            - ``scatterOnSky`` :  median on-sky separation between reference objects and sources in "matches"
                (an lsst.afw.geom.Angle)
            - ``matchMeta`` :  metadata needed to unpersist matches (an lsst.daf.base.PropertyList)

        Notes
        -----
        ignores config.forceKnownWcs
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        expMd = self._getExposureMetadata(exposure)

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            calib=expMd.calib,
            epoch=expMd.epoch,
        )
        matchMeta = self.refObjLoader.getMetadataBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            calib=expMd.calib,
            epoch=expMd.epoch,
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
        match_tolerance = None
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
                    match_tolerance=match_tolerance,
                )
            except Exception as e:
                # if we have had a succeessful iteration then use that; otherwise fail
                if i > 0:
                    self.log.info("Fit WCS iter %d failed; using previous iteration: %s" % (iterNum, e))
                    iterNum -= 1
                    break
                else:
                    raise

            match_tolerance = tryRes.match_tolerance
            tryMatchDist = self._computeMatchStatsOnSky(tryRes.matches)
            self.log.debug(
                "Match and fit WCS iteration %d: found %d matches with scatter = %0.3f +- %0.3f arcsec; "
                "max match distance = %0.3f arcsec",
                iterNum, len(tryRes.matches), tryMatchDist.distMean.asArcseconds(),
                tryMatchDist.distStdDev.asArcseconds(), tryMatchDist.maxMatchDist.asArcseconds())

            maxMatchDist = tryMatchDist.maxMatchDist
            res = tryRes
            wcs = res.wcs
            if maxMatchDist.asArcseconds() < self.config.minMatchDistanceArcSec:
                self.log.debug(
                    "Max match distance = %0.3f arcsec < %0.3f = config.minMatchDistanceArcSec; "
                    "that's good enough",
                    maxMatchDist.asArcseconds(), self.config.minMatchDistanceArcSec)
                break
            match_tolerance.maxMatchDist = maxMatchDist

        self.log.info(
            "Matched and fit WCS in %d iterations; "
            "found %d matches with scatter = %0.3f +- %0.3f arcsec" %
            (iterNum, len(tryRes.matches), tryMatchDist.distMean.asArcseconds(),
                tryMatchDist.distStdDev.asArcseconds()))
        for m in res.matches:
            if self.usedKey:
                m.second.set(self.usedKey, True)
        exposure.setWcs(res.wcs)

        return pipeBase.Struct(
            refCat=loadRes.refCat,
            matches=res.matches,
            scatterOnSky=res.scatterOnSky,
            matchMeta=matchMeta,
        )

    @pipeBase.timeMethod
    def _matchAndFitWcs(self, refCat, sourceCat, refFluxField, bbox, wcs, match_tolerance,
                        exposure=None):
        """Match sources to reference objects and fit a WCS

        Parameters
        ----------
        refCat :
            catalog of reference objects
        sourceCat :
            catalog of sources detected on the exposure (an lsst.afw.table.SourceCatalog)
        refFluxField :
            field of refCat to use for flux
        bbox :
            bounding box of exposure (an lsst.afw.geom.Box2I)
        wcs :
            initial guess for WCS of exposure (an lsst.afw.geom.Wcs)
        match_tolerance :
            a MatchTolerance object (or None) specifying
            internal tolerances to the matcher. See the MatchTolerance
            definition in the respective matcher for the class definition.
        exposure :
            exposure whose WCS is to be fit, or None; used only for the debug display

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            with these fields:

            - ``matches`` :  astrometric matches, a list of lsst.afw.table.ReferenceMatch
            - ``wcs`` :  the fit WCS (an lsst.afw.geom.Wcs)
            - ``scatterOnSky`` :  median on-sky separation between reference objects and sources in "matches"
                (an lsst.afw.geom.Angle)
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)
        matchRes = self.matcher.matchObjectsToSources(
            refCat=refCat,
            sourceCat=sourceCat,
            wcs=wcs,
            refFluxField=refFluxField,
            match_tolerance=match_tolerance,
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
            match_tolerance=matchRes.match_tolerance,
        )
