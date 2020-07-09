#
# LSST Data Management System
# See COPYRIGHT file at the top of the source tree.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program. If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ["AstrometryConfig", "AstrometryTask"]

import copy

import lsst.geom
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .ref_match import RefMatchTask, RefMatchConfig
from .fitTanSipWcs import FitTanSipWcsTask
from .display import displayAstrometry
from . import makeMatchStatistics
from .directMatch import DirectMatchTask


class AstrometryConfig(RefMatchConfig):
    """Config for AstrometryTask.
    """
    wcsFitter = pexConfig.ConfigurableField(
        target=FitTanSipWcsTask,
        doc="WCS fitter",
    )
    forceKnownWcs = pexConfig.Field(
        dtype=bool,
        doc="If True then load reference objects and match sources but do not fit a WCS; "
            "this simply controls whether 'run' calls 'solve' or 'loadAndMatch'",
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
    doCheckAfterRematch = pexConfig.Field(
        doc="do double-check WCS fit scatter, rematching to refCat",
        dtype=bool,
        default=False,
    )
    maxScatterArcsec = pexConfig.Field(
        doc="maximum permitted rematched scatter",
        dtype=float,
        default=10,
    )
    rematcher = pexConfig.ConfigurableField(
        target=DirectMatchTask,
        doc="simple matcher to confirm fit if doCheckAfterRematch",
    )

    def setDefaults(self):
        # Override the default source selector for astrometry tasks
        self.sourceFluxType = "Ap"

        self.sourceSelector.name = "matcher"
        self.sourceSelector["matcher"].sourceFluxType = self.sourceFluxType

        # Note that if the matcher is MatchOptimisticBTask, then the
        # default should be self.sourceSelector['matcher'].excludePixelFlags = False
        # However, there is no way to do this automatically.

        self.rematcher.referenceSelection = self.referenceSelector
        self.rematcher.sourceSelection.fluxLimit.fluxField = \
            'slot_%sFlux_instFlux' % (self.sourceFluxType)
        self.rematcher.sourceSelection.signalToNoise.fluxField = \
            'slot_%sFlux_instFlux' % (self.sourceFluxType)
        self.rematcher.sourceSelection.signalToNoise.errField = \
            'slot_%sFlux_instFluxErr' % (self.sourceFluxType)
        # defaults below chosen to mimic Dominique Boutigny's checkCcdAstrometry script
        # extra flags to consider for the future: 'base_PixelFlags_flag_cr'
        #                                         'base_SdssCentroid_flag'
        self.rematcher.sourceSelection.doFluxLimit = True
        self.rematcher.sourceSelection.doSignalToNoise = True
        self.rematcher.sourceSelection.fluxLimit.minimum = 0
        self.rematcher.sourceSelection.signalToNoise.minimum = 5
        self.rematcher.sourceSelection.doFlags = True
        sourceBadFlags = ['base_PixelFlags_flag_edge',
                          'base_PixelFlags_flag_saturated',
                          'base_PixelFlags_flag_interpolated']
        self.rematcher.sourceSelection.flags.bad = sourceBadFlags


class AstrometryTask(RefMatchTask):
    """Match an input source catalog with objects from a reference catalog and
    solve for the WCS.

    This task is broken into two main subasks: matching and WCS fitting which
    are very interactive. The matching here can be considered in part a first
    pass WCS fitter due to the fitter's sensitivity to outliers.

    Parameters
    ----------
    refObjLoader : `lsst.meas.algorithms.ReferenceLoader`
        A reference object loader object
    schema : `lsst.afw.table.Schema`
        Used to set "calib_astrometry_used" flag in output source catalog.
    **kwargs
        additional keyword arguments for pipe_base
        `lsst.pipe.base.Task.__init__`
    """
    ConfigClass = AstrometryConfig
    _DefaultName = "astrometricSolver"

    def __init__(self, refObjLoader, schema=None, **kwargs):
        RefMatchTask.__init__(self, refObjLoader, **kwargs)

        if schema is not None:
            self.usedKey = schema.addField("calib_astrometry_used", type="Flag",
                                           doc="set if source was used in astrometric calibration")
        else:
            self.usedKey = None

        self.makeSubtask("wcsFitter")

        if self.config.doCheckAfterRematch:
            self.makeSubtask("rematcher", refObjLoader=refObjLoader)

    @pipeBase.timeMethod
    def run(self, sourceCat, exposure):
        """Load reference objects, match sources and optionally fit a WCS.

        This is a thin layer around solve or loadAndMatch, depending on
        config.forceKnownWcs.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            exposure whose WCS is to be fit
            The following are read only:

            - bbox
            - photoCalib (may be absent)
            - filter (may be unset)
            - detector (if wcs is pure tangent; may be absent)

            The following are updated:

            - wcs (the initial value is used as an initial guess, and is
              required)

        sourceCat : `lsst.afw.table.SourceCatalog`
            catalog of sources detected on the exposure

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            with these fields:

            - ``refCat`` : reference object catalog of objects that overlap the
              exposure (with some margin) (`lsst.afw.table.SimpleCatalog`).
            - ``matches`` : astrometric matches
              (`list` of `lsst.afw.table.ReferenceMatch`).
            - ``scatterOnSky`` :  median on-sky separation between reference
              objects and sources in "matches"
              (`lsst.afw.geom.Angle`) or `None` if config.forceKnownWcs True
            - ``matchMeta`` :  metadata needed to unpersist matches
              (`lsst.daf.base.PropertyList`)
        """
        if self.refObjLoader is None:
            raise RuntimeError("Running matcher task with no refObjLoader set in __init__ or setRefObjLoader")
        if self.config.forceKnownWcs:
            res = self.loadAndMatch(exposure=exposure, sourceCat=sourceCat)
            res.scatterOnSky = None
        else:
            res = self.solve(exposure=exposure, sourceCat=sourceCat)
        return res

    def setRefObjLoader(self, refObjLoader):
        # Using copy instead of deepcopy to avoid a recursion error in copying
        # the config in the refObjLoader. The underlying data in refObjLoader
        # should not be modified, so this should be safe.
        # Future people: if this task is experiencing errors, look here first!
        copyRefObjLoader = copy.copy(refObjLoader)
        super().setRefObjLoader(refObjLoader)
        if self.config.doCheckAfterRematch:
            self.rematcher.setRefObjLoader(copyRefObjLoader)

    @pipeBase.timeMethod
    def solve(self, exposure, sourceCat):
        """Load reference objects overlapping an exposure, match to sources and
        fit a WCS

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``refCat`` : reference object catalog of objects that overlap the
              exposure (with some margin) (`lsst::afw::table::SimpleCatalog`).
            - ``matches`` :  astrometric matches
              (`list` of `lsst.afw.table.ReferenceMatch`).
            - ``scatterOnSky`` :  median on-sky separation between reference
              objects and sources in "matches" (`lsst.geom.Angle`)
            - ``matchMeta`` :  metadata needed to unpersist matches
              (`lsst.daf.base.PropertyList`)

        Notes
        -----
        ignores config.forceKnownWcs
        """
        if self.refObjLoader is None:
            raise RuntimeError("Running matcher task with no refObjLoader set in __init__ or setRefObjLoader")
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        expMd = self._getExposureMetadata(exposure)

        sourceSelection = self.sourceSelector.run(sourceCat)

        self.log.info("Purged %d sources, leaving %d good sources" %
                      (len(sourceCat) - len(sourceSelection.sourceCat),
                       len(sourceSelection.sourceCat)))

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            photoCalib=expMd.photoCalib,
            epoch=expMd.epoch,
        )

        refSelection = self.referenceSelector.run(loadRes.refCat)

        matchMeta = self.refObjLoader.getMetadataBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            photoCalib=expMd.photoCalib,
            epoch=expMd.epoch,
        )

        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=refSelection.sourceCat,
                sourceCat=sourceSelection.sourceCat,
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
                tryRes = self._matchAndFitWcs(
                    refCat=refSelection.sourceCat,
                    sourceCat=sourceCat,
                    goodSourceCat=sourceSelection.sourceCat,
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
            "found %d matches with clipped-mean scatter = %0.3f +- %0.3f arcsec" %
            (iterNum, len(tryRes.matches), tryMatchDist.distMean.asArcseconds(),
                tryMatchDist.distStdDev.asArcseconds()))
        for m in res.matches:
            if self.usedKey:
                m.second.set(self.usedKey, True)
        exposure.setWcs(res.wcs)

        # perform a rematching check
        if self.config.doCheckAfterRematch:

            sourceCatCopy = sourceSelection.sourceCat.copy(deep=True)
            rematchRes = self.rematcher.run(sourceCatCopy, exposure.getInfo().getFilter().getName())

            if not rematchRes.matches:
                errMsg = (f'WCS fit failed, no matches to refCat within a match tolerance of '
                          f'{self.rematcher.config.matchRadius} arcsec')
                raise pipeBase.TaskError(errMsg)

            distStats = makeMatchStatistics(rematchRes.matches, afwMath.MEANCLIP | afwMath.STDEVCLIP)
            distMean = distStats.getValue(afwMath.MEANCLIP)*lsst.geom.radians
            distStdDev = distStats.getValue(afwMath.STDEVCLIP)*lsst.geom.radians
            self.log.info(
                "Rematch check found %d matches with clipped-mean scatter = %0.3f +- %0.3f arcsec",
                len(rematchRes.matches),
                distMean.asArcseconds(),
                distStdDev.asArcseconds(),
            )

            if distMean.asArcseconds() > self.config.maxScatterArcsec:
                errMsg = (f'WCS fit failed, rematched scatter greater than '
                          f'{self.config.maxScatterArcsec} arcsec')
                raise pipeBase.TaskError(errMsg)

        return pipeBase.Struct(
            refCat=refSelection.sourceCat,
            matches=res.matches,
            scatterOnSky=res.scatterOnSky,
            matchMeta=matchMeta,
        )

    @pipeBase.timeMethod
    def _matchAndFitWcs(self, refCat, sourceCat, goodSourceCat, refFluxField, bbox, wcs, match_tolerance,
                        exposure=None):
        """Match sources to reference objects and fit a WCS.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            catalog of reference objects
        sourceCat : `lsst.afw.table.SourceCatalog`
            catalog of sources detected on the exposure
        goodSourceCat : `lsst.afw.table.SourceCatalog`
            catalog of down-selected good sources detected on the exposure
        refFluxField : 'str'
            field of refCat to use for flux
        bbox : `lsst.geom.Box2I`
            bounding box of exposure
        wcs : `lsst.afw.geom.SkyWcs`
            initial guess for WCS of exposure
        match_tolerance : `lsst.meas.astrom.MatchTolerance`
            a MatchTolerance object (or None) specifying
            internal tolerances to the matcher. See the MatchTolerance
            definition in the respective matcher for the class definition.
        exposure : `lsst.afw.image.Exposure`
            exposure whose WCS is to be fit, or None; used only for the debug
            display.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``matches``:  astrometric matches
              (`list` of `lsst.afw.table.ReferenceMatch`).
            - ``wcs``:  the fit WCS (lsst.afw.geom.SkyWcs).
            - ``scatterOnSky`` :  median on-sky separation between reference
              objects and sources in "matches" (`lsst.afw.geom.Angle`).
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        sourceFluxField = "slot_%sFlux_instFlux" % (self.config.sourceFluxType)

        matchRes = self.matcher.matchObjectsToSources(
            refCat=refCat,
            sourceCat=goodSourceCat,
            wcs=wcs,
            sourceFluxField=sourceFluxField,
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
