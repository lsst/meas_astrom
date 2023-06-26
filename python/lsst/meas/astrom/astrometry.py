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

import numpy as np
from astropy import units
import scipy.stats

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod
from .ref_match import RefMatchTask, RefMatchConfig
from .fitTanSipWcs import FitTanSipWcsTask
from .display import displayAstrometry


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
    maxMeanDistanceArcsec = pexConfig.RangeField(
        doc="Maximum mean on-sky distance (in arcsec) between matched source and rerference "
            "objects post-fit.  A mean distance greater than this threshold raises a TaskError "
            "and the WCS fit is considered a failure.  The default is set to the maximum tolerated "
            "by the external global calibration (e.g. jointcal) step for conceivable recovery.  "
            "Appropriate value will be dataset and workflow dependent.",
        dtype=float,
        default=0.5,
        min=0,
    )
    doMagnitudeOutlierRejection = pexConfig.Field(
        dtype=bool,
        doc=("If True then a rough zeropoint will be computed from matched sources "
             "and outliers will be rejected in the iterations."),
        default=False,
    )
    magnitudeOutlierRejectionNSigma = pexConfig.Field(
        dtype=float,
        doc=("Number of sigma (measured from the distribution) in magnitude "
             "for a potential reference/source match to be rejected during "
             "iteration."),
        default=3.0,
    )

    def setDefaults(self):
        # Override the default source selector for astrometry tasks
        self.sourceFluxType = "Ap"

        self.sourceSelector.name = "matcher"
        self.sourceSelector["matcher"].sourceFluxType = self.sourceFluxType

        # Note that if the matcher is MatchOptimisticBTask, then the default
        # should be self.sourceSelector['matcher'].excludePixelFlags = False
        # However, there is no way to do this automatically.


class AstrometryTask(RefMatchTask):
    """Match an input source catalog with objects from a reference catalog and
    solve for the WCS.

    This task is broken into two main subasks: matching and WCS fitting which
    are very interactive. The matching here can be considered in part a first
    pass WCS fitter due to the fitter's sensitivity to outliers.

    Parameters
    ----------
    refObjLoader : `lsst.meas.algorithms.ReferenceLoader`
        A reference object loader object; gen3 pipeline tasks will pass `None`
        and call `setRefObjLoader` in `runQuantum`.
    schema : `lsst.afw.table.Schema`
        Used to set "calib_astrometry_used" flag in output source catalog.
    **kwargs
        additional keyword arguments for pipe_base
        `lsst.pipe.base.Task.__init__`
    """
    ConfigClass = AstrometryConfig
    _DefaultName = "astrometricSolver"

    def __init__(self, refObjLoader=None, schema=None, **kwargs):
        RefMatchTask.__init__(self, refObjLoader=refObjLoader, **kwargs)

        if schema is not None:
            self.usedKey = schema.addField("calib_astrometry_used", type="Flag",
                                           doc="set if source was used in astrometric calibration")
        else:
            self.usedKey = None

        self.makeSubtask("wcsFitter")

    @timeMethod
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

    @timeMethod
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

        Raises
        ------
        TaskError
            If the measured mean on-sky distance between the matched source and
            reference objects is greater than
            ``self.config.maxMeanDistanceArcsec``.

        Notes
        -----
        ignores config.forceKnownWcs
        """
        if self.refObjLoader is None:
            raise RuntimeError("Running matcher task with no refObjLoader set in __init__ or setRefObjLoader")
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        expMd = self._getExposureMetadata(exposure)

        if self.config.doSourceSelection:
            sourceSelection = self.sourceSelector.run(sourceCat)
            self.log.info("Purged %d sources, leaving %d good sources",
                          len(sourceCat) - len(sourceSelection.sourceCat),
                          len(sourceSelection.sourceCat))
            catalog = sourceSelection.sourceCat
        else:
            catalog = sourceCat

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            epoch=expMd.epoch,
        )

        refSelection = self.referenceSelector.run(loadRes.refCat)

        matchMeta = self.refObjLoader.getMetadataBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            epoch=expMd.epoch,
        )

        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=refSelection.sourceCat,
                sourceCat=catalog,
                exposure=exposure,
                bbox=expMd.bbox,
                frame=frame,
                title="Reference catalog",
            )

        res = None
        wcs = expMd.wcs
        match_tolerance = None
        fitFailed = False
        for i in range(self.config.maxIter):
            if not fitFailed:
                iterNum = i + 1
                try:
                    tryRes = self._matchAndFitWcs(
                        refCat=refSelection.sourceCat,
                        sourceCat=sourceCat,
                        goodSourceCat=catalog,
                        refFluxField=loadRes.fluxField,
                        bbox=expMd.bbox,
                        wcs=wcs,
                        exposure=exposure,
                        match_tolerance=match_tolerance,
                    )
                except Exception as e:
                    # If we have had a succeessful iteration then use that;
                    # otherwise fail.
                    if i > 0:
                        self.log.info("Fit WCS iter %d failed; using previous iteration: %s", iterNum, e)
                        iterNum -= 1
                        break
                    else:
                        self.log.info("Fit WCS iter %d failed: %s" % (iterNum, e))
                        fitFailed = True

            if not fitFailed:
                match_tolerance = tryRes.match_tolerance
                tryMatchDist = self._computeMatchStatsOnSky(tryRes.matches)
                self.log.debug(
                    "Match and fit WCS iteration %d: found %d matches with on-sky distance mean and "
                    "scatter = %0.3f +- %0.3f arcsec; max match distance = %0.3f arcsec",
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

        if not fitFailed:
            self.log.info("Matched and fit WCS in %d iterations; "
                          "found %d matches with mean and scatter = %0.3f +- %0.3f arcsec" %
                          (iterNum, len(tryRes.matches), tryMatchDist.distMean.asArcseconds(),
                           tryMatchDist.distStdDev.asArcseconds()))
            if tryMatchDist.distMean.asArcseconds() > self.config.maxMeanDistanceArcsec:
                self.log.info("Assigning as a fit failure: mean on-sky distance = %0.3f arcsec > %0.3f "
                              "(maxMeanDistanceArcsec)" % (tryMatchDist.distMean.asArcseconds(),
                                                           self.config.maxMeanDistanceArcsec))
                fitFailed = True

        if fitFailed:
            self.log.warning("WCS fit failed.  Setting exposure's WCS to None and coord_ra & coord_dec "
                             "cols in sourceCat to nan.")
            sourceCat["coord_ra"] = np.nan
            sourceCat["coord_dec"] = np.nan
            exposure.setWcs(None)
            matches = None
            scatterOnSky = None
        else:
            for m in res.matches:
                if self.usedKey:
                    m.second.set(self.usedKey, True)
            exposure.setWcs(res.wcs)
            matches = res.matches
            scatterOnSky = res.scatterOnSky

        # If fitter converged, record the scatter in the exposure metadata
        # even if the fit was deemed a failure according to the value of
        # the maxMeanDistanceArcsec config.
        if res is not None:
            md = exposure.getMetadata()
            md['SFM_ASTROM_OFFSET_MEAN'] = tryMatchDist.distMean.asArcseconds()
            md['SFM_ASTROM_OFFSET_STD'] = tryMatchDist.distStdDev.asArcseconds()

        return pipeBase.Struct(
            refCat=refSelection.sourceCat,
            matches=matches,
            scatterOnSky=scatterOnSky,
            matchMeta=matchMeta,
        )

    @timeMethod
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

        if self.config.doMagnitudeOutlierRejection:
            matches = self._removeMagnitudeOutliers(sourceFluxField, refFluxField, matchRes.matches)
        else:
            matches = matchRes.matches

        self.log.debug("Fitting WCS")
        fitRes = self.wcsFitter.fitWcs(
            matches=matches,
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
                matches=matches,
                exposure=exposure,
                bbox=bbox,
                frame=frame + 2,
                title="Fit TAN-SIP WCS",
            )

        return pipeBase.Struct(
            matches=matches,
            wcs=fitWcs,
            scatterOnSky=scatterOnSky,
            match_tolerance=matchRes.match_tolerance,
        )

    def _removeMagnitudeOutliers(self, sourceFluxField, refFluxField, matchesIn):
        """Remove magnitude outliers, computing a simple zeropoint.

        Parameters
        ----------
        sourceFluxField : `str`
            Field in source catalog for instrumental fluxes.
        refFluxField : `str`
            Field in reference catalog for fluxes (nJy).
        matchesIn : `list` [`lsst.afw.table.ReferenceMatch`]
            List of source/reference matches input

        Returns
        -------
        matchesOut : `list` [`lsst.afw.table.ReferenceMatch`]
            List of source/reference matches with magnitude
            outliers removed.
        """
        nMatch = len(matchesIn)
        sourceMag = np.zeros(nMatch)
        refMag = np.zeros(nMatch)
        for i, match in enumerate(matchesIn):
            sourceMag[i] = -2.5*np.log10(match[1][sourceFluxField])
            refMag[i] = (match[0][refFluxField]*units.nJy).to_value(units.ABmag)

        deltaMag = refMag - sourceMag
        # Protect against negative fluxes and nans in the reference catalog.
        goodDelta, = np.where(np.isfinite(deltaMag))
        zp = np.median(deltaMag[goodDelta])
        # Use median absolute deviation (MAD) for zpSigma.
        # Also require a minimum scatter to prevent floating-point errors from
        # rejecting objects in zero-noise tests.
        zpSigma = np.clip(scipy.stats.median_abs_deviation(deltaMag[goodDelta], scale='normal'),
                          1e-3,
                          None)

        self.log.info("Rough zeropoint from astrometry matches is %.4f +/- %.4f.",
                      zp, zpSigma)

        goodStars = goodDelta[(np.abs(deltaMag[goodDelta] - zp)
                               <= self.config.magnitudeOutlierRejectionNSigma*zpSigma)]

        nOutlier = nMatch - goodStars.size
        self.log.info("Removed %d magnitude outliers out of %d total astrometry matches.",
                      nOutlier, nMatch)

        matchesOut = []
        for matchInd in goodStars:
            matchesOut.append(matchesIn[matchInd])

        return matchesOut
