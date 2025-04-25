# This file is part of meas_astrom.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["AstrometryConfig", "AstrometryTask"]

import numpy as np
from astropy import units
import scipy.stats

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod
from . import exceptions
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
        doc="Maximum mean on-sky distance (in arcsec) between matched source and reference "
            "objects post-fit. A mean distance greater than this threshold raises BadAstrometryFit "
            "and the WCS fit is considered a failure. The default is set to the maximum tolerated "
            "by the external global calibration (e.g. jointcal) step for conceivable recovery; "
            "the appropriate value will be dataset and workflow dependent.",
        dtype=float,
        default=0.5,
        min=0,
    )
    doMagnitudeOutlierRejection = pexConfig.Field(
        dtype=bool,
        doc=("If True then a rough zeropoint will be computed from matched sources "
             "and outliers will be rejected in the iterations."),
        default=True,
    )
    magnitudeOutlierRejectionNSigma = pexConfig.Field(
        dtype=float,
        doc=("Number of sigma (measured from the distribution) in magnitude "
             "for a potential reference/source match to be rejected during "
             "iteration."),
        default=3.0,
    )
    fiducialZeroPoint = pexConfig.DictField(
        keytype=str,
        itemtype=float,
        doc="Fiducial zeropoint values, keyed by band.",
        default=None,
        optional=True,
    )
    doFiducialZeroPointCull = pexConfig.Field(
        dtype=bool,
        doc="If True, use the obs_package-defined fiducial zeropoint values to compute the "
        "expected zeropoint for the current exposure.  This is for use in culling reference "
        "objects down to the approximate magnitude range of the detection source catalog "
        "used for astrometric calibration.",
        default=False,
    )
    cullMagBuffer = pexConfig.Field(
        dtype=float,
        doc="Generous buffer on the fiducial zero point culling to account for observing "
        "condition variations and uncertainty of the fiducial values.",
        default=0.3,
        optional=True,
    )
    doRandomDownsample = pexConfig.Field(
        dtype=bool,
        doc="Do random downsampling of reference catalog objects to send to the "
            "matcher task? This can only be used if doFiducialZeroPointCull=True "
            "to ensure that the random downsampling is using an appropriate "
            "magnitude range.",
        default=False,
    )
    randomDownsampleMaxObjects = pexConfig.Field(
        dtype=int,
        doc="Maximum number of objects to use in random downsampling.",
        default=2048,
    )

    def setDefaults(self):
        super().setDefaults()
        # Astrometry needs sources to be isolated, too.
        self.sourceSelector["science"].doRequirePrimary = True
        self.sourceSelector["science"].doIsolated = True
        self.referenceSelector.doCullFromMaskedRegion = True

    def validate(self):
        super().validate()
        if self.doFiducialZeroPointCull and self.fiducialZeroPoint is None:
            msg = "doFiducialZeroPointCull=True requires `fiducialZeroPoint`, a dict of the "
            "fiducial zeropoints measured for the camera/filter, be set."
            raise pexConfig.FieldValidationError(AstrometryConfig.fiducialZeroPoint, self, msg)
        if self.doRandomDownsample and not self.doFiducialZeroPointCull:
            msg = "doRandomDownsample=True requires doFiducialZeroPointCull=True."
            raise pexConfig.FieldValidationError(AstrometryConfig.fiducialZeroPoint, self, msg)


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
        Additional keyword arguments for pipe_base
        `lsst.pipe.base.Task.__init__`.
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
        BadAstrometryFit
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

        epoch = exposure.visitInfo.date.toAstropy()

        sourceSelection = self.sourceSelector.run(sourceCat)

        self.log.info("Purged %d sources, leaving %d good sources",
                      len(sourceCat) - len(sourceSelection.sourceCat),
                      len(sourceSelection.sourceCat))
        if len(sourceSelection.sourceCat) == 0:
            raise exceptions.AstrometryError(
                "No good sources selected for astrometry.",
                lenSourceSelectionCat=len(sourceSelection.sourceCat)
            )

        loadResult = self.refObjLoader.loadPixelBox(
            bbox=exposure.getBBox(),
            wcs=exposure.wcs,
            filterName=exposure.filter.bandLabel,
            epoch=epoch,
        )
        refSelection = self.referenceSelector.run(loadResult.refCat, exposure=exposure)
        refCat = refSelection.sourceCat

        if self.config.doFiducialZeroPointCull:
            nRefCatPreCull = len(refCat)
            # Compute rough limiting magnitudes of selected sources
            psfFlux = sourceSelection.sourceCat["base_PsfFlux_instFlux"]
            expTime = exposure.visitInfo.getExposureTime()
            fiducialZeroPoint = self.config.fiducialZeroPoint[exposure.filter.bandLabel]
            psfMag = -2.5*np.log10(psfFlux) + fiducialZeroPoint + 2.5*np.log10(expTime)
            sourceMagMin = min(psfMag) - self.config.cullMagBuffer
            sourceMagMax = max(psfMag) + self.config.cullMagBuffer
            self.log.info("sourceMagMin = %0.3f sourceMagMax = %0.3f", sourceMagMin, sourceMagMax)

            refCat = refCat[np.isfinite(refCat[loadResult.fluxField])]
            refMag = (refCat[loadResult.fluxField]*units.nJy).to_value(units.ABmag)
            refCat = refCat[(refMag > sourceMagMin) & (refMag < sourceMagMax)]
            self.log.info("Selected %d/%d reference sources based on fiducial zeropoint culling.",
                          len(refCat), nRefCatPreCull)

            if self.config.doRandomDownsample and len(refCat) > self.config.randomDownsampleMaxObjects:
                if not refCat.isContiguous():
                    refCat = refCat.copy(deep=True)

                if exposure.info.id is None:
                    # Fallback seed in case something is missing.
                    seed = 12345
                else:
                    seed = exposure.info.id & 0xFFFFFFFF
                rng = np.random.RandomState(seed=seed)

                sel = np.zeros(len(refCat), dtype=np.bool_)
                inds = rng.choice(len(refCat), size=self.config.randomDownsampleMaxObjects, replace=False)
                sel[inds] = True
                refCat = refCat[sel]

        # Some operations below require catalog contiguity, which is not
        # guaranteed from the source selector/culling/downsampling.
        if not refCat.isContiguous():
            refCat = refCat.copy(deep=True)

        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=refCat,
                sourceCat=sourceSelection.sourceCat,
                exposure=exposure,
                bbox=exposure.getBBox(),
                frame=frame,
                title="Reference catalog",
            )

        result = pipeBase.Struct(matchTolerance=None)
        maxMatchDistance = np.inf
        i = 0
        while (maxMatchDistance > self.config.minMatchDistanceArcSec and i < self.config.maxIter):
            i += 1
            try:
                result = self._matchAndFitWcs(
                    refCat=refCat,
                    sourceCat=sourceCat,
                    goodSourceCat=sourceSelection.sourceCat,
                    refFluxField=loadResult.fluxField,
                    bbox=exposure.getBBox(),
                    wcs=exposure.wcs,
                    exposure=exposure,
                    matchTolerance=result.matchTolerance,
                )
                exposure.setWcs(result.wcs)
            except exceptions.AstrometryError as e:
                e._metadata['iterations'] = i
                sourceCat["coord_ra"] = np.nan
                sourceCat["coord_dec"] = np.nan
                exposure.setWcs(None)
                self.log.error("Failure fitting astrometry. %s: %s", type(e).__name__, e)
                raise

            result.stats = self._computeMatchStatsOnSky(result.matches)
            maxMatchDistance = result.stats.maxMatchDist.asArcseconds()
            distMean = result.stats.distMean.asArcseconds()
            distStdDev = result.stats.distStdDev.asArcseconds()
            self.log.info("Astrometric fit iteration %d: found %d matches with mean separation "
                          "= %0.3f +- %0.3f arcsec; max match distance = %0.3f arcsec.",
                          i, len(result.matches), distMean, distStdDev, maxMatchDistance)

        # If fitter converged, record the scatter in the exposure metadata
        # even if the fit was deemed a failure according to the value of
        # the maxMeanDistanceArcsec config.
        md = exposure.getMetadata()
        md['SFM_ASTROM_OFFSET_MEAN'] = distMean
        md['SFM_ASTROM_OFFSET_STD'] = distStdDev
        md['SFM_ASTROM_OFFSET_MEDIAN'] = result.scatterOnSky.asArcseconds()

        if self.usedKey:
            for m in result.matches:
                m.second.set(self.usedKey, True)

        matchMeta = self.refObjLoader.getMetadataBox(
            bbox=exposure.getBBox(),
            wcs=exposure.wcs,
            filterName=exposure.filter.bandLabel,
            epoch=epoch,
        )

        return pipeBase.Struct(
            refCat=refCat,
            matches=result.matches,
            scatterOnSky=result.scatterOnSky,
            matchMeta=matchMeta,
        )

    def check(self, exposure, sourceCat, nMatches):
        """Validate the astrometric fit against the maxMeanDistance threshold.

        If the distMean metric does not satisfy the requirement of being less
        than the value set in config.maxMeanDistanceArcsec, the WCS on the
        exposure will be set to None and the coordinate values in the
        source catalog will be set to NaN to reflect a failed astrometric fit.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure whose astrometric fit is being evaluated.
        sourceCat : `lsst.afw.table.SourceCatalog`
            The catalog of sources associated with the exposure.
        nMatches : `int`
            The number of matches that were found and used during
            the astrometric fit (for logging purposes only).

        Raises
        ------
        BadAstrometryFit
            If the measured mean on-sky distance between the matched source and
            reference objects is greater than
            ``self.config.maxMeanDistanceArcsec``.
        """
        # Poor quality fits are a failure.
        md = exposure.getMetadata()
        distMean = md['SFM_ASTROM_OFFSET_MEAN']
        distMedian = md['SFM_ASTROM_OFFSET_MEDIAN']
        if distMean > self.config.maxMeanDistanceArcsec:
            exception = exceptions.BadAstrometryFit(nMatches=nMatches, distMean=distMean,
                                                    maxMeanDist=self.config.maxMeanDistanceArcsec,
                                                    distMedian=distMedian)
            exposure.setWcs(None)
            sourceCat["coord_ra"] = np.nan
            sourceCat["coord_dec"] = np.nan
            self.log.error(exception)
            raise exception
        return

    @timeMethod
    def _matchAndFitWcs(self, refCat, sourceCat, goodSourceCat, refFluxField, bbox, wcs, matchTolerance,
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
        matchTolerance : `lsst.meas.astrom.MatchTolerance`
            a MatchTolerance object (or None) specifying
            internal tolerances to the matcher. See the MatchTolerance
            definition in the respective matcher for the class definition.
        exposure : `lsst.afw.image.Exposure`, optional
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
            matchTolerance=matchTolerance,
            bbox=bbox,
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
                title=f"Fitter: {self.wcsFitter._DefaultName}",
            )

        return pipeBase.Struct(
            matches=matches,
            wcs=fitWcs,
            scatterOnSky=scatterOnSky,
            matchTolerance=matchRes.matchTolerance,
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

        matchesOut = [matchesIn[idx] for idx in goodStars]

        return matchesOut
