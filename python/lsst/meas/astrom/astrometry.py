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

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.geom as geom
from .ref_match import RefMatchTask, RefMatchConfig
from .fitTanSipWcs import FitTanSipWcsTask
from .display import displayAstrometry
from .fitAffineWcs import TransformedSkyWcsMaker


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
    shiftSize = pexConfig.RangeField(
        doc="Shift in pixels to impose on the data.",
        dtype=float,
        default=0,
        min=0,
        max=300,
    )
    rotsize = pexConfig.RangeField(
        doc="Angle in degrees to rotate the boresight.",
        dtype=float,
        default=0,
        min=-180,
        max=180,
    )
    affineXScale = pexConfig.RangeField(
        doc="Amount to distort unrotated x",
        dtype=float,
        default=1,
        min=1,
        max=10,
    )
    affineYScale = pexConfig.RangeField(
        doc="Amount to distort unrotated x",
        dtype=float,
        default=1,
        min=1,
        max=10,
    )
    affineXShear = pexConfig.RangeField(
        doc="Amount to shear the affine matrix",
        dtype=float,
        default=0,
        min=0,
        max=1,
    )
    affineYShear = pexConfig.RangeField(
        doc="Amount to shear the affine matrix",
        dtype=float,
        default=0,
        min=0,
        max=1,
    )

    def setDefaults(self):
        # Override the default source selector for astrometry tasks
        self.sourceFluxType = "Ap"

        self.sourceSelector.name = "matcher"
        self.sourceSelector["matcher"].sourceFluxType = self.sourceFluxType

        # Note that if the matcher is MatchOptimisticBTask, then the
        # default should be self.sourceSelector['matcher'].excludePixelFlags = False
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
        tmpCoord = wcs.pixelToSky(sourceCat[0].getCentroid())
        self.log.info("Before Ra: %.8f Dec: %.8f" %
                      (tmpCoord.getRa().asDegrees(), tmpCoord.getDec().asDegrees()))
        match_tolerance = None
        shift = wcs.getPixelScale().asArcseconds() * self.config.shiftSize
        angle = np.radians(self.config.rotsize)
        transWcsMaker = TransformedSkyWcsMaker(wcs)
        rotate = np.array([[np.cos(angle), -np.sin(angle)],
                           [np.sin(angle), np.cos(angle)]])
        scale = np.array([[self.config.affineXScale, 0],
                          [0, self.config.affineYScale]])
        shear = np.array([[1, self.config.affineXShear],
                          [self.config.affineYShear, 1]])
        matrix = np.dot(rotate, np.dot(scale, shear))
        self.log.info(f"Test shift in arcseconds: {shift}")
        self.log.info("Affine Matrix: [[%.5f, %.5f], [%.5f, %.5f]]"
                      % (matrix[0, 0], matrix[0, 1], matrix[1, 0], matrix[1, 1]))

        wcs = transWcsMaker.makeWcs(
            [np.random.uniform(-180, 180), shift],
            matrix)
        tmpCoord = wcs.pixelToSky(sourceCat[0].getCentroid())
        self.log.info("After Ra: %.8f Dec: %.8f" %
                      (tmpCoord.getRa().asDegrees(), tmpCoord.getDec().asDegrees()))
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
            "found %d matches with scatter = %0.3f +- %0.3f arcsec" %
            (iterNum, len(tryRes.matches), tryMatchDist.distMean.asArcseconds(),
                tryMatchDist.distStdDev.asArcseconds()))
        for m in res.matches:
            if self.usedKey:
                m.second.set(self.usedKey, True)
        exposure.setWcs(res.wcs)

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
