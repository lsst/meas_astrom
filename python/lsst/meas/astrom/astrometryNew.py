# Sketch of new AstrometryTask using new pipeline exceptions

import numpy as np
import astropy.time

from lsst.daf.base import DateTime
import lsst.pipe.base as pipeBase


class AstrometryError(pipeBase.RepeatableQuantumError):
    """Parent class for failures in astrometric fitting."""
    def __init__(self, *, msg, nMatches, iterations):
        super.__init__()
        self.nMatches = nMatches
        self.iterations = iterations

    @property
    def metadata(self):
        return {"nMatches": self.nMatches,
                "iterations": self.iterations
                }


class AstrometryFailure(AstrometryError):
    """Raised if the astrometry fitter fails."""
    def __init__(self, *args):
        super.__init__(msg="Failed to fit astrometry.", *args)


class BadAstrometryFit(AstrometryError):
    """Raised if the quality of the astrometric fit is worse than some
    threshold.
    """
    def __init__(self, *args, quality, **kwargs):
        msg = "Poor quality astrometric fit."
        super.__init__(msg=msg, *args, **kwargs)
        self.quality = quality

    @property
    def metadata(self):
        temp = super.metadata()
        temp["quality"] = self.quality
        return temp


class AstrometryTask(pipeBase.Task):
    def solve(self, exposure, sourceCat):
        sourceSelection = self.sourceSelector.run(sourceCat)

        epoch = astropy.time.Time(exposure.visitInfo.date.get(system=DateTime.MJD, scale=DateTime.TAI),
                                  scale="tai", format="mjd")
        loadResult = self.refObjLoader.loadPixelBox(
            bbox=exposure.getBBox(),
            wcs=exposure.wcs,
            filterName=exposure.filter.bandLabel,
            epoch=epoch,
        )
        refSelection = self.referenceSelector.run(loadResult.refCat)
        matchMeta = self.refObjLoader.getMetadataBox(
            bbox=exposure.getBBox(),
            wcs=exposure.wcs,
            filterName=exposure.filterLabel.band,
            epoch=epoch,
        )

        matchTolerance = None
        i = 0
        maxMatchDistance = np.inf
        while (maxMatchDistance > self.config.minMatchDistanceArcSec and i < self.config.maxIter):
            i += 1
            result = self._fit(i,
                               refSelection.sourceCat,
                               sourceSelection.sourceCat,
                               loadResult.fluxField,
                               exposure,
                               matchTolerance)
            matchTolerance = result.match_tolerance
            maxMatchDistance = result.stats.maxMatchDist.asArcseconds()

        # Poor quality fits are a failure.
        if result.stats.distMean.asArcseconds() > self.config.maxMeanDistanceArcsec:
            exception = BadAstrometryFit(nMatches=len(result.matches), iterations=i)
            self._fail(exception, exposure)

        self.log.info("Matched and fit WCS in %d iterations; "
                      "found %d matches with mean and scatter = %0.3f +- %0.3f arcsec",
                      i, len(result.matches), result.stats.distMean.asArcseconds(),
                      result.stats.distStdDev.asArcseconds())

        return pipeBase.Struct(
            refCat=refSelection.sourceCat,
            matches=result.matches,
            scatterOnSky=result.scatterOnSky,
            matchMeta=matchMeta,
        )

    def _fit(self, i, refs, sources, fluxField, exposure, match_tolerance):
        result = self._matchAndFitWcs(
            refCat=refs,
            sourceCat=sources,
            refFluxField=fluxField,
            bbox=exposure.getBBox(),
            wcs=exposure.wcs,
            exposure=exposure,
            match_tolerance=match_tolerance,
        )

        stats = self._computeMatchStatsOnSky(result.matches)
        self.log.debug("Match and fit WCS iteration %d: found %d matches with on-sky distance mean and "
                       "scatter = %0.3f +- %0.3f arcsec; max match distance = %0.3f arcsec",
                       i, len(result.matches), stats.distMean.asArcseconds(),
                       stats.distStdDev.asArcseconds(), stats.maxMatchDist.asArcseconds())

        result.stats = stats
        return result

    def _fail(self, exception, exposure):
        """Emit appropriate messages and clear the exposure WCS due to a
        failure, and raise.

        Parameters
        ----------
        exception : `AstrometryError`
            Description
        exposure : `lsst.afw.image.Exposure`
            Description

        Raises
        ------
        AstrometryError
            Passed in ``exception`` is raised after messages and cleanup are
            complete.
        """
        self.log.warning(exception.msg)
        exposure.setWcs(None)
        raise exception
