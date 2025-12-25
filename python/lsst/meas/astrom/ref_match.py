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

__all__ = ['RefMatchConfig', 'RefMatchTask']

import lsst.geom
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import ReferenceSourceSelectorTask
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry
from lsst.utils.timer import timeMethod
from .matchPessimisticB import MatchPessimisticBTask
from .display import displayAstrometry
from . import makeMatchStatistics


class RefMatchConfig(pexConfig.Config):
    matcher = pexConfig.ConfigurableField(
        target=MatchPessimisticBTask,
        doc="reference object/source matcher",
    )
    matchDistanceSigma = pexConfig.RangeField(
        doc="the maximum match distance is set to "
        " mean_match_distance + matchDistanceSigma*std_dev_match_distance; "
        "ignored if not fitting a WCS",
        dtype=float,
        default=2,
        min=0,
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources for cross-matching.",
        default="science",
    )
    referenceSelector = pexConfig.ConfigurableField(
        target=ReferenceSourceSelectorTask,
        doc="How to select reference objects for cross-matching."
    )
    sourceFluxType = pexConfig.Field(
        dtype=str,
        doc="Source flux type to use in source selection.",
        default='Psf'
    )

    def setDefaults(self):
        super().setDefaults()
        # Configured to match the deprecated "matcher" selector:
        # SN > 40, some bad flags, valid centroids.
        self.sourceSelector["science"].doSignalToNoise = True
        self.sourceSelector["science"].signalToNoise.minimum = 40
        self.sourceSelector["science"].signalToNoise.fluxField = f"slot_{self.sourceFluxType}Flux_instFlux"
        self.sourceSelector["science"].signalToNoise.errField = f"slot_{self.sourceFluxType}Flux_instFluxErr"
        self.sourceSelector["science"].doFlags = True
        self.sourceSelector["science"].flags.bad = ["base_PixelFlags_flag_edge",
                                                    "base_PixelFlags_flag_nodata",
                                                    "base_PixelFlags_flag_interpolatedCenter",
                                                    "base_PixelFlags_flag_saturated",
                                                    "base_SdssCentroid_flag",
                                                    ]


class RefMatchTask(pipeBase.Task):
    """Match an input source catalog with objects from a reference catalog.

    Parameters
    ----------
    refObjLoader : `lsst.meas.algorithms.ReferenceLoader`
        A reference object loader object; gen3 pipeline tasks will pass `None`
        and call `setRefObjLoader` in `runQuantum`.
    **kwargs
        Additional keyword arguments for pipe_base `lsst.pipe.base.Task`.
    """
    ConfigClass = RefMatchConfig
    _DefaultName = "calibrationBaseClass"

    def __init__(self, refObjLoader=None, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        if refObjLoader:
            self.refObjLoader = refObjLoader
        else:
            self.refObjLoader = None

        if self.config.sourceSelector.name == 'matcher':
            if self.config.sourceSelector['matcher'].sourceFluxType != self.config.sourceFluxType:
                raise RuntimeError("The sourceFluxType in the sourceSelector['matcher'] must match "
                                   "the configured sourceFluxType")

        self.makeSubtask("matcher")
        self.makeSubtask("sourceSelector")
        self.makeSubtask("referenceSelector")

    def setRefObjLoader(self, refObjLoader):
        """Sets the reference object loader for the task.

        Parameters
        ----------
        refObjLoader
            An instance of a reference object loader task or class.
        """
        self.refObjLoader = refObjLoader

    @timeMethod
    def loadAndMatch(self, exposure, sourceCat):
        """Load reference objects overlapping an exposure and match to sources
        detected on that exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            exposure that the sources overlap
        sourceCat : `lsst.afw.table.SourceCatalog.`
            catalog of sources detected on the exposure

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with Components:

            - ``refCat`` : reference object catalog of objects that overlap the
              exposure (`lsst.afw.table.SimpleCatalog`)
            - ``matches`` :  Matched sources and references
              (`list` of `lsst.afw.table.ReferenceMatch`)
            - ``matchMeta`` : metadata needed to unpersist matches
              (`lsst.daf.base.PropertyList`)

        Notes
        -----
        ignores config.matchDistanceSigma
        """
        if self.refObjLoader is None:
            raise RuntimeError("Running matcher task with no refObjLoader set in __ini__ or setRefObjLoader")
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        epoch = exposure.visitInfo.date.toAstropy()

        sourceSelection = self.sourceSelector.run(sourceCat)

        sourceFluxField = "slot_%sFlux_instFlux" % (self.config.sourceFluxType)

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=exposure.getBBox(),
            wcs=exposure.wcs,
            filterName=exposure.filter.bandLabel,
            epoch=epoch,
        )

        refSelection = self.referenceSelector.run(loadRes.refCat, exposure=exposure)

        matchMeta = self.refObjLoader.getMetadataBox(
            bbox=exposure.getBBox(),
            wcs=exposure.wcs,
            filterName=exposure.filter.bandLabel,
            epoch=epoch,
        )

        matchRes = self.matcher.matchObjectsToSources(
            refCat=refSelection.sourceCat,
            sourceCat=sourceSelection.sourceCat,
            wcs=exposure.wcs,
            sourceFluxField=sourceFluxField,
            refFluxField=loadRes.fluxField,
            matchTolerance=None,
            bbox=exposure.getBBox(),
        )

        distStats = self._computeMatchStatsOnSky(matchRes.matches)
        self.log.info(
            "Found %d matches with scatter = %0.3f +- %0.3f arcsec; ",
            len(matchRes.matches), distStats.distMean.asArcseconds(), distStats.distStdDev.asArcseconds()
        )

        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=refSelection.sourceCat,
                sourceCat=sourceSelection.sourceCat,
                matches=matchRes.matches,
                exposure=exposure,
                bbox=exposure.getBBox(),
                frame=frame,
                title="Matches",
            )

        return pipeBase.Struct(
            refCat=loadRes.refCat,
            refSelection=refSelection,
            sourceSelection=sourceSelection,
            matches=matchRes.matches,
            matchMeta=matchMeta,
        )

    def _computeMatchStatsOnSky(self, matchList):
        """Compute on-sky radial distance statistics for a match list

        Parameters
        ----------
        matchList : `list` of `lsst.afw.table.ReferenceMatch`
            list of matches between reference object and sources;
            the distance field is the only field read and it must be set to distance in radians

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``distMean`` : clipped mean of on-sky radial separation (`float`)
            - ``distStdDev`` : clipped standard deviation of on-sky radial
              separation (`float`)
            - ``maxMatchDist`` : distMean + self.config.matchDistanceSigma *
              distStdDev (`float`)
        """
        distStatsInRadians = makeMatchStatistics(matchList, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        distMean = distStatsInRadians.getValue(afwMath.MEANCLIP)*lsst.geom.radians
        distStdDev = distStatsInRadians.getValue(afwMath.STDEVCLIP)*lsst.geom.radians
        return pipeBase.Struct(
            distMean=distMean,
            distStdDev=distStdDev,
            maxMatchDist=distMean + self.config.matchDistanceSigma * distStdDev,
        )
