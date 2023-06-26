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

__all__ = ['RefMatchConfig', 'RefMatchTask']

import astropy.time

import lsst.geom
from lsst.daf.base import DateTime
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
    doSourceSelection = pexConfig.Field(
        dtype=bool,
        doc="Select sources to be matched with `sourceSelector`?"
            " Set to False if you want to use exactly the sources that are passed in.",
        default=True
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources for cross-matching.",
        default="science",
        optional=True
    )
    referenceSelector = pexConfig.ConfigurableField(
        target=ReferenceSourceSelectorTask,
        doc="How to select reference objects for cross-matching."
    )
    sourceFluxType = pexConfig.Field(
        dtype=str,
        doc="Source flux type to use in source selection.",
        default='Calib'
    )

    def setDefaults(self):
        self.sourceSelector['science'].fluxLimit.fluxField = \
            'slot_%sFlux_instFlux' % (self.sourceFluxType)
        self.sourceSelector['science'].signalToNoise.fluxField = \
            'slot_%sFlux_instFlux' % (self.sourceFluxType)
        self.sourceSelector['science'].signalToNoise.errField = \
            'slot_%sFlux_instFluxErr' % (self.sourceFluxType)


class RefMatchTask(pipeBase.Task):
    """Match an input source catalog with objects from a reference catalog.

    Parameters
    ----------
    refObjLoader : `lsst.meas.algorithms.ReferenceLoader`
        A reference object loader object; gen3 pipeline tasks will pass `None`
        and call `setRefObjLoader` in `runQuantum`.
    **kwargs
        additional keyword arguments for pipe_base `lsst.pipe.base.Task`
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
        if self.config.doSourceSelection:
            self.makeSubtask("sourceSelector")
        self.makeSubtask("referenceSelector")

    def setRefObjLoader(self, refObjLoader):
        """Sets the reference object loader for the task

        Parameters
        ----------
        refObjLoader
            An instance of a reference object loader task or class
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

        expMd = self._getExposureMetadata(exposure)

        if self.config.doSourceSelection:
            sourceSelection = self.sourceSelector.run(sourceCat)
            catalog = sourceSelection.sourceCat
        else:
            sourceSelection = None
            catalog = sourceCat

        sourceFluxField = "slot_%sFlux_instFlux" % (self.config.sourceFluxType)

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

        matchRes = self.matcher.matchObjectsToSources(
            refCat=refSelection.sourceCat,
            sourceCat=catalog,
            wcs=expMd.wcs,
            sourceFluxField=sourceFluxField,
            refFluxField=loadRes.fluxField,
            match_tolerance=None,
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
                sourceCat=catalog,
                matches=matchRes.matches,
                exposure=exposure,
                bbox=expMd.bbox,
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

    def _getExposureMetadata(self, exposure):
        """Extract metadata from an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``bbox`` : parent bounding box (`lsst.geom.Box2I`)
            - ``wcs`` : exposure WCS (`lsst.afw.geom.SkyWcs`)
            - ``photoCalib`` : photometric calibration (`lsst.afw.image.PhotoCalib`)
            - ``filterName`` : name of filter band (`str`)
            - ``epoch`` : date of exposure (`astropy.time.Time`)
        """
        filterLabel = exposure.info.getFilter()
        filterName = filterLabel.bandLabel if filterLabel is not None else None
        epoch = None
        if exposure.info.hasVisitInfo():
            epochTaiMjd = exposure.visitInfo.date.get(system=DateTime.MJD, scale=DateTime.TAI)
            epoch = astropy.time.Time(epochTaiMjd, scale="tai", format="mjd")

        return pipeBase.Struct(
            bbox=exposure.getBBox(),
            wcs=exposure.info.getWcs(),
            photoCalib=exposure.info.getPhotoCalib(),
            filterName=filterName,
            epoch=epoch,
        )
