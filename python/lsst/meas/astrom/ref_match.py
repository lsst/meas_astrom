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
        default='Calib'
    )

    def setDefaults(self):
        self.sourceSelector.name = "science"
        self.sourceSelector['science'].fluxLimit.fluxField = \
            'slot_%sFlux_instFlux' % (self.sourceFluxType)
        self.sourceSelector['science'].SignalToNoiseLimit.fluxField = \
            'slot_%sFlux_instFlux' % (self.sourceFluxType)
        self.sourceSelector['science'].SignalToNoiseLimit.errField = \
            'slot_%sFlux_instFluxErr' % (self.sourceFluxType)


class RefMatchTask(pipeBase.Task):
    """Match an input source catalog with objects from a reference catalog.

    Parameters
    ----------
    refObjLoader : `lsst.meas.algorithms.ReferenceLoader`
        A reference object loader object
    **kwargs
        additional keyword arguments for pipe_base `lsst.pipe.base.Task`
    """
    ConfigClass = RefMatchConfig
    _DefaultName = "calibrationBaseClass"

    def __init__(self, refObjLoader, **kwargs):
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
        """Sets the reference object loader for the task

        Parameters
        ----------
        refObjLoader
            An instance of a reference object loader task or class
        """
        self.refObjLoader = refObjLoader

    @pipeBase.timeMethod
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

        sourceSelection = self.sourceSelector.run(sourceCat)

        sourceFluxField = "slot_%sFlux_instFlux" % (self.config.sourceFluxType)

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            photoCalib=expMd.photoCalib,
        )

        refSelection = self.referenceSelector.run(loadRes.refCat)

        matchMeta = self.refObjLoader.getMetadataBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            photoCalib=expMd.photoCalib,
        )

        matchRes = self.matcher.matchObjectsToSources(
            refCat=refSelection.sourceCat,
            sourceCat=sourceSelection.sourceCat,
            wcs=expMd.wcs,
            sourceFluxField=sourceFluxField,
            refFluxField=loadRes.fluxField,
            match_tolerance=None,
        )

        distStats = self._computeMatchStatsOnSky(matchRes.matches)
        self.log.info(
            "Found %d matches with scatter = %0.3f +- %0.3f arcsec; " %
            (len(matchRes.matches), distStats.distMean.asArcseconds(), distStats.distStdDev.asArcseconds())
        )

        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=refSelection.sourceCat,
                sourceCat=sourceSelection.sourceCat,
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
            - ``filterName`` : name of filter (`str`)
            - ``epoch`` : date of exposure (`astropy.time.Time`)

        """
        exposureInfo = exposure.getInfo()
        filterName = exposureInfo.getFilter().getName() or None
        if filterName == "_unknown_":
            filterName = None
        epoch = None
        if exposure.getInfo().hasVisitInfo():
            epochTaiMjd = exposure.getInfo().getVisitInfo().getDate().get(system=DateTime.MJD,
                                                                          scale=DateTime.TAI)
            epoch = astropy.time.Time(epochTaiMjd, scale="tai", format="mjd")

        return pipeBase.Struct(
            bbox=exposure.getBBox(),
            wcs=exposureInfo.getWcs(),
            photoCalib=exposureInfo.getPhotoCalib(),
            filterName=filterName,
            epoch=epoch,
        )
