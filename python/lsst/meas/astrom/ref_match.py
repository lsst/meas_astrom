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
from lsst.meas.algorithms import ScienceSourceSelectorTask, ReferenceSourceSelectorTask
from .matchOptimisticB import MatchOptimisticBTask
from .display import displayAstrometry
from . import makeMatchStatistics


class RefMatchConfig(pexConfig.Config):
    matcher = pexConfig.ConfigurableField(
        target=MatchOptimisticBTask,
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
    sourceSelection = pexConfig.ConfigurableField(target=ScienceSourceSelectorTask,
                                                  doc="Selection of science sources")
    referenceSelection = pexConfig.ConfigurableField(target=ReferenceSourceSelectorTask,
                                                     doc="Selection of reference sources")

# The following block adds links to this task from the Task Documentation page.
## \addtogroup LSST_task_documentation
## \{
## \page measAlgorithms_RefMatchTask
## \ref RefMatchTask_ "RefMatchTask"
##         Basic functionality for all calibration tasks: i.e. a matcher
## \}


class RefMatchTask(pipeBase.Task):
    """Match an input source catalog with objects from a reference catalog

    Parameters
    ----------
    refObjLoader :
        A reference object loader object
    schema :
        ignored; available for compatibility with an older astrometry task
    kwargs :
        additional keyword arguments for pipe_base Task.\_\_init\_\_
    """
    ConfigClass = RefMatchConfig
    _DefaultName = "calibrationBaseClass"

    def __init__(self, refObjLoader, schema=None, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.refObjLoader = refObjLoader
        self.makeSubtask("matcher")
        self.makeSubtask("sourceSelection")
        self.makeSubtask("referenceSelection")

    @pipeBase.timeMethod
    def loadAndMatch(self, exposure, sourceCat):
        """Load reference objects overlapping an exposure and match to sources detected on that exposure

        Parameters
        ----------
        exposure :
            exposure that the sources overlap
        sourceCat :
            catalog of sources detected on the exposure (an lsst.afw.table.SourceCatalog)

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            with these fields:

            - ``refCat`` :  reference object catalog of objects that overlap the exposure (with some margin)
                (an lsst::afw::table::SimpleCatalog)
            - ``matches`` :  a list of lsst.afw.table.ReferenceMatch
            - ``matchMeta`` :  metadata needed to unpersist matches (an lsst.daf.base.PropertyList)

        Notes
        -----
        ignores config.matchDistanceSigma
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        expMd = self._getExposureMetadata(exposure)

        sourceSelection = self.sourceSelection.run(sourceCat)

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            calib=expMd.calib,
        )

        refSelection = self.referenceSelection.run(loadRes.refCat)

        matchMeta = self.refObjLoader.getMetadataBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            calib=expMd.calib,
        )

        matchRes = self.matcher.matchObjectsToSources(
            refCat=refSelection.sourceCat,
            sourceCat=sourceSelection.sourceCat,
            wcs=expMd.wcs,
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
        matchList :
            list of matches between reference object and sources;
            the distance field is the only field read and it must be set to distance in radians

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            a pipe_base Struct containing these fields:

            - ``distMean`` :  clipped mean of on-sky radial separation
            - ``distStdDev`` :  clipped standard deviation of on-sky radial separation
            - ``maxMatchDist`` :  distMean + self.config.matchDistanceSigma*distStdDev
        """
        distStatsInRadians = makeMatchStatistics(matchList, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        distMean = distStatsInRadians.getValue(afwMath.MEANCLIP)*lsst.geom.radians
        distStdDev = distStatsInRadians.getValue(afwMath.STDEVCLIP)*lsst.geom.radians
        return pipeBase.Struct(
            distMean=distMean,
            distStdDev=distStdDev,
            maxMatchDist=distMean + self.config.matchDistanceSigma*distStdDev,
        )

    def _getExposureMetadata(self, exposure):
        """Extract metadata from an exposure

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            containing the following exposure metadata:

            - ``bbox`` : parent bounding box
            - ``wcs`` : WCS (an lsst.afw.geom.Wcs)
            - ``calib`` : calibration (an lsst.afw.image.Calib), or None if unknown
            - ``filterName`` : name of filter, or None if unknown
            - ``epoch`` : date of exposure (an astropy.time.Time), or None

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
            calib=exposureInfo.getCalib() if exposureInfo.hasCalib() else None,
            filterName=filterName,
            epoch=epoch,
        )
