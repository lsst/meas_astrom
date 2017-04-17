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
from __future__ import absolute_import, division, print_function

__all__ = ['RefMatchConfig', 'RefMatchTask']

from lsst.afw.image.utils import getDistortedWcs
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .matchOptimisticB import MatchOptimisticBTask
from .display import displayAstrometry
from . import makeMatchStatistics
from .createMatchMetadata import createMatchMetadata


class RefMatchConfig(pexConfig.Config):
    matcher = pexConfig.ConfigurableField(
        target=MatchOptimisticBTask,
        doc="reference object/source matcher",
    )
    matchDistanceSigma = pexConfig.RangeField(
        doc="the maximum match distance is set to "
        " mean_match_distance + matchDistanceSigma*std_dev_match_distance; " +
        "ignored if not fitting a WCS",
        dtype=float,
        default=2,
        min=0,
    )

# The following block adds links to this task from the Task Documentation page.
## \addtogroup LSST_task_documentation
## \{
## \page measAlgorithms_RefMatchTask
## \ref RefMatchTask_ "RefMatchTask"
##         Basic functionality for all calibration tasks: i.e. a matcher
## \}


class RefMatchTask(pipeBase.Task):
    """!Match an input source catalog with objects from a reference catalog

    @anchor RefMatchTask_
    """
    ConfigClass = RefMatchConfig
    _DefaultName = "calibrationBaseClass"

    def __init__(self, refObjLoader, schema=None, **kwargs):
        """!Construct a RefMatchTask

        @param[in] refObjLoader A reference object loader object
        @param[in] schema  ignored; available for compatibility with an older astrometry task
        @param[in] kwargs  additional keyword arguments for pipe_base Task.\_\_init\_\_
        """
        pipeBase.Task.__init__(self, **kwargs)
        self.refObjLoader = refObjLoader
        self.makeSubtask("matcher")

    @pipeBase.timeMethod
    def loadAndMatch(self, exposure, sourceCat):
        """!Load reference objects overlapping an exposure and match to sources detected on that exposure

        @param[in] exposure  exposure that the sources overlap
        @param[in] sourceCat  catalog of sources detected on the exposure (an lsst.afw.table.SourceCatalog)

        @return an lsst.pipe.base.Struct with these fields:
        - refCat  reference object catalog of objects that overlap the exposure (with some margin)
            (an lsst::afw::table::SimpleCatalog)
        - matches  a list of lsst.afw.table.ReferenceMatch
        - matchMeta  metadata needed to unpersist matches (an lsst.daf.base.PropertyList)

        @note ignores config.matchDistanceSigma
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        matchMeta = createMatchMetadata(exposure, border=self.refObjLoader.config.pixelMargin)
        expMd = self._getExposureMetadata(exposure)

        loadRes = self.refObjLoader.loadPixelBox(
            bbox=expMd.bbox,
            wcs=expMd.wcs,
            filterName=expMd.filterName,
            calib=expMd.calib,
        )

        matchRes = self.matcher.matchObjectsToSources(
            refCat=loadRes.refCat,
            sourceCat=sourceCat,
            wcs=expMd.wcs,
            refFluxField=loadRes.fluxField,
            toleranceStruct=None,
        )

        distStats = self._computeMatchStatsOnSky(matchRes.matches)
        self.log.info(
            "Found %d matches with scatter = %0.3f +- %0.3f arcsec; " %
            (len(matchRes.matches), distStats.distMean.asArcseconds(), distStats.distStdDev.asArcseconds())
        )

        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat=loadRes.refCat,
                sourceCat=sourceCat,
                matches=matchRes.matches,
                exposure=exposure,
                bbox=expMd.bbox,
                frame=frame,
                title="Matches",
            )

        return pipeBase.Struct(
            refCat=loadRes.refCat,
            matches=matchRes.matches,
            matchMeta=matchMeta,
        )

    def _computeMatchStatsOnSky(self, matchList):
        """Compute on-sky radial distance statistics for a match list

        @param[in] matchList  list of matches between reference object and sources;
            the distance field is the only field read and it must be set to distance in radians

        @return a pipe_base Struct containing these fields:
        - distMean  clipped mean of on-sky radial separation
        - distStdDev  clipped standard deviation of on-sky radial separation
        - maxMatchDist  distMean + self.config.matchDistanceSigma*distStdDev
        """
        distStatsInRadians = makeMatchStatistics(matchList, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        distMean = distStatsInRadians.getValue(afwMath.MEANCLIP)*afwGeom.radians
        distStdDev = distStatsInRadians.getValue(afwMath.STDEVCLIP)*afwGeom.radians
        return pipeBase.Struct(
            distMean=distMean,
            distStdDev=distStdDev,
            maxMatchDist=distMean + self.config.matchDistanceSigma*distStdDev,
        )

    def _getExposureMetadata(self, exposure):
        """!Extract metadata from an exposure

        @return an lsst.pipe.base.Struct containing the following exposure metadata:
        - bbox: parent bounding box
        - wcs: WCS (an lsst.afw.image.Wcs)
        - calib calibration (an lsst.afw.image.Calib), or None if unknown
        - filterName: name of filter, or None if unknown
        """
        exposureInfo = exposure.getInfo()
        filterName = exposureInfo.getFilter().getName() or None
        if filterName == "_unknown_":
            filterName = None
        return pipeBase.Struct(
            bbox=exposure.getBBox(),
            wcs=getDistortedWcs(exposureInfo, log=self.log),
            calib=exposureInfo.getCalib() if exposureInfo.hasCalib() else None,
            filterName=filterName,
        )
