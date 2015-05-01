from __future__ import absolute_import, division, print_function

import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from lsst.meas.algorithms import setMatchDistance
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .sip import makeCreateWcsWithSip

__all__ = ["FitTanSipWcsTask", "FitTanSipWcsConfig"]

class FitTanSipWcsConfig(pexConfig.Config):
    order = pexConfig.RangeField(
        doc = "order of SIP polynomial (0 for pure TAN WCS)",
        dtype = int,
        default = 3,
        min = 0,
    )

# The following block adds links to this task from the Task Documentation page.
## \addtogroup LSST_task_documentation
## \{
## \page measAstrom_fitTanSipWcsTask
## \ref FitTanSipWcsTask "FitTanSipWcsTask"
##      Fit a TAN-SIP WCS given a list of reference object/source matches
## \}

class FitTanSipWcsTask(pipeBase.Task):
    """!Fit a TAN-SIP WCS given a list of reference object/source matches

    @anchor FitTanSipWcsTask_

    @section meas_astrom_fitTanSipWcs_Contents Contents

     - @ref meas_astrom_fitTanSipWcs_Purpose
     - @ref meas_astrom_fitTanSipWcs_Initialize
     - @ref meas_astrom_fitTanSipWcs_IO
     - @ref meas_astrom_fitTanSipWcs_Schema
     - @ref meas_astrom_fitTanSipWcs_Config
     - @ref meas_astrom_fitTanSipWcs_Example
     - @ref meas_astrom_fitTanSipWcs_Debug

    @section meas_astrom_fitTanSipWcs_Purpose  Description

    Fit a TAN-SIP WCS given a list of reference object/source matches.
    See CreateWithSip.h for information about the fitting algorithm.

    @section meas_astrom_fitTanSipWcs_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_fitTanSipWcs_IO       Invoking the Task

    @copydoc fitWcs

    @section meas_astrom_fitTanSipWcs_Config       Configuration parameters

    See @ref FitTanSipWcsConfig

    @section meas_astrom_fitTanSipWcs_Example  A complete example of using FitTanSipWcsTask

    FitTanSipWcsTask is a subtask of AstrometryTask, which is called by PhotoCalTask.
    See \ref meas_photocal_photocal_Example.

    @section meas_astrom_fitTanSipWcs_Debug        Debug variables

    FitTanSipWcsTask does not support any debug variables.
    """
    ConfigClass = FitTanSipWcsConfig
    _DefaultName = "fitWcs"

    @pipeBase.timeMethod
    def fitWcs(self, matches, initWcs, bbox=None, refCat=None, sourceCat=None):
        """!Fit a TAN-SIP WCS from a list of reference object/source matches

        @param[in,out] matches  a list of reference object/source matches
            (an lsst::afw::table::ReferenceMatchVector)
            The following fields are read:
            - match.first (reference object) coord
            - match.second (source) centroid
            The following fields are written:
            - match.first (reference object) centroid,
            - match.second (source) centroid
            - match.distance (on sky separation, in radians)
        @param[in] initWcs  initial WCS
        @param[in] bbox  the region over which the WCS will be valid (an lsst:afw::geom::Box2I);
            if None or an empty box then computed from matches
        @param[in,out] refCat  reference object catalog, or None.
            If provided then all centroids are updated with the new WCS,
            otherwise only the centroids for ref objects in matches are updated.
            Required fields are "centroid_x", "centroid_y" and "coord".
        @param[in,out] sourceCat  source catalog, or None.
            If provided then coords are updated with the new WCS;
            otherwise only the coords for sources in matches are updated.
            Required fields are "slot_Centroid_x", "slot_Centroid_y" and "coord".

        @return an lsst.pipe.base.Struct with the following fields:
        - wcs  the fit WCS as an lsst.afw.image.Wcs
        - scatterOnSky  median on-sky separation between reference objects and sources in "matches",
            as an lsst.afw.geom.Angle
        """
        if bbox is None:
            bbox = afwGeom.Box2I()
        sipObject = makeCreateWcsWithSip(matches, initWcs, self.config.order, bbox)
        wcs = sipObject.getNewWcs()

        if refCat is not None:
            self.log.info("Updating centroids in refCat")
            self.updateRefCentroids(wcs, refList=refCat)
        else:
            self.log.warning("Updating reference object centroids in match list; refCat is None")
            self.updateRefCentroids(wcs, refList=[match.first for match in matches])

        if sourceCat is not None:
            self.log.info("Updating coords in sourceCat")
            self.updateSourceCoords(wcs, sourceList=sourceCat)
        else:
            self.log.warning("Updating source coords in match list; sourceCat is None")
            self.updateSourceCoords(wcs, sourceList=[match.second for match in matches])

        self.log.info("Updating distance in match list")
        setMatchDistance(matches)

        return pipeBase.Struct(
            wcs = wcs,
            scatterOnSky = sipObject.getScatterOnSky(),
        )

    @staticmethod
    def updateRefCentroids(wcs, refList):
        """Update centroids in a collection of reference objects, given a WCS
        """
        if len(refList) < 1:
            return
        schema = refList[0].schema
        coordKey = schema["coord"].asKey()
        centroidKey = afwTable.Point2DKey(schema["centroid"])
        for refObj in refList:
            refObj.set(centroidKey, wcs.skyToPixel(refObj.get(coordKey)))

    @staticmethod
    def updateSourceCoords(wcs, sourceList):
        """Update coords in a collection of sources, given a WCS
        """
        if len(sourceList) < 1:
            return
        schema = sourceList[1].schema
        srcCoordKey = schema["coord"].asKey()
        for src in sourceList:
            src.set(srcCoordKey, wcs.pixelToSky(src.getCentroid()))

