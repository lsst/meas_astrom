from __future__ import absolute_import, division, print_function

import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
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
            The centroids of the reference objects in matches may be updated for the new WCS,
            but it is strongly recommended that you provide the refCat and sourceCat arguments
            so that all matches are fully updated, as well as any sources or reference objects
            not in matches.
        @param[in] initWcs  initial WCS
        @param[in] bbox  the region over which the WCS will be valid (an lsst:afw::geom::Box2I);
            if None or an empty box then computed from matches
        @param[in,out] refCat  reference object catalog, or None.
            If provided then all centroids are updated with the new WCS,
            otherwise the centroids in matches might be updated but the others will not be touched.
            Required fields are "centroid_x", "centroid_y" and "coord".
        @param[in,out] sourceCat  source catalog, or None.
            If provided then coords are updated with the new WCS.
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
            self.updateRefCat(wcs, refCat)

        if sourceCat is not None:
            self.updateSrcCat(wcs, sourceCat)

        return pipeBase.Struct(
            wcs = wcs,
            scatterOnSky = sipObject.getScatterOnSky(),
        )

    @staticmethod
    def updateRefCat(wcs, refCat):
        """Update centroids in a reference catalog, given a WCS
        """
        coordKey = refCat.schema["coord"].asKey()
        centroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
        for refObj in refCat:
            refObj.set(centroidKey, wcs.skyToPixel(refObj.get(coordKey)))

    @staticmethod
    def updateSrcCat(wcs, sourceCat):
        """Update coords in a source catalog, given a WCS
        """
        srcCoordKey = sourceCat.schema["coord"].asKey()
        for src in sourceCat:
            src.set(srcCoordKey, wcs.pixelToSky(src.getCentroid()))

