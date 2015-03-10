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
        default = 5,
        min = 0,
    )

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
    def fitWcs(self, matches, initWcs, bbox=None, refCat=None):
        """!Fit a TAN-SIP WCS from a list of reference object/source matches

        @param[in] matches  a list of reference object/source matches
            (an lsst::afw::table::ReferenceMatchVector)
        @param[in] initWcs  initial WCS
        @param[in] bbox  the region over which the WCS will be valid (an lsst:afw::geom::Box2I);
            if None or an empty box then computed from matches
        @param[in,out] refCat  reference object catalog, or None.
            If provided then all centroids are updated with the new WCS,
            otherwise only the centroids in matches are updated.
            Required fields are "coord", "centroid_x" and "centroid_y".

        @return an lsst.pipe.base.Struct with the following fields:
        - wcs  the fit WCS as an lsst.afw.image.Wcs
        """
        if bbox is None:
            bbox = afwGeom.Box2I()
        sipObject = makeCreateWcsWithSip(matches, initWcs, self.config.order, bbox)
        wcs = sipObject.getNewWcs()

        if refCat is not None:
            coordKey = refCat.schema["coord"].asKey()
            centroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
            for refObj in refCat:
                refObj.set(centroidKey, wcs.skyToPixel(refObj.get(coordKey)))

        return pipeBase.Struct(
            wcs = wcs,
        )
