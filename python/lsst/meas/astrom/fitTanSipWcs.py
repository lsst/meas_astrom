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
    """!Fit a TAN-SIP WCS

    See CreateWithSip.h for information about the fitting algorithm.
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
