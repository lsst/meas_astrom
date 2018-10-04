#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ["approximateWcs"]

import numpy as np

import lsst.geom
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.astrom.sip import makeCreateWcsWithSip
from lsst.afw.geom.utils import assertWcsAlmostEqualOverBBox


class _MockTestCase:
    """A fake unit test case class that will enable us to call
    assertWcsAlmostEqualOverBBox from the method approximateWcs"""

    def fail(self, msgStr):
        raise UserWarning("WCS fitting failed " + msgStr)


def approximateWcs(wcs, bbox, order=3, nx=20, ny=20, iterations=3,
                   skyTolerance=0.001*lsst.geom.arcseconds, pixelTolerance=0.02, useTanWcs=False):
    """Approximate an existing WCS as a TAN-SIP WCS

    The fit is performed by evaluating the WCS at a uniform grid of points within a bounding box.

    Parameters
    ----------
    wcs :
        wcs to approximate
    bbox :
        the region over which the WCS will be fit
    order :
        order of SIP fit
    nx :
        number of grid points along x
    ny :
            number of grid points along y
    iterations :
        number of times to iterate over fitting
    skyTolerance :
        maximum allowed difference in world coordinates between
        input wcs and approximate wcs (default is 0.001 arcsec)
    pixelTolerance :
        maximum allowed difference in pixel coordinates between
        input wcs and approximate wcs (default is 0.02 pixels)
    useTanWcs :
        send a TAN version of wcs to the fitter? It is documented to require that,
        but I don't think the fitter actually cares

    Returns
    -------
    fitWcs :
        the fit TAN-SIP WCS
    """
    if useTanWcs:
        crpix = wcs.getPixelOrigin()
        crval = wcs.getSkyOrigin()
        cdMatrix = wcs.getCdMatrix(crpix)
        tanWcs = afwGeom.makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cdMatrix)
    else:
        tanWcs = wcs

    # create a matchList consisting of a grid of points covering the bbox
    refSchema = afwTable.SimpleTable.makeMinimalSchema()
    refCoordKey = afwTable.CoordKey(refSchema["coord"])
    refCat = afwTable.SimpleCatalog(refSchema)

    sourceSchema = afwTable.SourceTable.makeMinimalSchema()
    SingleFrameMeasurementTask(schema=sourceSchema)  # expand the schema
    sourceCentroidKey = afwTable.Point2DKey(sourceSchema["slot_Centroid"])

    sourceCat = afwTable.SourceCatalog(sourceSchema)

    matchList = []

    bboxd = lsst.geom.Box2D(bbox)
    for x in np.linspace(bboxd.getMinX(), bboxd.getMaxX(), nx):
        for y in np.linspace(bboxd.getMinY(), bboxd.getMaxY(), ny):
            pixelPos = lsst.geom.Point2D(x, y)
            skyCoord = wcs.pixelToSky(pixelPos)

            refObj = refCat.addNew()
            refObj.set(refCoordKey, skyCoord)

            source = sourceCat.addNew()
            source.set(sourceCentroidKey, pixelPos)

            matchList.append(afwTable.ReferenceMatch(refObj, source, 0.0))

    # The TAN-SIP fitter is fitting x and y separately, so we have to iterate to make it converge
    for indx in range(iterations):
        sipObject = makeCreateWcsWithSip(matchList, tanWcs, order, bbox)
        tanWcs = sipObject.getNewWcs()
    fitWcs = sipObject.getNewWcs()

    mockTest = _MockTestCase()
    assertWcsAlmostEqualOverBBox(mockTest, wcs, fitWcs, bbox, maxDiffSky=skyTolerance,
                                 maxDiffPix=pixelTolerance)

    return fitWcs
