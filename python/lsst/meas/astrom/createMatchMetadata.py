from __future__ import absolute_import, division, print_function

from lsst.daf.base import PropertyList
from lsst.afw.geom import Box2D
from lsst.afw.image.utils import getDistortedWcs

__all__ = ["createMatchMetadata"]

def createMatchMetadata(exposure, border=0):
    """Create metadata required for unpersisting a match list

    @param[in] exposure  exposure for which to create metadata
    @param[in] border    number of pixels by which to grow the bbox in all directions

    @return metadata about the field (a daf_base PropertyList)
    """
    matchMeta = PropertyList()
    bboxd = Box2D(exposure.getBBox())
    bboxd.grow(border)
    wcs = getDistortedWcs(exposure.getInfo())
    ctrCoord = wcs.pixelToSky(bboxd.getCenter()).toIcrs()
    approxRadius = max(ctrCoord.angularSeparation(wcs.pixelToSky(pp).toIcrs()) for pp in bboxd.getCorners())

    matchMeta.add('RA', ctrCoord.getRa().asDegrees(), 'field center in degrees')
    matchMeta.add('DEC', ctrCoord.getDec().asDegrees(), 'field center in degrees')
    matchMeta.add('RADIUS', approxRadius.asDegrees(), 'field radius in degrees, approximate')
    matchMeta.add('SMATCHV', 1, 'SourceMatchVector version number')
    filterName = exposure.getFilter().getName() or None
    if filterName is not None and filterName not in ("_unknmown_", ""):
        matchMeta.add('FILTER', filterName, 'filter name for tagalong data')
    return matchMeta
