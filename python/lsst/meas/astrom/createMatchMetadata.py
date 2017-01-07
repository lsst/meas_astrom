from __future__ import absolute_import, division, print_function

from lsst.daf.base import PropertyList
from lsst.afw.geom import Box2D
from lsst.afw.image.utils import getDistortedWcs

__all__ = ["MatchMetadata", "createMatchMetadata"]


class MatchMetadata(PropertyList):
    """Metadata required for unpersisting a match list"""

    def __init__(self, ctrCoord, radius, filterName):
        """!Ctor

        @param[in] ctrCoord: Coordinates of center (lsst.afw.coord.IcrsCoord)
        @param[in] radius: Minimum radius for selecting sources (lsst.afw.geom.Angle)
        @param[in] filterName: Name of filter (str) or None
        """
        PropertyList.__init__(self)
        ctrCoord = ctrCoord.toIcrs()
        self.add('RA', ctrCoord.getRa().asDegrees(), 'field center in degrees')
        self.add('DEC', ctrCoord.getDec().asDegrees(), 'field center in degrees')
        self.add('RADIUS', radius.asDegrees(), 'field radius in degrees, minimum')
        self.add('SMATCHV', 1, 'SourceMatchVector version number')
        filterName = "UNKNOWN" if filterName is None else str(filterName)
        self.add('FILTER', filterName, 'filter name for photometric data')


def createMatchMetadata(exposure, border=0):
    """Create metadata required for unpersisting a match list

    @param[in] exposure  exposure for which to create metadata
    @param[in] border    number of pixels by which to grow the bbox in all directions

    @return metadata about the field (a daf_base PropertyList)
    """
    bboxd = Box2D(exposure.getBBox())
    bboxd.grow(border)
    wcs = getDistortedWcs(exposure.getInfo())
    ctrCoord = wcs.pixelToSky(bboxd.getCenter()).toIcrs()
    approxRadius = max(ctrCoord.angularSeparation(wcs.pixelToSky(pp).toIcrs()) for pp in bboxd.getCorners())
    return MatchMetadata(ctrCoord, approxRadius, exposure.getFilter().getName())
