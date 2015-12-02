# 
# LSST Data Management System
#
# Copyright 2008-2015 AURA/LSST.
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
import lsst.afw.table as afwTable
from lsst.meas.astrom.anetBasicAstrometry import ANetBasicAstrometryTask, ANetBasicAstrometryConfig

__all__ = ["readMatches"]

def readMatches(butler, dataId, sourcesName='icSrc', matchesName='icMatch',
                config=ANetBasicAstrometryConfig(), sourcesFlags=afwTable.SOURCE_IO_NO_FOOTPRINTS):
    """Read matches, sources and catalogue; combine.
    \param[in] butler Data butler
    \param[in] dataId Data identifier for butler
    \param[in] sourcesName Name for sources from butler
    \param[in] matchesName Name for matches from butler
    \param[in] sourcesFlags Flags to pass for source retrieval
    \returns Matches
    """
    sources = butler.get(sourcesName, dataId, flags=sourcesFlags)
    packedMatches = butler.get(matchesName, dataId)
    astrom = ANetBasicAstrometryTask(config)
    return astrom.joinMatchListWithCatalog(packedMatches, sources)
