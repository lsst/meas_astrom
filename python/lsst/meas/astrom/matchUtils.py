from __future__ import absolute_import, division, print_function
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable

__all__ = ["setMatchDistance", "joinMatchListWithCatalog"]

def setMatchDistance(matches):
    """Set the distance field of the matches in a match list to the distance in radians on the sky

    @warning the coord field of the source in each match must be correct

    @param[in,out] matches  a list of matches, an instance of lsst.afw.table.ReferenceMatch
        reads the coord field of the source and reference object of each match
        writes the distance field of each match
    """
    if len(matches) < 1:
        return

    sourceCoordKey = afwTable.CoordKey(matches[0].first.schema["coord"])
    refObjCoordKey = afwTable.CoordKey(matches[0].second.schema["coord"])
    for match in matches:
        sourceCoord = match.first.get(sourceCoordKey)
        refObjCoord = match.second.get(refObjCoordKey)
        match.distance = refObjCoord.angularSeparation(sourceCoord).asRadians()


def joinMatchListWithCatalog(packedMatches, sourceCat, astrom=None):
    """
    This function is required to reconstitute a ReferenceMatchVector after being
    unpersisted.  The persisted form of a ReferenceMatchVector is the
    normalized Catalog of IDs produced by afw.table.packMatches(), with the result of
    InitialAstrometry.getMatchMetadata() in the associated tables\' metadata.

    The "live" form of a matchlist has links to
    the real record objects that are matched; it is "denormalized".
    This function takes a normalized match catalog, along with the catalog of
    sources to which the match catalog refers.  It fetches the reference
    sources that are within range, and then denormalizes the matches
    -- sets the "matches[*].first" and "matches[*].second" entries
    to point to the sources in the "sourceCat" argument, and to the
    reference sources fetched from the astrometry_net_data files.

    @param[in] packedMatches  Unpersisted match list (an lsst.afw.table.BaseCatalog).
                              packedMatches.table.getMetadata() must contain the
                              values from InitialAstrometry.getMatchMetadata()
    @param[in,out] sourceCat  Source catalog used for the 'second' side of the matches
                              (an lsst.afw.table.SourceCatalog).  As a side effect,
                              the catalog will be sorted by ID.
    @param[in] astrom         Initialized Astrometry task (must have an refObjLoader.loadSkyCircle()
                              function).  If None, meas.astrom.AstrometryTask will be initialized.

    @return An lsst.afw.table.ReferenceMatchVector of denormalized matches.
    """
    if astrom is None:
        from .astrometry import AstrometryTask, AstrometryConfig
        astrom = AstrometryTask(AstrometryConfig)

    matchmeta = packedMatches.table.getMetadata()
    version = matchmeta.getInt('SMATCHV')
    if version != 1:
        raise ValueError('SourceMatchVector version number is %i, not 1.' % version)
    filterName = matchmeta.getString('FILTER').strip()
    ctrCoord = afwCoord.IcrsCoord(
        matchmeta.getDouble('RA') * afwGeom.degrees,
        matchmeta.getDouble('DEC') * afwGeom.degrees,
        )
    rad = matchmeta.getDouble('RADIUS') * afwGeom.degrees
    refCat = astrom.refObjLoader.loadSkyCircle(ctrCoord, rad, filterName).refCat
    refCat.sort()
    sourceCat.sort()
    return afwTable.unpackMatches(packedMatches, refCat, sourceCat)
