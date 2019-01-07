# This file is part of meas_astrom.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["denormalizeMatches"]

import lsst.afw.table as afwTable


def denormalizeMatches(matches, matchMeta=None):
    """Generate a denormalized Catalog of matches

    Parameters
    ----------
    matches : `list` of `lsst.afw.table.ReferenceMatch`
        List of matches between reference catalog and source catalog.
    matchMeta : `lsst.daf.base.PropertyList`
        Matching metadata to write in catalog.

    Returns
    -------
    catalog : `lsst.afw.table.BaseCatalog`
        Catalog containing matchlist entries.

    Notes
    -----
    This is intended for writing matches in a convenient way.
    Normally we write matches in a 'normalized' form: recording only the join
    table (reference ID, source ID) to minimise space (the reference and source
    catalogs should both be available separately, so the only extra information
    we need is how to join them). However, using that can be a pain, since it
    requires reading each catalog and doing the join.

    This function generates a Catalog containing all the information in the
    matches. The reference catalog entries are in columns with 'ref'
    prepended, while the source catalog entries are in columns with 'src'
    prepended (including any alias mappings). The distance between the
    matches is in a column named "distance".

    See Also
    --------
    lsst.afw.table.packMatches
    """
    # TODO: DM-16863 Current this link is removed due to the conversion of
    # afw.table not yet being complete and causing an error on build.
    # """
    if len(matches) == 0:
        raise RuntimeError("No matches provided.")

    refSchema = matches[0].first.getSchema()
    srcSchema = matches[0].second.getSchema()

    refMapper, srcMapper = afwTable.SchemaMapper.join([refSchema, srcSchema], ["ref_", "src_"])
    schema = refMapper.editOutputSchema()

    schema = afwTable.catalogMatches.copyAliasMapWithPrefix(srcSchema, schema, prefix="src_")
    schema = afwTable.catalogMatches.copyAliasMapWithPrefix(refSchema, schema, prefix="ref_")

    distKey = schema.addField("distance", type=float, doc="Distance between ref and src")

    catalog = afwTable.BaseCatalog(schema)
    catalog.reserve(len(matches))
    for mm in matches:
        row = catalog.addNew()
        row.assign(mm.first, refMapper)
        row.assign(mm.second, srcMapper)
        row.set(distKey, mm.distance)

    if matchMeta is not None:
        catalog.getTable().setMetadata(matchMeta)

    return catalog
