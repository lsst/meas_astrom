from __future__ import absolute_import, division, print_function

__all__ = ["denormalizeMatches"]

import lsst.afw.table


def denormalizeMatches(matches, matchMeta=None):
    """Generate a denormalized Catalog of matches

    This is intended for writing matches in a convenient way.
    Normally we write matches in a 'normalized' form: recording only the join
    table (reference ID, source ID) to minimise space (the reference and source
    catalogs should both be available separately, so the only extra information
    we need is how to join them). However, using that can be a pain, since it
    requires reading each catalog and doing the join.

    This function generates a Catalog containing all the information in the
    matches. The reference catalog entries are in columns with "ref_"
    prepended, while the source catalog entries are in columns with "src_"
    prepended. The distance between the matches is in a column named
    "distance".

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

    See also
    --------
    `lsst.afw.table.packMatches`
    """
    if len(matches) == 0:
        raise RuntimeError("No matches provided.")

    refSchema = matches[0].first.getSchema()
    srcSchema = matches[0].second.getSchema()

    refMapper, srcMapper = lsst.afw.table.SchemaMapper.join([refSchema, srcSchema], ["ref_", "src_"])
    schema = refMapper.editOutputSchema()
    distKey = schema.addField("distance", type=float, doc="Distance between ref and src")

    catalog = lsst.afw.table.BaseCatalog(schema)
    catalog.reserve(len(matches))
    for mm in matches:
        row = catalog.addNew()
        row.assign(mm.first, refMapper)
        row.assign(mm.second, srcMapper)
        row.set(distKey, mm.distance)

    if matchMeta is not None:
        catalog.getTable().setMetadata(matchMeta)

    return catalog
