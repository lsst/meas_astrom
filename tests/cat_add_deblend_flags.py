# Add deblend and primary fields, and centroid flag to the the astrometry test catalog.
import numpy as np
from lsst.afw.table import SourceCatalog, SchemaMapper

catalog = SourceCatalog.readFits("cat.xy.fits")
mapper = SchemaMapper(catalog.schema)
mapper.addMinimalSchema(catalog.schema)
# nChild field can remain 0; we don't want to change existing test behavior.
mapper.addOutputField("deblend_nChild", type=np.int32,
                      doc='Number of children this object has (defaults to 0)')
# We want all sources to be primary sources.
mapper.addOutputField("detect_isPrimary", type=np.int32,
                      doc="true if source has no children and is not a sky source")
# We want all sources to have valid centroids (i.e. flag=0).
mapper.addOutputField("base_SdssCentroid_flag", type="Flag",
                      doc="General centroid failure flag")

schema = mapper.getOutputSchema()
schema.setAliasMap(catalog.schema.getAliasMap())

new_catalog = SourceCatalog(schema)
new_catalog.extend(catalog, mapper=mapper)
new_catalog['detect_isPrimary'] = 1

new_catalog.writeFits("cat.xy.fits")
