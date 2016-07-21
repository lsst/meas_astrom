import lsst.meas.algorithms.ingestIndexReferenceTask
assert type(config)==lsst.meas.algorithms.ingestIndexReferenceTask.IngestIndexedReferenceConfig, 'config is of type %s.%s instead of lsst.meas.algorithms.ingestIndexReferenceTask.IngestIndexedReferenceConfig' % (type(config).__module__, type(config).__name__)
import astropy.utils.xml.validate
import astropy._erfa.core
import astropy.io.fits.card
import astropy.table
import astropy.io.fits.py3compat
import astropy.io.fits.column
import astropy.io.votable.table
import astropy.table.np_utils
import astropy.io.fits.hdu.compressed
import astropy.io.fits.diff
import astropy.time.formats
import astropy.io.fits
import astropy.io.ascii.misc
import astropy.io.votable.exceptions
import astropy.io.ascii.connect
import astropy.io.fits.hdu.groups
import astropy.table.table
import numpy.lib.recfunctions
import astropy.utils.compat.numpy.lib
import astropy.time.utils
import astropy.io.ascii.core
import astropy.time.core
import astropy.io.fits.verify
import astropy.io.ascii.basic
import astropy.io.fits._numpy_hacks
import astropy.io.ascii.daophot
import astropy.io.ascii.ecsv
import astropy.io.ascii.fixedwidth
import astropy.utils.xml.writer
import astropy.io.ascii.html
import astropy.io.ascii.ui
import astropy.table.operations
import astropy.utils.xml._iterparser
import astropy.io.fits.file
import astropy.utils.compat.gzip
import astropy.io.misc.pickle_helpers
import astropy.io.votable
import astropy.io.fits.header
import astropy.table.info
import astropy.io.fits.hdu.streaming
import astropy.io.misc.hdf5
import astropy.io.fits.hdu.hdulist
import astropy.io.votable.converters
import astropy.utils.collections
import lsst.meas.algorithms.readFitsCatalogTask
import astropy.table.meta
import astropy._erfa
import astropy.time.erfa_time
import astropy.io.fits.compression
import astropy.io.votable.connect
import astropy.io.votable.ucd
import astropy.io.misc
import astropy.io.fits.fitsrec
import astropy.utils.data
import astropy.io.fits.connect
import astropy.io.misc.connect
import astropy.io.fits.hdu.base
import astropy.table.pprint
import astropy.table.row
import astropy.io.fits.hdu.nonstandard
import astropy.utils.compat.numpy.lib.stride_tricks
import astropy.io.votable.tablewriter
import astropy.io.ascii.ipac
import astropy.utils.xml
import astropy.utils.xml.iterparser
import astropy.io.ascii.latex
import astropy.io.registry
import astropy.io.ascii
import astropy.io
import astropy.io.ascii.cparser
import astropy.table._np_utils
import astropy.io.fits.util
import backports_abc
import astropy.table.column
import astropy.table._column_mixins
import astropy.utils.compat.numpy
import astropy.io.fits.hdu.table
import astropy.table.sorted_array
import astropy.utils.xml.check
import astropy.io.fits.convenience
import astropy.io.ascii.cds
import astropy.io.ascii.fastbasic
import numpy.ma.mrecords
import astropy.table.groups
import astropy.time
import astropy.io.fits.hdu
import astropy.io.fits.hdu.image
import astropy.utils.metadata
import astropy._erfa._core
import astropy.io.ascii.sextractor
import astropy.table.index
import astropy.io.votable.tree
import astropy.io.votable.xmlutil
import astropy.io.votable.util
import astropy.table.jsviewer
import astropy.table.bst
# Name of column stating if satisfactory for photometric calibration (optional).
config.is_photometric_name='photometric'

# Default HTM level.  Level 8 gives ~0.08 sq deg per trixel.
config.level=8

# Name of RA column
config.ra_name='ra'

# Name of Dec column
config.dec_name='dec'

import lsst.meas.algorithms.readFitsCatalogTask
config.file_reader.retarget(target=lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogTask, ConfigClass=lsst.meas.algorithms.readFitsCatalogTask.ReadFitsCatalogConfig)
# HDU containing the desired binary table, 0-based but a binary table never occurs in HDU 0
config.file_reader.hdu=1

# Mapping of input column name: output column name; each specified column must exist, but additional columns in the input data are written using their original name. 
config.file_reader.column_map={}

# Name of column stating if the object is resolved (optional).
config.is_resolved_name='resolved'

# Name of column to use as an identifier (optional).
config.id_name='id'

# The values in the reference catalog are assumed to be in AB magnitudes. List of column names to use for photometric information.  At least one entry is required.
config.mag_column_list=['u', 'g', 'r', 'i', 'z']

# A map of magnitude column name (key) to magnitude error column (value).
config.mag_err_column_map={'i': 'i_err', 'r': 'r_err', 'u': 'u_err', 'z': 'z_err', 'g': 'g_err'}

# String to pass to the butler to retrieve persisted files.
config.ref_dataset_name='cal_ref_cat'

# Extra columns to add to the reference catalog.
config.extra_col_names=[]

# Name of column stating if the object is measured to be variable (optional).
config.is_variable_name=None

