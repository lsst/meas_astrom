from __future__ import absolute_import, division, print_function

import os

import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import LoadReferenceObjectsTask, getRefFluxField
from . import astrometry_net as astromNet
from .astrometryNetDataConfig import AstrometryNetDataConfig

__all__ = ["LoadAstrometryNetObjectsTask", "LoadAstrometryNetObjectsConfig"]

LoadAstrometryNetObjectsConfig = LoadReferenceObjectsTask.ConfigClass

# The following block adds links to this task from the Task Documentation page.
## \addtogroup LSST_task_documentation
## \{
## \page measAstrom_loadAstrometryNetObjectsTask
## \ref LoadAstrometryNetObjectsTask "LoadAstrometryNetObjectsTask"
##      Load reference objects from astrometry.net index files
## \}

class LoadAstrometryNetObjectsTask(LoadReferenceObjectsTask):
    """!Load reference objects from astrometry.net index files

    @anchor LoadAstrometryNetObjectsTask_

    @section meas_astrom_loadAstrometryNetObjects_Contents Contents

     - @ref meas_astrom_loadAstrometryNetObjects_Purpose
     - @ref meas_astrom_loadAstrometryNetObjects_Initialize
     - @ref meas_astrom_loadAstrometryNetObjects_IO
     - @ref meas_algorithms_loadReferenceObjects_Schema
     - @ref meas_astrom_loadAstrometryNetObjects_Config
     - @ref meas_astrom_loadAstrometryNetObjects_Example
     - @ref meas_astrom_loadAstrometryNetObjects_Debug

    @section meas_astrom_loadAstrometryNetObjects_Purpose  Description

    Load reference objects from astrometry.net index files.

    @section meas_astrom_loadAstrometryNetObjects_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_loadAstrometryNetObjects_IO       Invoking the Task

    @copydoc loadObjectsInBBox

    @section meas_astrom_loadAstrometryNetObjects_Config       Configuration parameters

    See @ref LoadAstrometryNetObjectsConfig

    @section meas_astrom_loadAstrometryNetObjects_Example  A complete example of using LoadAstrometryNetObjectsTask

    LoadAstrometryNetObjectsTask is a subtask of AstrometryTask, which is called by PhotoCalTask.
    See \ref meas_photocal_photocal_Example.

    @section meas_astrom_loadAstrometryNetObjects_Debug        Debug variables

    LoadAstrometryNetObjectsTask does not support any debug variables.
    """
    ConfigClass = LoadAstrometryNetObjectsConfig

    def __init__(self, config, andConfig=None, **kwargs):
        """!Create a LoadAstrometryNetObjectsTask

        @param[in] config  configuration (an instance of self.ConfigClass)
        @param[in] andConfig  astrometry.net data config (an instance of AstromNetDataConfig, or None);
            if None then use andConfig.py in the astrometry_net_data product (which must be setup)

        @throw RuntimeError if andConfig is None and the configuration cannot be found,
            either because astrometry_net_data is not setup in eups
            or because the setup version does not include the file "andConfig.py"
        """
        LoadReferenceObjectsTask.__init__(self, config=config, **kwargs)
        self.andConfig = andConfig
        self.haveIndexFiles = False # defer reading index files until we know they are needed
            # because astrometry may not be used, in which case it may not be properly configured

    @pipeBase.timeMethod
    def loadSkyCircle(self, ctrCoord, radius, filterName=None):
        """!Load reference objects that overlap a circular sky region

        @param[in] ctrCoord  center of search region (an afwGeom.Coord)
        @param[in] radius  radius of search region (an afwGeom.Angle)
        @param[in] filterName  name of filter, or None for the default filter;
            used for flux values in case we have flux limits (which are not yet implemented)

        @return an lsst.pipe.base.Struct containing:
        - refCat a catalog of reference objects with the
            \link meas_algorithms_loadReferenceObjects_Schema standard schema \endlink
            as documented in LoadReferenceObjects, including photometric, resolved and variable;
            hasCentroid is False for all objects.
        - fluxField = name of flux field for specified filterName
        """
        self._readIndexFiles()

        names = []
        mcols = []
        ecols = []
        for col, mcol in self.andConfig.magColumnMap.items():
            names.append(col)
            mcols.append(mcol)
            ecols.append(self.andConfig.magErrorColumnMap.get(col, ''))
        margs = (names, mcols, ecols)

        solver = self._getSolver()

        # Find multi-index files within range
        multiInds = self._getMIndexesWithinRange(ctrCoord, radius)

        # compute solver.getCatalog arguments that follow the list of star kd-trees:
        # - center equatorial angle (e.g. RA) in deg
        # - center polar angle (e.g. Dec) in deg
        # - radius, in deg
        # - idColumn
        # - (margs)
        # - star-galaxy column
        # - variability column
        fixedArgTuple = (
            ctrCoord,
            radius,
            self.andConfig.idColumn,
        ) + margs + (
            self.andConfig.starGalaxyColumn,
            self.andConfig.variableColumn,
            True, # eliminate duplicate IDs
        )

        self.log.info("search for objects at %s with radius %s deg" % (ctrCoord, radius.asDegrees()))
        with LoadMultiIndexes(multiInds):
            # We just want to pass the star kd-trees, so just pass the
            # first element of each multi-index.
            inds = tuple(mi[0] for mi in multiInds)
            refCat = solver.getCatalog(inds, *fixedArgTuple)

        self._addFluxAliases(schema=refCat.schema)

        fluxField = getRefFluxField(schema=refCat.schema, filterName=filterName)

        self.log.info("found %d objects" % (len(refCat),))
        return pipeBase.Struct(
            refCat = refCat,
            fluxField = fluxField,
        )

    @pipeBase.timeMethod
    def _readIndexFiles(self):
        """!Read all astrometry.net index files, if not already read
        """
        if self.haveIndexFiles:
            return

        self.log.info("read index files")

        self.multiInds = []
        self.haveIndexFiles = True # just try once

        if self.andConfig is None:
            # use andConfig.py in the astrometry_net product setup in eups
            anDir = lsst.utils.getPackageDir('astrometry_net_data')
            if anDir is None:
                raise RuntimeError("astrometry_net_data is not setup")

            andConfig = AstrometryNetDataConfig()
            andConfigPath = os.path.join(anDir, "andConfig.py")
            if not os.path.exists(andConfigPath):
                raise RuntimeError("astrometry_net_data config file \"%s\" required but not found" %
                    andConfigPath)
            andConfig.load(andConfigPath)
            self.andConfig = andConfig

        # merge indexFiles and multiIndexFiles; we'll treat both as multiindex for simplicity.
        mifiles = [(True, [fn,fn]) for fn  in self.andConfig.indexFiles] + \
            [(False, fns) for fns in self.andConfig.multiIndexFiles]

        nMissing = 0
        for single, fns in mifiles:
            # First filename in "fns" is star kdtree, the rest are index files.
            fn = fns[0]
            if single:
                self.log.log(self.log.DEBUG, 'Adding index file %s' % fns[0])
            else:
                self.log.log(self.log.DEBUG, 'Adding multiindex files %s' % str(fns))
            fn2 = self._getIndexPath(fn)
            if fn2 is None:
                if single:
                    self.log.logdebug('Unable to find index file %s' % fn)
                else:
                    self.log.logdebug('Unable to find star part of multiindex file %s' % fn)
                nMissing += 1
                continue
            fn = fn2
            self.log.log(self.log.DEBUG, 'Path: %s' % fn)

            mi = astromNet.multiindex_new(fn)
            if mi is None:
                raise RuntimeError('Failed to read objects from multiindex filename "%s"' % fn)
            for i,fn in enumerate(fns[1:]):
                self.log.log(self.log.DEBUG, 'Reading index from multiindex file "%s"' % fn)
                fn2 = self._getIndexPath(fn)
                if fn2 is None:
                    self.log.logdebug('Unable to find index part of multiindex file %s' % fn)
                    nMissing += 1
                    continue
                fn = fn2
                self.log.log(self.log.DEBUG, 'Path: %s' % fn)
                if astromNet.multiindex_add_index(mi, fn, astromNet.INDEX_ONLY_LOAD_METADATA):
                    raise RuntimeError('Failed to read index from multiindex filename "%s"' % fn)
                ind = mi[i]
                self.log.log(self.log.DEBUG, '  index %i, hp %i (nside %i), nstars %i, nquads %i' %
                                (ind.indexid, ind.healpix, ind.hpnside, ind.nstars, ind.nquads))
            astromNet.multiindex_unload_starkd(mi)
            self.multiInds.append(mi)

        if len(self.multiInds) == 0:
            self.log.warn('Unable to find any index files')
        elif nMissing > 0:
            self.log.warn('Unable to find %d index files' % (nMissing,))

    def _getIndexPath(self, fn):
        """!Get the path to the specified astrometry.net index file

        @param[in] fn  path to index file; if relative, then relative to astrometry_net_data
            if that product is setup, else relative to the current working directory
        @return the absolute path to the index file, or None if the file was not found
        """
        if os.path.isabs(fn):
            absFn = fn
        else:
            anDir = lsst.utils.getPackageDir('astrometry_net_data')
            if anDir is not None:
                absFn = os.path.join(anDir, fn)

        if os.path.exists(absFn):
            return os.path.abspath(absFn)
        else:
            return None

    def _getMIndexesWithinRange(self, ctrCoord, radius):
        """!Get list of muti-index objects within range

        @param[in] ctrCoord  center of search region (an afwGeom.Coord)
        @param[in] radius  radius of search region (an afwGeom.Angle)

        @return list of multiindex objects
        """
        longDeg  = ctrCoord.getLongitude().asDegrees()
        latDeg = ctrCoord.getLatitude().asDegrees()
        multiIndexList = []
        for mi in self.multiInds:
            if mi.isWithinRange(longDeg, latDeg, radius.asDegrees()):
                multiIndexList.append(mi)
        return multiIndexList

    def _getSolver(self):
        solver = astromNet.solver_new()
        # HACK, set huge default pixel scale range.
        lo,hi = 0.01, 3600.
        solver.setPixelScaleRange(lo, hi)
        return solver


class LoadMultiIndexes(object):
    """Context manager for loading astrometry.net multi-index files, then unloading and finalizing a.net
    """
    def __init__(self, multiInds):
        self.multiInds = multiInds
    def __enter__(self):
        for mi in self.multiInds:
            mi.reload()
        return self.multiInds
    def __exit__(self, typ, val, trace):
        for mi in self.multiInds:
            mi.unload()
        astromNet.finalize()
