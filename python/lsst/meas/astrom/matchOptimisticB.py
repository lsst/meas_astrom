from __future__ import absolute_import, division, print_function
import math

import numpy as np

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .setMatchDistance import setMatchDistance
from .astromLib import matchOptimisticB, MatchOptimisticBControl

__all__ = ["MatchOptimisticBTask", "MatchOptimisticBConfig", "SourceInfo"]

class MatchOptimisticBConfig(pexConfig.Config):
    """Configuration for MatchOptimisticBTask
    """
    sourceFluxType = pexConfig.Field(
        doc = "Type of source flux; typically one of Ap or Psf",
        dtype = str,
        default = "Ap",
    )
    maxMatchDistArcSec = pexConfig.RangeField(
        doc = "Maximum distance between reference objects and sources (arcsec)",
        dtype = float,
        default = 3,
        min = 0,
    )
    numBrightStars = pexConfig.RangeField(
        doc = "Number of bright stars to use",
        dtype = int,
        default = 50,
        min = 2,
    )
    minMatchedPairs = pexConfig.RangeField(
        doc = "Minimum number of matched pairs; see also minFracMatchedPairs",
        dtype = int,
        default = 30,
        min = 2,
    )
    minFracMatchedPairs = pexConfig.RangeField(
        doc = "Minimum number of matched pairs as a fraction of the smaller of "
            "the number of reference stars or the number of good sources; "
            "the actual minimum is the smaller of this value or minMatchedPairs",
        dtype = float,
        default = 0.3,
        min = 0,
        max = 1,
    )
    maxOffsetPix = pexConfig.RangeField(
        doc = "Maximum allowed shift of WCS, due to matching (pixel)",
        dtype = int,
        default = 300,
        max = 4000,
    )
    rotationAllowedInRad = pexConfig.RangeField(
        doc = "Rotation angle allowed between sources and position reference objects (radian)",
        dtype = float,
        default = 0.02,
        max = 0.1,
    )
    angleDiffFrom90 = pexConfig.RangeField(
        doc = "Difference of angle between x and y from 90 degree allowed (degree)",
        dtype = float,
        default = 0.2,
        max = 45.0,
    )
    numPointsForShape = pexConfig.Field(
        doc = "number of points to define a shape for matching",
        dtype = int,
        default = 6,
    )
    maxDeterminant = pexConfig.Field(
        doc = "maximum determinant of linear transformation matrix for a usable solution",
        dtype = float,
        default = 0.02,
    )

class SourceInfo(object):
    """Provide usability tests and catalog keys for sources in a source catalog

    Fields set include:
    - centroidKey  key for centroid
    - centroidFlagKey  key for flag that is True if centroid is valid
    - edgeKey  key for field that is True if source is near an edge
    - saturatedKey  key for field that is True if source has any saturated pixels
    - fluxField  name of flux field

    @throw RuntimeError if schema version unsupported or a needed field is not found
    """
    def __init__(self, schema, fluxType="Ap"):
        """Construct a SourceInfo

        @param[in] schema  source catalog schema
        @param[in] fluxType  flux type: typically one of "Ap" or "Psf"

        @throw RuntimeError if the flux field is not found
        """
        version = schema.getVersion()
        if version == 1:
            self.centroidKey = afwTable.Point2DKey(schema["slot_Centroid"])
            self.centroidFlagKey = schema["slot_Centroid_flag"].asKey()
            self.edgeKey = schema["base_PixelFlags_flag_edge"].asKey()
            self.saturatedKey = schema["base_PixelFlags_flag_saturated"].asKey()
            self.fluxField = "slot_%sFlux_flux" % (fluxType,)
            # extra keys that might be useful
            self.parentKey = schema["parent"].asKey()
            self.interpolatedKey = schema["base_PixelFlags_flag_interpolated"].asKey()
        else:
            raise RuntimeError("Version %r of sourceCat schema not supported" % (version,))

        if self.fluxField not in schema:
            raise RuntimeError("Could not find flux field %s in source schema" % (self.fluxField,))

    def _isResolved(self, source):
        """Return True if source is resolved
        """
        if source.get(self.parentKey) != 0:
            return True
        footprint = source.getFootprint()
        if footprint is not None and len(footprint.getPeaks()) > 1:
            return True
        return False

    def hasCentroid(self, source):
        """Return True if the source has a valid centroid
        """
        centroid = source.get(self.centroidKey)
        return np.all(np.isfinite(centroid)) and not source.getCentroidFlag()

    def isUsable(self, source):
        """Return True if the source is usable for matching, even if possibly saturated

        For a source to be usable it must:
        - have a valid centroid 
        - not be interpolated
        - be not too near the edge
        - not be extended
        - not include cosmic rays
        """
        return (
            self.hasCentroid(source)
            and not source.get(self.edgeKey)
            # and not source.get(self.interpolatedKey)
            # and not self._isResolved(source)
        )

    def isGood(self, source):
        """Return True if source is usable for matching (as per isUsable) and is not saturated
        """
        return self.isUsable(source) and not source.get(self.saturatedKey)


# The following block adds links to this task from the Task Documentation page.
## \addtogroup LSST_task_documentation
## \{
## \page measAstrom_matchOptimisticBTask
## \ref MatchOptimisticBTask "MatchOptimisticBTask"
##      Match sources to reference objects
## \}

class MatchOptimisticBTask(pipeBase.Task):
    """!Match sources to reference objects

    @anchor MatchOptimisticBTask_

    @section meas_astrom_matchOptimisticB_Contents Contents

     - @ref meas_astrom_matchOptimisticB_Purpose
     - @ref meas_astrom_matchOptimisticB_Initialize
     - @ref meas_astrom_matchOptimisticB_IO
     - @ref meas_astrom_matchOptimisticB_Config
     - @ref meas_astrom_matchOptimisticB_Example
     - @ref meas_astrom_matchOptimisticB_Debug

    @section meas_astrom_matchOptimisticB_Purpose  Description

    Match sources to reference objects. This is often done as a preliminary step to fitting an astrometric
    or photometric solution. For details about the matching algorithm see matchOptimisticB.h

    @section meas_astrom_matchOptimisticB_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_matchOptimisticB_IO       Invoking the Task

    @copydoc matchObjectsToSources

    @section meas_astrom_matchOptimisticB_Config       Configuration parameters

    See @ref MatchOptimisticBConfig

    To modify the tests for usable sources and good sources, subclass SourceInfo and
    set MatchOptimisticBTask.SourceInfoClass to your subclass.

    @section meas_astrom_matchOptimisticB_Example  A complete example of using MatchOptimisticBTask

    MatchOptimisticBTask is a subtask of AstrometryTask, which is called by PhotoCalTask.
    See \ref meas_photocal_photocal_Example.

    @section meas_astrom_matchOptimisticB_Debug        Debug variables

    The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink interface supports a
    flag @c -d to import @b debug.py from your @c PYTHONPATH; see @ref baseDebug for more about
    @b debug.py files.

    The available variables in MatchOptimisticBTask are:
    <DL>
      <DT> @c verbose (bool)
      <DD> If True then the matcher prints debug messages to stdout
    </DL>

    To investigate the @ref meas_astrom_matchOptimisticB_Debug, put something like
    @code{.py}
        import lsstDebug
        def DebugInfo(name):
            debug = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
            if name == "lsst.pipe.tasks.astrometry":
                debug.verbose = True

            return debug

        lsstDebug.Info = DebugInfo
    @endcode
    into your debug.py file and run this task with the @c --debug flag.
    """
    ConfigClass = MatchOptimisticBConfig
    _DefaultName = "matchObjectsToSources"
    SourceInfoClass = SourceInfo

    def filterStars(self, refCat):
        """Extra filtering pass; subclass if desired
        """
        return refCat

    @pipeBase.timeMethod
    def matchObjectsToSources(self, refCat, sourceCat, wcs, refFluxField, maxMatchDistArcSec=None):
        """!Match sources to position reference stars

        @param[in] refCat  catalog of reference objects that overlap the exposure; reads fields for:
            - coord
            - the specified flux field
        @param[in] sourceCat  catalog of sources found on an exposure; reads fields for:
            - centroid
            - centroid flag
            - edge flag
            - saturated flag
            - aperture flux, if found, else PSF flux
        @param[in] wcs  estimated WCS
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] maxMatchDistArcSec  maximum distance between reference objects and sources (arcsec);
            if specified then min(config.maxMatchDistArcSec, maxMatchDistArcSec) is used
            if None then config.maxMatchDistArcSec is used
        @return an lsst.pipe.base.Struct with fields:
        - matches  a list of matches, an instance of lsst.afw.table.ReferenceMatch
        - usableSourcCat  a catalog of sources potentially usable for matching.
            For this fitter usable sources include unresolved sources not too near the edge.
            It includes saturated sources, even those these are removed from the final match list,
            because saturated sources may be used to determine the match list.
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        preNumObj = len(refCat)
        refCat = self.filterStars(refCat)
        numRefObj = len(refCat)

        if self.log: self.log.info("filterStars purged %d reference stars, leaving %d stars" % \
            (preNumObj - numRefObj, numRefObj))

        sourceInfo = self.SourceInfoClass(
            schema = sourceCat.schema,
            fluxType = self.config.sourceFluxType,
        )

        # usableSourceCat: sources that are good but may be saturated
        numSources = len(sourceCat)
        usableSourceCat = afwTable.SourceCatalog(sourceCat.table)
        usableSourceCat.extend(s for s in sourceCat if sourceInfo.isUsable(s))
        numUsableSources = len(usableSourceCat)
        self.log.info("Purged %d unsuable sources, leaving %d usable sources" % \
            (numSources - numUsableSources, numUsableSources))

        if len(usableSourceCat) == 0:
            raise pipeBase.TaskError("No sources are usable")

        del sourceCat # avoid accidentally using sourceCat; use usableSourceCat or goodSourceCat from now on

        minMatchedPairs = min(self.config.minMatchedPairs,
                            int(self.config.minFracMatchedPairs * min([len(refCat), len(usableSourceCat)])))

        # match usable (possibly saturated) sources and then purge saturated sources from the match list
        matches = self._doMatch(
            refCat = refCat,
            sourceCat = usableSourceCat,
            wcs = wcs,
            refFluxField = refFluxField,
            numUsableSources = numUsableSources,
            minMatchedPairs = minMatchedPairs,
            maxMatchDistArcSec = maxMatchDistArcSec,
            sourceInfo = sourceInfo,
            verbose = debug.verbose,
        )

        if len(matches) == 0:
            raise RuntimeError("Unable to match sources")

        if self.log: self.log.info("Matched %d sources" % len(matches))
        if len(matches) < minMatchedPairs:
            self.log.warn("Number of matches is smaller than request")

        return pipeBase.Struct(
            matches = matches,
            usableSourceCat = usableSourceCat,
        )

    @pipeBase.timeMethod
    def _doMatch(self, refCat, sourceCat, wcs, refFluxField, numUsableSources, minMatchedPairs,
        maxMatchDistArcSec, sourceInfo, verbose):
        """!Implementation of matching sources to position reference stars

        Unlike matchObjectsToSources, this method does not check if the sources are suitable.

        @param[in] refCat  catalog of position reference stars that overlap an exposure
        @param[in] sourceCat  catalog of sources found on the exposure
        @param[in] wcs  estimated WCS of exposure
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] numUsableSources  number of usable sources (sources with known centroid
            that are not near the edge, but may be saturated)
        @param[in] minMatchedPairs  minimum number of matches
        @param[in] maxMatchDistArcSec  maximum distance between reference objects and sources (arcsec);
            if specified then min(config.maxMatchDistArcSec, maxMatchDistArcSec) is used
            if None then config.maxMatchDistArcSec is used
        @param[in] sourceInfo  SourceInfo for the sourceCat
        @param[in] verbose  true to print diagnostic information to std::cout

        @return a list of matches, an instance of lsst.afw.table.ReferenceMatch
        """
        numSources = len(sourceCat)
        posRefBegInd = numUsableSources - numSources
        if maxMatchDistArcSec is None:
            maxMatchDistArcSec = self.config.maxMatchDistArcSec
        else:
            maxMatchDistArcSec = min(maxMatchDistArcSec, self.config.maxMatchDistArcSec)
        configMatchDistPix = maxMatchDistArcSec/wcs.pixelScale().asArcseconds()

        matchControl = MatchOptimisticBControl()
        matchControl.refFluxField = refFluxField
        matchControl.sourceFluxField = sourceInfo.fluxField
        matchControl.numBrightStars = self.config.numBrightStars
        matchControl.minMatchedPairs = self.config.minMatchedPairs
        matchControl.maxOffsetPix = self.config.maxOffsetPix
        matchControl.numPointsForShape = self.config.numPointsForShape
        matchControl.maxDeterminant = self.config.maxDeterminant

        for maxRotInd in range(4):
            matchControl.maxRotationRad = self.config.rotationAllowedInRad * math.pow(2.0, 0.5*maxRotInd)
            for matchRadInd in range(3):
                matchControl.matchingAllowancePix = configMatchDistPix * math.pow(1.25, matchRadInd)
                                  
                for angleDiffInd in range(3):
                    matchControl.angleDiffFrom90 = self.config.angleDiffFrom90*(angleDiffInd+1)
                    matches = matchOptimisticB(
                        refCat,
                        sourceCat,
                        matchControl,
                        posRefBegInd,
                        verbose,
                    )
                    if matches is not None and len(matches) > 0:
                        setMatchDistance(matches)
                        return matches
        return matches
