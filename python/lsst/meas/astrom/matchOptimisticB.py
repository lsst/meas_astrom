from __future__ import absolute_import, division, print_function
import math
import numpy

import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .astromLib import matchOptimisticB, MatchOptimisticBControl

__all__ = ["MatchOptimisticBTask", "MatchOptimisticBConfig"]

class MatchOptimisticBConfig(pexConfig.Config):
    sourceFluxType = pexConfig.Field(
        doc = "Type of source flux; typically one of Ap or Psf",
        dtype = str,
        default = "Ap",
    )
    maxMatchDistArcSec = pexConfig.RangeField(
        doc = "Maximum radius for matching sources to reference objects (arcsec)",
        dtype = float,
        default = 1,
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
        doc = "Maximum offset between sources and position reference objects (pixel)",
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
    """Provide information about sources in a source catalog for various versions of the schema

    Fields set include:
    - edgeKey  key for field that is True if source is near an edge
    - saturatedKey  key for field that is True if source has any saturated pixels
    - centroidKey  key for centroid
    - centroidFlagKey  key for flag that is True if centroid is valid
    - fluxField  name of flux field

    @throw RuntimeError if schema version unsupported or a needed field not found
    """
    def __init__(self, schema, fluxType="Ap"):
        """Construct a SourceInfo

        @param[in] schema  source catalog schema
        @param[in] fluxType  flux type: typically one of "Ap" or "Psf"

        @throw RuntimeError if the flux field is not found
        """
        version = schema.getVersion()
        if version == 0:
            self.edgeKey = schema["flags.pixel.edge"].asKey()
            self.saturatedKey = schema["flags.pixel.saturated.any"].asKey()
            self.centroidKey = afwTable.Point2DKey(schema["slot.Centroid"])
            self.centroidFlagKey = schema["slot.Centroid.flags"].asKey()
            self.fluxField = "slot.%sFlux" % (fluxType,)
        elif version == 1:
            self.edgeKey = schema["base_PixelFlags_flag_edge"].asKey()
            self.saturatedKey = schema["base_PixelFlags_flag_saturated"].asKey()
            self.centroidKey = afwTable.Point2DKey(schema["slot_Centroid"])
            self.centroidFlagKey = schema["slot_Centroid_flag"].asKey()
            self.fluxField = "slot_%sFlux_flux" % (fluxType,)
        else:
            raise RuntimeError("Version %r of sourceCat schema not supported" % (version,))

        if self.fluxField not in schema:
            raise RuntimeError("Could not find flux field %s in source schema" % (self.fluxField,))

    def hasCentroid(self, source):
        """Return True if the source has a valid centroid
        """
        centroid = source.get(self.centroidKey)
        return numpy.all(numpy.isfinite(centroid)) and not source.getCentroidFlag()

    def isCleanSource(self, source):
        """Return True if the source has a valid centroid and is not near the edge
        """
        return self.hasCentroid(source) and not source.get(self.edgeKey)

    def isGoodSource(self, source):
        """Return True if source is clean (as per isCleanSource) and is not saturated
        """
        return self.isCleanSource(source) and not source.get(self.saturatedKey)


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

    def filterStars(self, refCat):
        """Extra filtering pass; subclass if desired
        """
        return refCat

    @pipeBase.timeMethod
    def matchObjectsToSources(self, refCat, sourceCat, wcs, refFluxField):
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
        @return an lsst.pipe.base.Struct with fields:
        - matches  a list of matches, an instance of lsst.afw.table.ReferenceMatch
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        preNumObj = len(refCat)
        refCat = self.filterStars(refCat)
        numRefObj = len(refCat)

        if self.log: self.log.info("filterStars purged %d reference stars, leaving %d stars" % \
            (preNumObj - numRefObj, numRefObj))

        sourceInfo = SourceInfo(schema=sourceCat.schema, fluxType=self.config.sourceFluxType)

        # cleanSourceCat: sources that are good but may be saturated
        numSources = len(sourceCat)
        cleanSourceCat = afwTable.SourceCatalog(sourceCat.table)
        cleanSourceCat.extend(s for s in sourceCat if sourceInfo.isCleanSource(s))
        numCleanSources = len(cleanSourceCat)
        self.log.info("Purged %d unsuable sources, leaving %d clean sources" % \
            (numSources - numCleanSources, numCleanSources))

        # goodSourceCat: sources that are good and not saturated
        goodSources = afwTable.SourceCatalog(sourceCat.table)
        goodSources.extend(s for s in cleanSourceCat if sourceInfo.isGoodSource(s))
        numGoodSources = len(goodSources)
        self.log.info("Purged %d saturated sources, leaving %d good sources" % \
            (numCleanSources - numGoodSources, numGoodSources))

        del sourceCat # avoid accidentally using sourceCat; use cleanSourceCat or goodSources from now on
        
        self.log.info("Matching to %d/%d good input sources" % (len(goodSources), len(cleanSourceCat)))

        minMatchedPairs = min(self.config.minMatchedPairs,
                            int(self.config.minFracMatchedPairs * min([len(refCat), len(goodSources)])))

        # match clean (possibly saturated) sources and then purge saturated sources from the match list
        matches0 = self._doMatch(
            refCat = refCat,
            sourceCat = cleanSourceCat,
            wcs = wcs,
            refFluxField = refFluxField,
            numCleanSources = numCleanSources,
            minMatchedPairs = minMatchedPairs,
            sourceInfo = sourceInfo,
            verbose = debug.verbose,
        )
        if matches0 is not None:
            matches0 = [m for m in matches0 if sourceInfo.isGoodSource(m.second)]
        else:
            matches0 = []

        # match good (unsaturated) sources (i.e. prefilter saturated sources)
        if len(refCat) > len(cleanSourceCat) - len(goodSources):
            matches1 = self._doMatch(
                refCat = refCat,
                sourceCat = goodSources,
                wcs = wcs,
                refFluxField = refFluxField,
                numCleanSources = numCleanSources,
                minMatchedPairs = minMatchedPairs,
                sourceInfo = sourceInfo,
                verbose = debug.verbose,
            )
        else:
            matches1 = []
        if matches1 == None:
            matches1 = []

        if len(matches0) == 0 and len(matches1) == 0:
            raise RuntimeError("Unable to match sources")

        # Adopt matches with more matches
        if len(matches0) > len(matches1):
            matches = matches0
        else:
            matches = matches1

        if self.log: self.log.info("Matched %d sources" % len(matches))
        if len(matches) < minMatchedPairs:
            self.log.warn("Number of matches is smaller than request")

        return pipeBase.Struct(
            matches = matches,
        )

    @pipeBase.timeMethod
    def _doMatch(self, refCat, sourceCat, wcs, refFluxField, numCleanSources, minMatchedPairs,
        sourceInfo, verbose):
        """!Implementation of matching sources to position reference stars

        Unlike matchObjectsToSources, this method does not check if the sources are suitable.

        @param[in] refCat  catalog of position reference stars that overlap an exposure
        @param[in] sourceCat  catalog of sources found on the exposure
        @param[in] wcs  estimated WCS of exposure
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] numCleanSources  number of clean sources (sources with known centroid
            that are not near the edge, but may be saturated)
        @param[in] minMatchedPairs  minimum number of matches
        @param[in] sourceInfo  SourceInfo for the sourceCat
        @param[in] verbose  true to print diagnostic information to std::cout
        """
        numSources = len(sourceCat)
        posRefBegInd = numCleanSources - numSources
        configMatchDistPix = self.config.maxMatchDistArcSec/wcs.pixelScale().asArcseconds()

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
                    if matches is not None and len(matches) != 0:
                        return matches
