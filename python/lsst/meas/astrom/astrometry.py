from __future__ import absolute_import, division, print_function

from lsst.daf.base import PropertyList
from lsst.afw.image.utils import getDistortedWcs
from lsst.afw.geom import Box2D, radToArcsec
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .loadAstrometryNetObjects import LoadAstrometryNetObjectsTask
from .matchOptimisticB import MatchOptimisticBTask
from .fitTanSipWcs import FitTanSipWcsTask
from .display import displayAstrometry

class AstrometryConfig(pexConfig.Config):
    refObjLoader = pexConfig.ConfigurableField(
        target = LoadAstrometryNetObjectsTask,
        doc = "reference object loader",
    )
    matcher = pexConfig.ConfigurableField(
        target = MatchOptimisticBTask,
        doc = "reference object/source matcher",
    )
    wcsFitter = pexConfig.ConfigurableField(
        target = FitTanSipWcsTask,
        doc = "WCS fitter",
    )
    forceKnownWcs = pexConfig.Field(
        dtype = bool,
        doc= "Assume that the input image's WCS is correct, without comparing it to any external reality " +
            "  (but still match reference objects to sources)",
        default = False,
    )
    maxIter = pexConfig.RangeField(
        doc = "maximum number of iterations of match sources and fit WCS; " +
            "ignored if forceKnownWcs True",
        dtype = int,
        default = 6,
        min = 1,
    )

# The following block adds links to this task from the Task Documentation page.
## \addtogroup LSST_task_documentation
## \{
## \page measAstrom_astrometryTask
## \ref AstrometryTask_ "AstrometryTask"
##      Match an input source catalog with objects from a reference catalog and solve for the WCS
## \}

class AstrometryTask(pipeBase.Task):
    """!Match an input source catalog with objects from a reference catalog and solve for the WCS

    @anchor AstrometryTask_

    @section meas_astrom_astrometry_Contents Contents

     - @ref meas_astrom_astrometry_Purpose
     - @ref meas_astrom_astrometry_Initialize
     - @ref meas_astrom_astrometry_IO
     - @ref meas_astrom_astrometry_Config
     - @ref meas_astrom_astrometry_Example
     - @ref meas_astrom_astrometry_Debug

    @section meas_astrom_astrometry_Purpose  Description

    Match input sourceCat with a reference catalog and solve for the Wcs

    There are three steps, each performed by different subtasks:
    - Find position reference stars that overlap the exposure
    - Match sourceCat to position reference stars
    - Fit a WCS based on the matches

    @section meas_astrom_astrometry_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_astrometry_IO       Invoking the Task

    @copydoc run

    @section meas_astrom_astrometry_Config       Configuration parameters

    See @ref AstrometryConfig

    @section meas_astrom_astrometry_Example  A complete example of using AstrometryTask

    See \ref meas_photocal_photocal_Example.

    @section meas_astrom_astrometry_Debug        Debug variables

    The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink interface supports a
    flag @c -d to import @b debug.py from your @c PYTHONPATH; see @ref baseDebug for more about
    @b debug.py files.

    The available variables in AstrometryTask are:
    <DL>
      <DT> @c display (bool)
      <DD> If True display information at three stages: after finding reference objects,
        after matching sources to reference objects, and after fitting the WCS; defaults to False
      <DT> @c frame (int)
      <DD> ds9 frame to use to display the reference objects; the next two frames are used
            to display the match list and the results of the final WCS; defaults to 0
    </DL>

    To investigate the @ref meas_astrom_astrometry_Debug, put something like
    @code{.py}
    import lsstDebug
    def DebugInfo(name):
        debug = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
        if name == "lsst.meas.astrom.astrometry":
            debug.display = True

        return debug

    lsstDebug.Info = DebugInfo
    @endcode
    into your debug.py file and run this task with the @c --debug flag.
    """
    ConfigClass = AstrometryConfig
    _DefaultName = "astrometricSolver"

    def __init__(self, schema=None, **kwargs):
        """!Construct an AstrometryTask

        @param[in] schema  ignored; available for compatibility with an older astrometry task
        """
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("refObjLoader")
        self.makeSubtask("matcher")
        self.makeSubtask("wcsFitter")

    @pipeBase.timeMethod
    def run(self, exposure, sourceCat):
        """!Fit a WCS given a source catalog and exposure

        @param[in,out] exposure  exposure whose WCS is to be fit
            The following are read only:
            - bbox
            - calib (may be absent)
            - filter (may be unset)
            - detector (if wcs is pure tangent; may be absent)
            The following are updated:
            - wcs (the initial value is used as an initial guess, and is required)
        @param[in] sourceCat  catalog of sourceCat detected on the exposure (an lsst.afw.table.SourceCatalog)
        @return the same struct as the "solve" method
        """
        bbox = exposure.getBBox()
        exposureInfo = exposure.getInfo()
        initWcs = getDistortedWcs(exposureInfo, log=self.log)
        calib = exposureInfo.getCalib() if exposureInfo.hasCalib() else None
        filterName = exposureInfo.getFilter().getName() or None

        retVal = self.solve(sourceCat=sourceCat, bbox=bbox, initWcs=initWcs, filterName=filterName,
            calib=calib, exposure=exposure)
        exposure.setWcs(retVal.wcs)
        return retVal

    @pipeBase.timeMethod
    def solve(self, sourceCat, bbox, initWcs, filterName=None, calib=None, exposure=None):
        """!Fit a WCS given a source catalog and exposure matchMetadata

        @param[in] sourceCat  catalog of sourceCat detected on the exposure (an lsst.afw.table.SourceCatalog)
        @param[in] bbox  bounding box of exposure (an lsst.afw.geom.Box2I)
        @param[in] wcs  initial guess for WCS of exposure (an lsst.afw.image.Wcs)
        @param[in] filterName  filter name, or None or "" if unknown (a string)
        @param[in] calib  calibration for exposure, or None if unknown (an lsst.afw.image.Calib)
        @param[in] exposure  exposure whose WCS is to be fit, or None; used only for the debug display

        @return an lsst.pipe.base.Struct with these fields:
        - refCat  reference object catalog of objects that overlap the exposure (with some margin)
            (an lsst::afw::table::SimpleCatalog)
        - matches  list of reference object/source matches (an lsst.afw.table.ReferenceMatchVector)
        - initWcs  initial WCS from exposure (possibly tweaked) (an lsst.afw.image.Wcs)
        - wcs  fit WCS (an lsst.afw.image.Wcs)
        - matchMeta  metadata about the field
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        loadRes = self.refObjLoader.loadPixelBox(
            bbox = bbox,
            wcs = initWcs,
            filterName = filterName,
            calib = calib,
        )
        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat = loadRes.refCat,
                sourceCat = sourceCat,
                exposure = exposure,
                bbox = bbox,
                frame = frame,
                title="Reference catalog",
            )

        res = None
        wcs = initWcs
        maxMatchDistArcSec = None
        for i in range(self.config.maxIter):
            try:
                tryRes = self._matchAndFitWcs( # refCat, sourceCat, refFluxField, bbox, wcs, exposure=None
                    refCat = loadRes.refCat,
                    sourceCat = sourceCat,
                    refFluxField = loadRes.fluxField,
                    bbox = bbox,
                    wcs = wcs,
                    exposure = exposure,
                    maxMatchDistArcSec = maxMatchDistArcSec,
                )
            except Exception as e:
                # if we have had a succeessful iteration then use that; otherwise fail
                if i > 0:
                    self.log.info("Fit WCS iter %d failed; using previous iteration: %s" % (i, e))
                    break
                else:
                    raise

            if self.config.forceKnownWcs:
                # run just once; note that the number of matches has already logged
                res = tryRes
                break

            distRadList = [match.distance for match in tryRes.matches]
            distRadStats = afwMath.makeStatistics(distRadList, afwMath.MEANCLIP | afwMath.STDEVCLIP)
            distArcsecMean = radToArcsec(distRadStats.getValue(afwMath.MEANCLIP))
            distArcsecStdDev = radToArcsec(distRadStats.getValue(afwMath.STDEVCLIP))

            self.log.info("Fit WCS iter %d: %d matches; scatter = %0.3f +- %0.3f arcsec" %
                    (i, len(tryRes.matches), distArcsecMean, distArcsecStdDev))

            newMaxMatchDistArcSec = distArcsecMean + 2.0*distArcsecStdDev
            if maxMatchDistArcSec is not None:
                if newMaxMatchDistArcSec >= maxMatchDistArcSec:
                    self.log.warn(
                        "Iteration %d had no better mean + 2sigma scatter; using previous iteration" % (i,))
                    break

            maxMatchDistArcSec = newMaxMatchDistArcSec
            res = tryRes
            wcs = res.wcs
            if newMaxMatchDistArcSec < 0.001:
                self.log.info(
                    "Iteration %d had mean + 2sigma scatter < 0.001 arcsec; that's good enough" % (i,))
                break

        return pipeBase.Struct(
            refCat = loadRes.refCat,
            matches = res.matches,
            initWcs = initWcs,
            wcs = res.wcs,
            scatterOnSky = res.scatterOnSky,
            matchMeta = self._createMatchMetadata(bbox=bbox, wcs=res.wcs, filterName=filterName)
        )

    @pipeBase.timeMethod
    def _matchAndFitWcs(self, refCat, sourceCat, refFluxField, bbox, wcs, maxMatchDistArcSec=None,
        exposure=None):
        """!Match sources to reference objects and fit a WCS

        @param[in] refCat  catalog of reference objects
        @param[in] sourceCat  catalog of sourceCat detected on the exposure (an lsst.afw.table.SourceCatalog)
        @param[in] bbox  bounding box of exposure (an lsst.afw.geom.Box2I)
        @param[in] wcs  initial guess for WCS of exposure (an lsst.afw.image.Wcs)
        @param[in] maxMatchDistArcSec  maximum distance between reference objects and sources (arcsec);
            if None then use the matcher's default
        @param[in] exposure  exposure whose WCS is to be fit, or None; used only for the debug display

        @return an lsst.pipe.base.Struct with these fields:
        - matches  list of reference object/source matches (an lsst.afw.table.ReferenceMatchVector)
        - wcs  the fit WCS as an lsst.afw.image.Wcs
        - scatterOnSky  median on-sky separation between reference objects and sources in "matches",
            as an lsst.afw.geom.Angle, or None if config.forceKnownWcs is True
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)
        matchRes = self.matcher.matchObjectsToSources(
            refCat = refCat,
            sourceCat = sourceCat,
            wcs = wcs,
            refFluxField = refFluxField,
            maxMatchDistArcSec = maxMatchDistArcSec,
        )
        self.log.info("Found %s matches" % (len(matchRes.matches),))
        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat = refCat,
                sourceCat = matchRes.usableSourceCat,
                matches = matchRes.matches,
                exposure = exposure,
                bbox = bbox,
                frame = frame + 1,
                title="Initial WCS",
            )

        if not self.config.forceKnownWcs:
            self.log.info("Fitting WCS")
            fitRes = self.wcsFitter.fitWcs(
                matches = matchRes.matches,
                initWcs = wcs,
                bbox = bbox,
                refCat = refCat,
                sourceCat = sourceCat,
            )
            fitWcs = fitRes.wcs
            scatterOnSky = fitRes.scatterOnSky
        else:
            self.log.info("Not fitting WCS (forceKnownWcs true)")
            fitWcs = wcs
            scatterOnSky = None
        if debug.display:
            frame = int(debug.frame)
            displayAstrometry(
                refCat = refCat,
                sourceCat = matchRes.usableSourceCat,
                matches = matchRes.matches,
                exposure = exposure,
                bbox = bbox,
                frame = frame + 2,
                title="Fit TAN-SIP WCS",
            )

        return pipeBase.Struct(
            matches = matchRes.matches,
            wcs = fitWcs,
            scatterOnSky = scatterOnSky,
        )

    @staticmethod
    def _createMatchMetadata(bbox, wcs, filterName):
        """Create matchMeta metadata required for regenerating the catalog

        This is copied from Astrom and I'm not sure why it is needed.
        I did not put this in any subtask because none have all the necessary info:
        - The matcher does not have the fit wcs
        - The fitter does not have the filter name

        @param bbox  bounding box of exposure (an lsst.afw.geom.Box2I or Box2D)
        @param wcs  WCS of exposure
        @param filterName Name of filter, used for magnitudes
        @return Metadata
        """
        matchMeta = PropertyList()
        bboxd = Box2D(bbox)
        ctrPos = bboxd.getCenter()
        ctrCoord = wcs.pixelToSky(ctrPos).toIcrs()
        llCoord = wcs.pixelToSky(bboxd.getMin())
        approxRadius = ctrCoord.angularSeparation(llCoord)
        matchMeta.add('RA', ctrCoord.getRa().asDegrees(), 'field center in degrees')
        matchMeta.add('DEC', ctrCoord.getDec().asDegrees(), 'field center in degrees')
        matchMeta.add('RADIUS', approxRadius.asDegrees(), 'field radius in degrees, approximate')
        matchMeta.add('SMATCHV', 1, 'SourceMatchVector version number')
        if filterName is not None:
            matchMeta.add('FILTER', filterName, 'filter name for tagalong data')
        return matchMeta
