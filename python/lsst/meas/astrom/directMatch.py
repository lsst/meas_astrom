
__all__ = ["DirectMatchConfig", "DirectMatchTask", "DirectMatchConfigWithoutLoader"]

from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.base import Task, Struct
from lsst.meas.algorithms import (LoadIndexedReferenceObjectsTask, ScienceSourceSelectorTask,
                                  ReferenceSourceSelectorTask)
import lsst.afw.table as afwTable
from lsst.afw.geom import arcseconds, averageSpherePoint


class DirectMatchConfigWithoutLoader(Config):
    """Configuration for DirectMatchTask when an already-initialized
    refObjLoader will be passed to this task."""
    matchRadius = Field(dtype=float, default=0.25, doc="Matching radius, arcsec")
    sourceSelection = ConfigurableField(target=ScienceSourceSelectorTask,
                                        doc="Selection of science sources")
    referenceSelection = ConfigurableField(target=ReferenceSourceSelectorTask,
                                           doc="Selection of reference sources")


class DirectMatchConfig(DirectMatchConfigWithoutLoader):
    """Configuration for DirectMatchTask"""
    refObjLoader = ConfigurableField(target=LoadIndexedReferenceObjectsTask, doc="Load reference objects")


class DirectMatchTask(Task):
    r"""!Simple matching of a source catalog to a reference catalog

    @anchor DirectMatchTask_

    @section meas_astrom_match_Contents Contents

     - @ref meas_astrom_match_Purpose
     - @ref meas_astrom_match_Initialize
     - @ref meas_astrom_match_IO
     - @ref meas_astrom_match_Config
     - @ref meas_astrom_match_Example

    @section meas_astrom_match_Purpose  Description

    Match sources to reference objects. The matching permits no rotation or scaling,
    but uses the existing sky positions in the source catalog. This is often useful
    for QA, as it allows validating the pipeline astrometry and photometry against
    the reference catalog.

    Note that this DirectMatchTask is not currently suitable for use within the
    AstrometryTask, as it has a different interface and serves a different purpose.

    @section meas_astrom_match_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_match_IO       Invoking the Task

    @copydoc run

    @section meas_astrom_match_Config       Configuration parameters

    See @ref DirectMatchConfig

    @section meas_astrom_match_Example  A complete example of using DirectMatchTask

    config = DirectMatchConfig()
    task = DirectMatchTask(butler=butler, config=config)
    matchResults = task.run(catalog)

    """

    ConfigClass = DirectMatchConfig
    _DefaultName = "directMatch"

    def __init__(self, butler=None, refObjLoader=None, **kwargs):
        """!Ctor

        Either a 'butler' or 'refObjLoader' is required.

        @param butler  Data butler, or None
        @param refObjLoader  For loading reference objects (lsst.meas.algorithms.LoadReferenceObjectsTask), or
            None
        @param kwargs  Other keyword arguments required for instantiating a Task (e.g., 'config')
        """
        Task.__init__(self, **kwargs)
        if not refObjLoader:
            if not isinstance(self.config, DirectMatchConfig):
                raise RuntimeError("DirectMatchTask must be initialized with DirectMatchConfig "
                                   "if a refObjLoader is not supplied at initialization")
            self.makeSubtask("refObjLoader", butler=butler)
        else:
            self.refObjLoader = refObjLoader
        self.makeSubtask("sourceSelection")
        self.makeSubtask("referenceSelection")

    def run(self, catalog, filterName=None, epoch=None):
        """!Load reference objects and match to them

        @param[in] catalog  Catalog to match to (lsst.afw.table.SourceCatalog)
        @param[in] filterName  Name of filter, for loading fluxes (str)
        @param[in] epoch  Epoch for proper motion and parallax correction
                          (an astropy.time.Time), or None
        @return Struct with matches (lsst.afw.table.SourceMatchVector) and
            matchMeta (lsst.meas.astrom.MatchMetadata)
        """
        circle = self.calculateCircle(catalog)
        matchMeta = self.refObjLoader.getMetadataCircle(circle.center, circle.radius, filterName, epoch=epoch)
        emptyResult = Struct(matches=[], matchMeta=matchMeta)
        sourceSelection = self.sourceSelection.run(catalog)
        if len(sourceSelection.sourceCat) == 0:
            self.log.warn("No objects selected from %d objects in source catalog", len(catalog))
            return emptyResult
        refData = self.refObjLoader.loadSkyCircle(circle.center, circle.radius, filterName, epoch=epoch)
        refCat = refData.refCat
        refSelection = self.referenceSelection.run(refCat)
        if len(refSelection.sourceCat) == 0:
            self.log.warn("No objects selected from %d objects in reference catalog", len(refCat))
            return emptyResult
        matches = afwTable.matchRaDec(refSelection.sourceCat, sourceSelection.sourceCat,
                                      self.config.matchRadius*arcseconds)
        self.log.info("Matched %d from %d/%d input and %d/%d reference sources" %
                      (len(matches), len(sourceSelection.sourceCat), len(catalog),
                       len(refSelection.sourceCat), len(refCat)))
        return Struct(matches=matches, matchMeta=matchMeta, refCat=refCat, sourceSelection=sourceSelection,
                      refSelection=refSelection)

    def calculateCircle(self, catalog):
        """!Calculate a circle enclosing the catalog

        @param[in] catalog  Catalog we will encircle (lsst.afw.table.SourceCatalog)
        @return Struct with ICRS center (lsst.afw.geom.SpherePoint) and radius (lsst.afw.geom.Angle)
        """
        coordList = [src.getCoord() for src in catalog]
        center = averageSpherePoint(coordList)
        radius = max(center.separation(coord) for coord in coordList)
        return Struct(center=center, radius=radius + self.config.matchRadius*arcseconds)
