
__all__ = ["DirectMatchConfig", "DirectMatchTask", "DirectMatchConfigWithoutLoader"]

from lsst.pex.config import Config, Field, ConfigurableField, ConfigField
from lsst.pipe.base import Task, Struct
from lsst.meas.algorithms import (LoadReferenceObjectsConfig, ScienceSourceSelectorTask,
                                  ReferenceSourceSelectorTask)
import lsst.afw.table as afwTable
from lsst.geom import arcseconds, averageSpherePoint


class DirectMatchConfigWithoutLoader(Config):
    """Configuration for `DirectMatchTask` when an already-initialized
    ``refObjLoader`` will be passed to this task.
    """
    matchRadius = Field(dtype=float, default=0.25, doc="Matching radius, arcsec")
    doSourceSelection = Field(
        dtype=bool,
        doc="Select sources to be matched with `sourceSelector`?"
        " Set to False if you want to use exactly the sources that are passed in.",
        default=True,
    )
    sourceSelection = ConfigurableField(target=ScienceSourceSelectorTask,
                                        doc="Selection of science sources")
    referenceSelection = ConfigurableField(target=ReferenceSourceSelectorTask,
                                           doc="Selection of reference sources")


class DirectMatchConfig(DirectMatchConfigWithoutLoader):
    """Configuration for `DirectMatchTask`.
    """
    refObjLoader = ConfigField(dtype=LoadReferenceObjectsConfig,
                               doc="Configuration of reference object loader")


class DirectMatchTask(Task):
    """Simple, brute force matching of a source catalog to a reference catalog.

    Parameters
    ----------
    butler : `None`
        Compatibility parameter. Should not be used.
    refObjLoader : `lsst.meas.algorithms.ReferenceObjectLoader` or `None`
        A reference object loader object; gen3 pipeline tasks will pass `None`
        and call `setRefObjLoader` in `runQuantum`.
    **kwargs
        Other keyword arguments required for instantiating a Task (such as
        ``config``).
    """
    ConfigClass = DirectMatchConfig
    _DefaultName = "directMatch"

    def __init__(self, refObjLoader=None, **kwargs):
        Task.__init__(self, **kwargs)
        self.refObjLoader = refObjLoader
        if self.config.doSourceSelection:
            self.makeSubtask("sourceSelection")
        self.makeSubtask("referenceSelection")

    def setRefObjLoader(self, refObjLoader):
        """Set the reference object loader for the task.

        Parameters
        ----------
        refObjLoader : `lsst.meas.algorithms.ReferenceObjectLoader`
            An instance of a reference object loader.
        """
        self.refObjLoader = refObjLoader

    def run(self, catalog, filterName=None, epoch=None):
        """Load reference objects and match to them.

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog`
            Catalog to match.
        filterName : `str`
            Name of filter loading fluxes.
        epoch : `astropy.time.Time` or `None`
            Epoch to which to correct proper motion and parallax, or `None` to
            not apply such corrections.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            ``matches``
                Matched sources with associated reference
                (`lsst.afw.table.SourceMatchVector`).
            ``matchMeta``
                Match metadata (`lsst.meas.astrom.MatchMetadata`).
        """
        if self.refObjLoader is None:
            raise RuntimeError("Running matcher task with no refObjLoader set in __ini__ or setRefObjLoader")
        circle = self.calculateCircle(catalog)
        matchMeta = self.refObjLoader.getMetadataCircle(circle.center, circle.radius, filterName, epoch=epoch)

        emptyResult = Struct(matches=[], matchMeta=matchMeta)
        if self.config.doSourceSelection:
            sourceSelection = self.sourceSelection.run(catalog)
            if len(sourceSelection.sourceCat) == 0:
                self.log.warning("No objects selected from %d objects in source catalog", len(catalog))
                return emptyResult
            sourceCat = sourceSelection.sourceCat
        else:
            sourceSelection = None
            sourceCat = catalog

        refData = self.refObjLoader.loadSkyCircle(circle.center, circle.radius, filterName, epoch=epoch)
        refCat = refData.refCat
        refSelection = self.referenceSelection.run(refCat)
        if len(refSelection.sourceCat) == 0:
            self.log.warning("No objects selected from %d objects in reference catalog", len(refCat))
            return emptyResult
        matches = afwTable.matchRaDec(refSelection.sourceCat, sourceCat, self.config.matchRadius*arcseconds)
        self.log.info("Matched %d from %d/%d input and %d/%d reference sources",
                      len(matches), len(sourceCat), len(catalog),
                      len(refSelection.sourceCat), len(refCat))
        return Struct(matches=matches, matchMeta=matchMeta, refCat=refCat, sourceSelection=sourceSelection,
                      refSelection=refSelection)

    def calculateCircle(self, catalog):
        """Calculate a circle enclosing the catalog.

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog`
            Catalog to encircle.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            ``center``
                ICRS center coordinate (`lsst.afw.geom.SpherePoint`).
            ``radius``
                Radius of the circle (`lsst.geom.Angle`).
        """
        coordList = [src.getCoord() for src in catalog]
        center = averageSpherePoint(coordList)
        radius = max(center.separation(coord) for coord in coordList)
        return Struct(center=center, radius=radius + self.config.matchRadius*arcseconds)
