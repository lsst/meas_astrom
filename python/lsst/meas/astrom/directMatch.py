
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
    """Simple, brute force matching of a source catalog to a reference catalog.

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler`
        Data butler containing the relevant reference catalog data.
    refObjLoader : `lsst.meas.algorithms.LoadReferenceObjectsTask` or `None`
        For loading reference objects
    **kwargs :
        Other keyword arguments required for instantiating a Task (e.g., 'config')
    """
    ConfigClass = DirectMatchConfig
    _DefaultName = "directMatch"

    def __init__(self, butler=None, refObjLoader=None, **kwargs):
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
        """Load reference objects and match to them.

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog`
            Catalog to match.
        filterName : `str`
            Name of filter loading fluxes
        epoch : `astropy.time.Time` or `None`
            Epoch to which to correct proper motion and parallax,
            or `None` to not apply such corrections.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - matches : Matched sources with associated reference
              (`lsst.afw.table.SourceMatchVector`)
            - matchMeta : Match metadata (`lsst.meas.astrom.MatchMetadata`)
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
        """Calculate a circle enclosing the catalog

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog`
            Catalog to encircle

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - center : ICRS center coordinate (`lsst.afw.geom.SpherePoint`)
            - radius : Radius of the circle (`lsst.geom.Angle`)
        """
        coordList = [src.getCoord() for src in catalog]
        center = averageSpherePoint(coordList)
        radius = max(center.separation(coord) for coord in coordList)
        return Struct(center=center, radius=radius + self.config.matchRadius*arcseconds)
