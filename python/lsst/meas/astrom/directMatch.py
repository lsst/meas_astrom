from __future__ import absolute_import, division, print_function

from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.base import Task, Struct
from .createMatchMetadata import MatchMetadata
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
from lsst.afw.geom import arcseconds


__all__ = ["DirectMatchConfig", "DirectMatchTask"]


class DirectMatchConfig(Config):
    """Configuration for DirectMatchTask"""
    refObjLoader = ConfigurableField(target=LoadIndexedReferenceObjectsTask, doc="Load reference objects")
    matchRadius = Field(dtype=float, default=0.25, doc="Matching radius, arcsec")


class DirectMatchTask(Task):
    """!Simple matching of a source catalog to a reference catalog

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
            self.makeSubtask("refObjLoader", butler=butler)
        else:
            self.refObjLoader = refObjLoader

    def run(self, catalog, filterName=None):
        """!Load reference objects and match to them

        @param[in] catalog  Catalog to match to (lsst.afw.table.SourceCatalog)
        @param[in] filterName  Name of filter, for loading fluxes (str)
        @return Struct with matches (lsst.afw.table.SourceMatchVector) and
            matchMeta (lsst.meas.astrom.MatchMetadata)
        """
        circle = self.calculateCircle(catalog)
        matchMeta = MatchMetadata(circle.center, circle.radius, filterName)
        refData = self.refObjLoader.loadSkyCircle(circle.center, circle.radius, filterName)
        matches = afwTable.matchRaDec(refData.refCat, catalog, self.config.matchRadius*arcseconds)
        self.log.info("Matched %d from %d input and %d reference sources" %
                      (len(matches), len(catalog), len(refData.refCat)))
        return Struct(matches=matches, matchMeta=matchMeta)

    def calculateCircle(self, catalog):
        """!Calculate a circle enclosing the catalog

        @param[in] catalog  Catalog we will encircle (lsst.afw.table.SourceCatalog)
        @return Struct with center (lsst.afw.coord.Coord) and radius (lsst.afw.geom.Angle)
        """
        coordList = [src.getCoord() for src in catalog]
        center = afwCoord.averageCoord(coordList)
        radius = max(center.angularSeparation(coord) for coord in coordList)
        return Struct(center=center, radius=radius + self.config.matchRadius*arcseconds)
