#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
from lsst.afw.table import SourceCatalog
from lsst.meas.algorithms import BaseStarSelectorTask, starSelectorRegistry
from lsst.pipe.base import Struct
import lsst.pex.config as pexConfig
import lsst.afw.display.ds9 as ds9


class CatalogStarSelectorConfig(BaseStarSelectorTask.ConfigClass):
    fluxLim = pexConfig.RangeField(
        doc="specify the minimum psfFlux for good Psf Candidates",
        dtype=float,
        default=0.0,
        min=0.0,
    )
    fluxMax = pexConfig.RangeField(
        doc="specify the maximum psfFlux for good Psf Candidates (ignored if == 0)",
        dtype=float,
        default=0.0,
        min=0.0,
    )

    def setDefaults(self):
        BaseStarSelectorTask.ConfigClass.setDefaults(self)
        self.badFlags = [
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_saturatedCenter",
        ]


class CheckSource(object):
    """A functor to check whether a source has any flags set that should cause it to be labeled bad."""

    def __init__(self, table, fluxLim, fluxMax, badFlags):
        self.keys = [table.getSchema().find(name).key for name in badFlags]
        self.keys.append(table.getCentroidFlagKey())
        self.fluxLim = fluxLim
        self.fluxMax = fluxMax

    def __call__(self, source):
        for k in self.keys:
            if source.get(k):
                return False
        if self.fluxLim is not None and source.getPsfFlux() < self.fluxLim:  # ignore faint objects
            return False
        if self.fluxMax != 0.0 and source.getPsfFlux() > self.fluxMax:  # ignore bright objects
            return False
        return True

## \addtogroup LSST_task_documentation
## \{
## \page CatalogStarSelectorTask
## \ref CatalogStarSelectorTask_ "CatalogStarSelectorTask"
## \copybrief CatalogStarSelectorTask
## \}


class CatalogStarSelectorTask(object):
    """!Select stars based on a reference catalog

    @anchor CatalogStarSelectorTask_

    @section meas_astrom_catalogStarSelector_Contents  Contents

     - @ref meas_astrom_catalogStarSelector_Purpose
     - @ref meas_astrom_catalogStarSelector_Initialize
     - @ref meas_astrom_catalogStarSelector_IO
     - @ref meas_astrom_catalogStarSelector_Config
     - @ref meas_astrom_catalogStarSelector_Debug

    @section meas_astrom_catalogStarSelector_Purpose  Description

    Select stars using a match list: select sources where the matching reference object is unresolved,
    plus the source passes the following tests:
    - no flag from config.badFlags is set
    - psf flux >= config.fluxLim
    - psf flux <= config.fluxMax (not checked if fluxMax == 0)

    @section meas_astrom_catalogStarSelector_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_catalogStarSelector_IO  Invoking the Task

    Like all star selectors, the main method is `run`. Unlike most star selectors,
    this one requires the `matches` argument (the `usesMatches` property is true).

    @section meas_astrom_catalogStarSelector_Config  Configuration parameters

    See @ref CatalogStarSelectorConfig

    @section meas_astrom_catalogStarSelector_Debug  Debug variables

    CatalogStarSelectorTask has a debug dictionary with the following keys:
    <dl>
    <dt>display
    <dd>bool; if True display debug information
    <dt>pauseAtEnd
    <dd>bool; if True wait after displaying everything and wait for user input
    </dl>

    For example, put something like:
    @code{.py}
        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)  # N.b. lsstDebug.Info(name) would call us recursively
            if name.endswith("catalogStarSelector"):
                di.display = True

            return di

        lsstDebug.Info = DebugInfo
    @endcode
    into your `debug.py` file and run your task with the `--debug` flag.
    """
    ConfigClass = CatalogStarSelectorConfig
    usesMatches = True  # `run` and `selectStars` require the `matches` argument

    def selectStars(self, exposure, sourceCat, matches=None):
        """!Return a list of PSF candidates that represent likely stars

        A list of PSF candidates may be used by a PSF fitter to construct a PSF.

        @param[in] exposure  the exposure containing the sources
        @param[in] sourceCat  catalog of sources that may be stars (an lsst.afw.table.SourceCatalog)
        @param[in] matches  a match vector as produced by meas_astrom; required
                            (defaults to None to match the StarSelector API and improve error handling)

        @return an lsst.pipe.base.Struct containing:
        - starCat  catalog of selected stars (a subset of sourceCat)
        """
        import lsstDebug
        debugInfo = lsstDebug.Info(__name__)
        display = debugInfo.display
        pauseAtEnd = debugInfo.pauseAtEnd  # pause when done

        if matches is None:
            raise RuntimeError("CatalogStarSelectorTask requires matches")

        mi = exposure.getMaskedImage()

        if display:
            frame = 1
            ds9.mtv(mi, frame=frame, title="PSF candidates")

        isGoodSource = CheckSource(sourceCat, self.config.fluxLim, self.config.fluxMax, self.config.badFlags)

        starCat = SourceCatalog(sourceCat.schema)
        with ds9.Buffering():
            for ref, source, d in matches:
                if not ref.get("resolved"):
                    if not isGoodSource(source):
                        symb, ctype = "+", ds9.RED
                    else:
                        starCat.append(source)
                        symb, ctype = "+", ds9.GREEN

                        if display:
                            ds9.dot(symb, source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                                    size=4, frame=frame, ctype=ctype)

        if display and pauseAtEnd:
            raw_input("Continue? y[es] p[db] ")

        return Struct(
            starCat=starCat,
        )

starSelectorRegistry.register("catalog", CatalogStarSelectorTask)
