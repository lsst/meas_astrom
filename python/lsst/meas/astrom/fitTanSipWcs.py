# This file is part of meas_astrom.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["FitTanSipWcsTask", "FitTanSipWcsConfig"]


import numpy as np

import lsst.geom
import lsst.sphgeom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod
from .setMatchDistance import setMatchDistance
from .sip import makeCreateWcsWithSip


class FitTanSipWcsConfig(pexConfig.Config):
    """Config for FitTanSipWcsTask."""
    order = pexConfig.RangeField(
        doc="order of SIP polynomial",
        dtype=int,
        default=2,
        min=0,
    )
    numIter = pexConfig.RangeField(
        doc="number of iterations of fitter (which fits X and Y separately, and so benefits from "
        "a few iterations",
        dtype=int,
        default=3,
        min=1,
    )
    numRejIter = pexConfig.RangeField(
        doc="number of rejection iterations",
        dtype=int,
        default=1,
        min=0,
    )
    rejSigma = pexConfig.RangeField(
        doc="Number of standard deviations for clipping level",
        dtype=float,
        default=3.0,
        min=0.0,
    )
    maxScatterArcsec = pexConfig.RangeField(
        doc="maximum median scatter of a WCS fit beyond which the fit fails (arcsec); "
        "be generous, as this is only intended to catch catastrophic failures",
        dtype=float,
        default=10,
        min=0,
    )


class FitTanSipWcsTask(pipeBase.Task):
    """Fit a TAN-SIP WCS given a list of reference object/source matches.
    """
    ConfigClass = FitTanSipWcsConfig
    _DefaultName = "fitWcs"

    @timeMethod
    def fitWcs(self, matches, initWcs, bbox=None, refCat=None, sourceCat=None, exposure=None):
        """Fit a TAN-SIP WCS from a list of reference object/source matches

        Parameters
        ----------
        matches : `list` of `lsst.afw.table.ReferenceMatch`
            The following fields are read:

            - match.first (reference object) coord
            - match.second (source) centroid

            The following fields are written:

            - match.first (reference object) centroid,
            - match.second (source) centroid
            - match.distance (on sky separation, in radians)

        initWcs : `lsst.afw.geom.SkyWcs`
            initial WCS
        bbox : `lsst.geom.Box2I`
            the region over which the WCS will be valid (an lsst:afw::geom::Box2I);
            if None or an empty box then computed from matches
        refCat : `lsst.afw.table.SimpleCatalog`
            reference object catalog, or None.
            If provided then all centroids are updated with the new WCS,
            otherwise only the centroids for ref objects in matches are updated.
            Required fields are "centroid_x", "centroid_y", "coord_ra", and "coord_dec".
        sourceCat : `lsst.afw.table.SourceCatalog`
            source catalog, or None.
            If provided then coords are updated with the new WCS;
            otherwise only the coords for sources in matches are updated.
            Required fields are "slot_Centroid_x", "slot_Centroid_y", and "coord_ra", and "coord_dec".
        exposure : `lsst.afw.image.Exposure`
            Ignored; present for consistency with FitSipDistortionTask.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            with the following fields:

            - ``wcs`` :  the fit WCS (`lsst.afw.geom.SkyWcs`)
            - ``scatterOnSky`` :  median on-sky separation between reference
              objects and sources in "matches" (`lsst.afw.geom.Angle`)
        """
        if bbox is None:
            bbox = lsst.geom.Box2I()

        import lsstDebug
        debug = lsstDebug.Info(__name__)

        wcs = self.initialWcs(matches, initWcs)
        rejected = np.zeros(len(matches), dtype=bool)
        for rej in range(self.config.numRejIter):
            sipObject = self._fitWcs([mm for i, mm in enumerate(matches) if not rejected[i]], wcs)
            wcs = sipObject.getNewWcs()
            rejected = self.rejectMatches(matches, wcs, rejected)
            if rejected.sum() == len(rejected):
                raise RuntimeError("All matches rejected in iteration %d" % (rej + 1,))
            self.log.debug(
                "Iteration {0} of astrometry fitting: rejected {1} outliers, "
                "out of {2} total matches.".format(
                    rej, rejected.sum(), len(rejected)
                )
            )
            if debug.plot:
                print("Plotting fit after rejection iteration %d/%d" % (rej + 1, self.config.numRejIter))
                self.plotFit(matches, wcs, rejected)
        # Final fit after rejection
        sipObject = self._fitWcs([mm for i, mm in enumerate(matches) if not rejected[i]], wcs)
        wcs = sipObject.getNewWcs()
        if debug.plot:
            print("Plotting final fit")
            self.plotFit(matches, wcs, rejected)

        if refCat is not None:
            self.log.debug("Updating centroids in refCat")
            afwTable.updateRefCentroids(wcs, refList=refCat)
        else:
            self.log.warn("Updating reference object centroids in match list; refCat is None")
            afwTable.updateRefCentroids(wcs, refList=[match.first for match in matches])

        if sourceCat is not None:
            self.log.debug("Updating coords in sourceCat")
            afwTable.updateSourceCoords(wcs, sourceList=sourceCat)
        else:
            self.log.warn("Updating source coords in match list; sourceCat is None")
            afwTable.updateSourceCoords(wcs, sourceList=[match.second for match in matches])

        self.log.debug("Updating distance in match list")
        setMatchDistance(matches)

        scatterOnSky = sipObject.getScatterOnSky()

        if scatterOnSky.asArcseconds() > self.config.maxScatterArcsec:
            raise pipeBase.TaskError(
                "Fit failed: median scatter on sky = %0.3f arcsec > %0.3f config.maxScatterArcsec" %
                (scatterOnSky.asArcseconds(), self.config.maxScatterArcsec))

        return pipeBase.Struct(
            wcs=wcs,
            scatterOnSky=scatterOnSky,
        )

    def initialWcs(self, matches, wcs):
        """Generate a guess Wcs from the astrometric matches

        We create a Wcs anchored at the center of the matches, with the scale
        of the input Wcs.  This is necessary because matching returns only
        matches with no estimated Wcs, and the input Wcs is a wild guess.
        We're using the best of each: positions from the matches, and scale
        from the input Wcs.

        Parameters
        ----------
        matches : `list` of `lsst.afw.table.ReferenceMatch`
            List of sources matched to references.
        wcs : `lsst.afw.geom.SkyWcs`
            Current WCS.

        Returns
        -------
        newWcs : `lsst.afw.geom.SkyWcs`
            Initial WCS guess from estimated crpix and crval.
        """
        crpix = lsst.geom.Extent2D(0, 0)
        crval = lsst.sphgeom.Vector3d(0, 0, 0)
        for mm in matches:
            crpix += lsst.geom.Extent2D(mm.second.getCentroid())
            crval += mm.first.getCoord().getVector()
        crpix /= len(matches)
        crval /= len(matches)
        newWcs = afwGeom.makeSkyWcs(crpix=lsst.geom.Point2D(crpix),
                                    crval=lsst.geom.SpherePoint(crval),
                                    cdMatrix=wcs.getCdMatrix())
        return newWcs

    def _fitWcs(self, matches, wcs):
        """Fit a Wcs based on the matches and a guess Wcs.

        Parameters
        ----------
        matches : `list` of `lsst.afw.table.ReferenceMatch`
            List of sources matched to references.
        wcs : `lsst.afw.geom.SkyWcs`
            Current WCS.

        Returns
        -------
        sipObject : `lsst.meas.astrom.sip.CreateWcsWithSip`
            Fitted SIP object.
        """
        for i in range(self.config.numIter):
            sipObject = makeCreateWcsWithSip(matches, wcs, self.config.order)
            wcs = sipObject.getNewWcs()
        return sipObject

    def rejectMatches(self, matches, wcs, rejected):
        """Flag deviant matches

        We return a boolean numpy array indicating whether the corresponding
        match should be rejected.  The previous list of rejections is used
        so we can calculate uncontaminated statistics.

        Parameters
        ----------
        matches : `list` of `lsst.afw.table.ReferenceMatch`
            List of sources matched to references.
        wcs : `lsst.afw.geom.SkyWcs`
            Fitted WCS.
        rejected : array-like of `bool`
            Array of matches rejected from the fit. Unused.

        Returns
        -------
        rejectedMatches : `ndarray` of type `bool`
            Matched objects found to be outside of tolerance.
        """
        fit = [wcs.skyToPixel(m.first.getCoord()) for m in matches]
        dx = np.array([ff.getX() - mm.second.getCentroid().getX() for ff, mm in zip(fit, matches)])
        dy = np.array([ff.getY() - mm.second.getCentroid().getY() for ff, mm in zip(fit, matches)])
        good = np.logical_not(rejected)
        return (dx > self.config.rejSigma*dx[good].std()) | (dy > self.config.rejSigma*dy[good].std())

    def plotFit(self, matches, wcs, rejected):
        """Plot the fit

        We create four plots, for all combinations of (dx, dy) against
        (x, y).  Good points are black, while rejected points are red.

        Parameters
        ----------
        matches : `list` of `lsst.afw.table.ReferenceMatch`
            List of sources matched to references.
        wcs : `lsst.afw.geom.SkyWcs`
            Fitted WCS.
        rejected : array-like of `bool`
            Array of matches rejected from the fit.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as e:
            self.log.warn("Unable to import matplotlib: %s", e)
            return

        fit = [wcs.skyToPixel(m.first.getCoord()) for m in matches]
        x1 = np.array([ff.getX() for ff in fit])
        y1 = np.array([ff.getY() for ff in fit])
        x2 = np.array([m.second.getCentroid().getX() for m in matches])
        y2 = np.array([m.second.getCentroid().getY() for m in matches])

        dx = x1 - x2
        dy = y1 - y2

        good = np.logical_not(rejected)

        figure = plt.figure()
        axes = figure.add_subplot(2, 2, 1)
        axes.plot(x2[good], dx[good], 'ko')
        axes.plot(x2[rejected], dx[rejected], 'ro')
        axes.set_xlabel("x")
        axes.set_ylabel("dx")

        axes = figure.add_subplot(2, 2, 2)
        axes.plot(x2[good], dy[good], 'ko')
        axes.plot(x2[rejected], dy[rejected], 'ro')
        axes.set_xlabel("x")
        axes.set_ylabel("dy")

        axes = figure.add_subplot(2, 2, 3)
        axes.plot(y2[good], dx[good], 'ko')
        axes.plot(y2[rejected], dx[rejected], 'ro')
        axes.set_xlabel("y")
        axes.set_ylabel("dx")

        axes = figure.add_subplot(2, 2, 4)
        axes.plot(y2[good], dy[good], 'ko')
        axes.plot(y2[rejected], dy[rejected], 'ro')
        axes.set_xlabel("y")
        axes.set_ylabel("dy")

        plt.show()
