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

__all__ = ["FitAffineWcsTask", "FitAffineWcsConfig", "TransformedSkyWcsMaker"]


import astshim
import numpy as np
from scipy.linalg import lstsq

from lsst.afw.geom import makeSkyWcs, SkyWcs
import lsst.afw.math
from lsst.geom import Point2D, radians
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod

from ._measAstromLib import makeMatchStatisticsInRadians
from .setMatchDistance import setMatchDistance


# Keeping this around for now in case any of the fit parameters need to be
# configurable. Likely the maximum allowed shift magnitude (parameter 2 in the
# fit.)
class FitAffineWcsConfig(pexConfig.Config):
    """Config for FitTanSipWcsTask."""
    pass


class FitAffineWcsTask(pipeBase.Task):
    """Fit a TAN-SIP WCS given a list of reference object/source matches.

    This WCS fitter should be used on top of a cameraGeom distortion model as
    the model assumes that only a shift the WCS center position and a small
    affine transform are required.
    """
    ConfigClass = FitAffineWcsConfig
    _DefaultName = "fitAffineWcs"

    @timeMethod
    def fitWcs(self,
               matches,
               initWcs,
               bbox=None,
               refCat=None,
               sourceCat=None,
               exposure=None):
        """Fit a simple Affine transform with a shift to the matches and update
        the WCS.

        This method assumes that the distortion model of the telescope is
        applied correctly and is accurate with only a slight rotation,
        rotation, and "squish" required to fit to the reference locations.

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
            Ignored; present for consistency with FitSipDistortionTask.
        refCat : `lsst.afw.table.SimpleCatalog`
            reference object catalog, or None.
            If provided then all centroids are updated with the new WCS,
            otherwise only the centroids for ref objects in matches are
            updated. Required fields are "centroid_x", "centroid_y",
            "coord_ra", and "coord_dec".
        sourceCat : `lsst.afw.table.SourceCatalog`
            source catalog, or None.
            If provided then coords are updated with the new WCS;
            otherwise only the coords for sources in matches are updated.
            Required fields are "slot_Centroid_x", "slot_Centroid_y", and
            "coord_ra", and "coord_dec".
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
        # Create a data-structure that decomposes the input Wcs frames and
        # appends the new transform.
        wcsMaker = TransformedSkyWcsMaker(initWcs)

        # Grab the initial transformations going back from sky coordinates
        # and forward from pixel coordinates.
        back = wcsMaker.frameDict.getMapping(wcsMaker.frameMax, wcsMaker.frameMax-1)
        forward = wcsMaker.lastMapBeforeSky

        # Create containers for the data going into the Affine fit. This will
        # be done by approximating the solution to  Ax=b where x will be the
        # affine parameters and a linear shift. The approximate solution is
        # calculated using a least squares minimization of the Ax=b equation.
        #
        # This is looking to find the affine transform of the following form:
        # [x', y'] = [[a, b], [c, d]] [x, y] + [s, t]
        #
        # where a,b,c,d are the parameters of the affine transform, and s,t
        # are linear shift parameters.
        #
        # To solve for these unknown parameters the unknown matrix x in the
        # equation Ax=b will be of the form:
        # [a, b, c, d, s, t].
        #
        # This implies that each constraining point will correspond to two rows
        # in the A matrix in the following form:
        # [x_i, y_i, 0, 0, 1, 0]
        # [0, 0, x_i, y_i, 0, 1].
        #
        # The corresponding output points in the b vector will have the form:
        # [x'_i, y'_i, x'_(i+i), y'_(i+1)....]
        A = np.zeros((len(matches)*2, 6), dtype=float)
        b = np.empty(len(matches)*2, dtype=float)

        # Constant terms related to the shift in x and and y parameters.
        A[::2, 4] = 1
        A[1::2, 5] = 1

        # loop over each of the matches and populate the matrices.
        for i, match in enumerate(matches):
            refCoord = match.first.getCoord()
            b[i*2:i*2+2] = back.applyForward(refCoord)

            srcCentroid = match.second.getCentroid()
            val = forward.applyForward(srcCentroid)
            A[i*2, :2] = val
            A[i*2+1, 2:4] = val

        # solve for the affine and shift parameters
        # The lapack_driver parameter is set to the quickest routine tested for
        # this application at the time of writing.
        fit = lstsq(A, b, lapack_driver='gelsy')[0]

        self.log.debug("Linear shift in x: %.3f, y: %.3f, "
                       "Affine matrix: [[%.6f, %.6f], [%.6f, %.6f]]...",
                       fit[4], fit[5],
                       fit[0], fit[1], fit[2], fit[3])

        # create the final wcs
        wcs = wcsMaker.makeWcs(fit[4:], fit[:4].reshape((2, 2)))

        # Copied from other fit*WcsTasks.
        if refCat is not None:
            self.log.debug("Updating centroids in refCat")
            lsst.afw.table.updateRefCentroids(wcs, refList=refCat)
        else:
            self.log.warning("Updating reference object centroids in match list; refCat is None")
            lsst.afw.table.updateRefCentroids(
                wcs,
                refList=[match.first for match in matches])

        if sourceCat is not None:
            self.log.debug("Updating coords in sourceCat")
            lsst.afw.table.updateSourceCoords(wcs, sourceList=sourceCat)
        else:
            self.log.warning("Updating source coords in match list; sourceCat is None")
            lsst.afw.table.updateSourceCoords(
                wcs,
                sourceList=[match.second for match in matches])
        setMatchDistance(matches)

        stats = makeMatchStatisticsInRadians(wcs,
                                             matches,
                                             lsst.afw.math.MEDIAN)
        scatterOnSky = stats.getValue() * radians

        self.log.debug("In fitter scatter %.4f", scatterOnSky.asArcseconds())

        return lsst.pipe.base.Struct(
            wcs=wcs,
            scatterOnSky=scatterOnSky,
        )


class TransformedSkyWcsMaker():
    """Convenience class for appending a shifting an input SkyWcs on sky and
    appending an affine transform.

    The class assumes that all frames are sequential and are mapped one to the
    next.

    Parameters
    ----------
    input_sky_wcs : `lsst.afw.geom.SkyWcs`
        WCS to decompose and append affine matrix and shift in on sky
        location to.
    """

    def __init__(self, inputSkyWcs):
        self.frameDict = inputSkyWcs.getFrameDict()

        # Grab the order of the frames by index.
        # TODO: DM-20825
        #    Change the frame the transform is appended to to be explicitly
        #    the FIELD_ANGLE->IWC transform. Requires related tickets to be
        #    completed.
        domains = self.frameDict.getAllDomains()
        self.frameIdxs = np.sort([self.frameDict.getIndex(domain)
                                  for domain in domains])
        self.frameMin = np.min(self.frameIdxs)
        self.frameMax = np.max(self.frameIdxs)

        # Find frame just before the final mapping to sky and store those
        # indices and mappings for later.
        self.mapFrom = self.frameMax - 2
        if self.mapFrom < self.frameMin:
            self.mapFrom = self.frameMin
        self.mapTo = self.frameMax - 1
        if self.mapTo <= self.mapFrom:
            self.mapTo = self.frameMax
        self.lastMapBeforeSky = self.frameDict.getMapping(
            self.mapFrom, self.mapTo)

        # Get the original WCS sky location.

        self.origin = inputSkyWcs.getSkyOrigin()

    def makeWcs(self, linearShift, affMatrix):
        """Apply a shift and affine transform to the WCS internal to this
        class.

        A new SkyWcs with these transforms applied is returns.

        Parameters
        ----------
        linearShift : `numpy.ndarray`, (2,)
            A linear shift to apply at the same time as applying the affine
            matrix transform.
        aff_matrix : 'numpy.ndarray', (3, 3)
            Affine matrix to apply to the mapping/transform to add to the
            WCS.

        Returns
        -------
        outputWcs : `lsst.afw.geom.SkyWcs`
            Wcs with a final shift and affine transform applied.
        """
        # Create a WCS that only maps from IWC to Sky with the shifted
        # Sky origin position. This is simply the final undistorted tangent
        # plane to sky. The PIXELS to SKY map will be become our IWC to SKY
        # map and gives us our final shift position.
        iwcsToSkyWcs = makeSkyWcs(
            Point2D(0., 0.),
            self.origin,
            np.array([[1., 0.], [0., 1.]]))
        iwcToSkyMap = iwcsToSkyWcs.getFrameDict().getMapping("PIXELS", "SKY")

        # Append a simple affine Matrix transform to the current to the
        # second to last frame mapping. e.g. the one just before IWC to SKY.
        newMapping = self.lastMapBeforeSky.then(astshim.MatrixMap(affMatrix))
        newMapping = newMapping.then(astshim.ShiftMap(linearShift))

        # Create a new frame dict starting from the input_sky_wcs's first
        # frame. Append the correct mapping created above and our new on
        # sky location.
        outputFrameDict = astshim.FrameDict(
            self.frameDict.getFrame(self.frameMin))
        for frameIdx in self.frameIdxs:
            if frameIdx == self.mapFrom:
                outputFrameDict.addFrame(
                    self.mapFrom,
                    newMapping,
                    self.frameDict.getFrame(self.mapTo))
            elif frameIdx >= self.mapTo:
                continue
            else:
                outputFrameDict.addFrame(
                    frameIdx,
                    self.frameDict.getMapping(frameIdx, frameIdx + 1),
                    self.frameDict.getFrame(frameIdx + 1))
        # Append the final sky frame to the frame dict.
        outputFrameDict.addFrame(
            self.frameMax - 1,
            iwcToSkyMap,
            iwcsToSkyWcs.getFrameDict().getFrame("SKY"))

        return SkyWcs(outputFrameDict)
