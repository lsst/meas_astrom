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

__all__ = ["FitAffineWcsTask", "FitAffineWcsConfig"]


import astshim
import numpy as np
from scipy.optimize import least_squares, minimize, Bounds

from lsst.afw.geom import (makeSkyWcs, degrees, arcseconds, arcminutes, radians, SkyWcs, SpherePoint) 
import lsst.afw.math
from lsst.geom import Point2D
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.sphgeom as sphgeom

from .makeMatchStatistics import makeMatchStatisticsInRadians


def _chi_func(x, ref_points, src_pixels, wcs_maker):
    """Function to minimize to fit the shift and rotation in the WCS.

    Parameters
    ----------
    ref_points : `list` of `lsst.afw.geom.SpherePoint`
        Reference object on Sky locations.
    src_pixels : `list` of `lsst.geom.Point2D`
        Source object positions on the pixels.
    wcs_maker : `TransformedSkyWcsMaker`
        Container class for producing the updated Wcs.

    Returns
    -------
    output_separations : `list` of `float`
        Separation between predicted source location and reference location in
        radians.
    """
    wcs = wcs_maker.makeWcs(x[:2], x[2:].reshape((2, 2)))

    output_separations = []
    # Fit both sky to pixel and pixel to sky to avoid any non-invertible
    # affine matrixes.
    for ref, src in zip(ref_points, src_pixels):
        sky_sep = ref.getTangentPlaneOffset(wcs.pixelToSky(src))
        output_separations.append(sky_sep[0].asArcseconds())
        output_separations.append(sky_sep[1].asArcseconds())
        xy_sep = src - wcs.skyToPixel(ref)
        # Convert the pixel separations to units, arcseconds to match units
        # of sky separation.
        output_separations.append(xy_sep[0] * wcs.getPixelScale(src).asArcseconds())
        output_separations.append(xy_sep[1] * wcs.getPixelScale(src).asArcseconds())

    return output_separations


# Keeping this around for now in case any of the fit parameters need to be
# configurable. Likely the maximum allowed shift magnitude (parameter 2 in the
# fit.)
class FitAffineWcsConfig(pexConfig.Config):
    """Config for FitTanSipWcsTask."""
    pass


class FitAffineWcsTask(pipeBase.Task):
    """Fit a TAN-SIP WCS given a list of reference object/source matches.
    """
    ConfigClass = FitAffineWcsConfig
    _DefaultName = "fitAffineWcs"

    @pipeBase.timeMethod
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
        wcs_maker = TransformedSkyWcsMaker(initWcs)

        ref_points = []
        src_pixels = []
        offset_dir = 0
        offset_dist = 0
        # Grab reference coordinates and source centroids. Compute the average
        # direction and separation between the reference and the sources.
        # I'm not sure if bearingTo should be computed from the src to ref
        # or ref to source.
        for match in matches:
            ref_coord = match.first.getCoord()
            ref_points.append(ref_coord)
            src_centroid = match.second.getCentroid()
            src_pixels.append(src_centroid)
            src_coord = initWcs.pixelToSky(src_centroid)
            offset_dir += src_coord.bearingTo(ref_coord).asDegrees()
            offset_dist += src_coord.separation(ref_coord).asArcminutes()
        offset_dir /= len(src_pixels)
        offset_dist /= len(src_pixels)
        if offset_dir > 180:
            offset_dir = offset_dir - 360

        # Best performing fitter in scipy tried so far (vs. default settings in
        # minimize). Fits all current test cases with a scatter of a most 0.15
        # arcseconds. exits early because of the xTol value which cannot be
        # disabled in scipy1.2.1.
        fit = least_squares(_chi_func,
                            x0=[offset_dir, offset_dist, 1., 0., 0., 1.],
                            args=(ref_points,
                                  src_pixels,
                                  wcs_maker),
                            method='dogbox',
                            bounds=[[-180, 0, -np.inf, -np.inf, -np.inf, -np.inf],
                                    [180, np.inf, np.inf, np.inf, np.inf, np.inf]],
                            ftol=2.3e-16,
                            gtol=2.31e-16,
                            xtol=2.3e-16)

        wcs = wcs_maker.makeWcs(fit.x[:2], fit.x[2:].reshape((2, 2)))

        # Copied from other fit*WcsTasks.
        if refCat is not None:
            self.log.debug("Updating centroids in refCat")
            lsst.afw.table.updateRefCentroids(wcs, refList=refCat)
        else:
            self.log.warn("Updating reference object centroids in match list; "
                          "refCat is None")
            lsst.afw.table.updateRefCentroids(
                wcs,
                refList=[match.first for match in matches])

        if sourceCat is not None:
            self.log.debug("Updating coords in sourceCat")
            lsst.afw.table.updateSourceCoords(wcs, sourceList=sourceCat)
        else:
            self.log.warn("Updating source coords in match list; sourceCat is "
                          "None")
            lsst.afw.table.updateSourceCoords(
                wcs,
                sourceList=[match.second for match in matches])

        stats = makeMatchStatisticsInRadians(wcs,
                                             matches,
                                             lsst.afw.math.MEDIAN)
        scatterOnSky = stats.getValue() * radians

        return lsst.pipe.base.Struct(
            wcs=wcs,
            scatterOnSky=scatterOnSky,
        )


class TransformedSkyWcsMaker(object):
    """Container class for appending a shift/rotation to an input SkyWcs.

    The class assumes that all frames are sequential and are mapped one to the
    next.

    Parameters
    ----------
    input_sky_wcs : `lsst.afw.geom.SkyWcs`
        WCS to decompose and append rotation matrix and shift in on sky
        location to.
    """

    def __init__(self, input_sky_wcs):
        self.frame_dict = input_sky_wcs.getFrameDict()

        # Grab the order of the frames by index.
        domains = self.frame_dict.getAllDomains()
        self.frame_idxs = np.sort([self.frame_dict.getIndex(domain)
                                   for domain in domains])
        self.frame_min = np.min(self.frame_idxs)
        self.frame_max = np.max(self.frame_idxs)

        # Find frame just before the final mapping to sky and store those
        # indices and mappings for later.
        self.map_from = self.frame_max - 2
        if self.map_from < self.frame_min:
            self.map_from = self.frame_min
        self.map_to = self.frame_max - 1
        if self.map_to <= self.map_from:
            self.map_to = self.frame_max
        self.last_map_before_sky = self.frame_dict.getMapping(
            self.map_from, self.map_to)

        # Get the original WCS sky location.

        self.origin = input_sky_wcs.getSkyOrigin()

    def makeWcs(self, crval_offset, rot_matrix):
        """Apply a shift and rotation to the WCS internal to this class,
        a new Wcs with these transforms applied.

        Parameters
        ----------
        crval_shift : `numpy.ndarray`, (2,)
            Shift in radians to apply to the Wcs origin/crvals.
        rot_matrix : 'numpy.ndarray', (3, 3)
            Rotation matrix to apply to the mapping/transform to add to the
            WCS. Rotation matrix need not be unitary and can therefor
            stretch/squish as well as rotate.

        Returns
        -------
        output_wcs : `lsst.afw.geom.SkyWcs`
            Wcs with a final shift and rotation applied.
        """
        # Create a WCS that only maps from IWC to Sky with the shifted
        # Sky origin position. This is simply the final undistorted tangent
        # plane to sky. The PIXELS to SKY map will be become our IWC to SKY
        # map and gives us our final shift position.
        iwcs_to_sky_wcs = makeSkyWcs(
            Point2D(0., 0.),
            self.origin.offset(crval_offset[0] * degrees,
                               crval_offset[1] * arcminutes),
            np.array([[1., 0.], [0., 1.]]))
        iwc_to_sky_map = iwcs_to_sky_wcs.getFrameDict().getMapping("PIXELS",
                                                                   "SKY")

        # Append a simple rotation Matrix transform to the current to the
        # second to last frame mapping. e.g. the one just before IWC to SKY.
        new_mapping = self.last_map_before_sky.then(
            astshim.MatrixMap(rot_matrix))

        # Create a new frame dict starting from the input_sky_wcs's first
        # frame. Append the correct mapping created above and our new on
        # sky location.
        output_frame_dict = astshim.FrameDict(
            self.frame_dict.getFrame(self.frame_min))
        for frame_idx in self.frame_idxs:
            if frame_idx == self.map_from:
                output_frame_dict.addFrame(
                    self.map_from,
                    new_mapping,
                    self.frame_dict.getFrame(self.map_to))
            elif frame_idx >= self.map_to:
                continue
            else:
                output_frame_dict.addFrame(
                    frame_idx,
                    self.frame_dict.getMapping(frame_idx, frame_idx + 1),
                    self.frame_dict.getFrame(frame_idx + 1))
        # Append the final sky frame to the frame dict.
        output_frame_dict.addFrame(
            self.frame_max - 1,
            iwc_to_sky_map,
            iwcs_to_sky_wcs.getFrameDict().getFrame("SKY"))

        return SkyWcs(output_frame_dict)
