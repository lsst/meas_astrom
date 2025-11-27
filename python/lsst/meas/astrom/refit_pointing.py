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

from __future__ import annotations

__all__ = ("RefitPointingConfig", "RefitPointingTask", "NoVisitWcs")

import math

import numpy as np

from lsst.afw.geom import getIntermediateWorldCoordsToSky, SipApproximation, makeTanSipWcs
from lsst.afw.cameraGeom import FIELD_ANGLE, PIXELS
from lsst.geom import (
    Angle,
    Box2D,
    Extent2I,
    Point2D,
    SpherePoint,
    SphereTransform,
    degrees,
    radians,
)
from lsst.pex.config import Config, Field
from lsst.pipe.base import AlgorithmError, Struct, Task
from lsst.obs.base.visit_geometry import VisitGeometry
from lsst.obs.base.utils import createInitialSkyWcsFromBoresight
from lsst.sphgeom import ConvexPolygon


class NoVisitWcs(AlgorithmError):
    """Exception raised when there are no WCSs for any detectors in a visit."""

    @property
    def metadata(self):
        return {}


class RefitPointingConfig(Config):
    pointing_grid_spacing = Field[float](
        doc=(
            "Spacing (in pixels) between grid points used to evaluate the WCS when fitting the pointing. "
            "This can be a very sparse grid (there are only three degrees of freedom).  "
            "If the spacing does not divide the detector bounding box evenly, it is decreased slightly."
        ),
        dtype=float,
        default=512.0,
    )
    sip_grid_spacing = Field[float](
        doc=(
            "Spacing (in pixels) between grid points used to evaluate the WCS when fitting the "
            "SIP approximation.  This needs to be a fine grid to fit the number of degrees of freedom in "
            "the configured polynomial order. "
            "If the spacing does not divide the detector bounding box evenly, it is decreased slightly."
        ),
        dtype=float,
        default=32.0,
    )
    sip_order = Field[int](
        doc="Polynomial order for the SIP approximation.",
        dtype=int,
        default=5,
    )
    add_wcs_fallbacks = Field[bool](
        doc=(
            "If True, add a fallback WCS (a raw-like camera geometry one, with an updated pointing) "
            "for any detector in the given summary catalog that does not already have a WCS.",
        ),
        dtype=bool,
        default=False,
    )


class RefitPointingTask(Task):
    """A task that uses the available WCSs of the detectors in a visit to
    re-fit the pointing for that visit, add SIP approximations to all WCSs,
    and compute new visit regions for the butler.
    """

    _DefaultName = "refitPointing"
    ConfigClass = RefitPointingConfig

    def __init__(self, config=None, *, schema=None, **kwargs):
        super().__init__(config, **kwargs)
        self._flag_key = None
        self._delta1_key = None
        self._delta2_key = None
        if schema is not None:
            if self.config.add_wcs_fallbacks:
                self._flag_key = schema.addField(
                    "wcs_is_fallback", type="Flag",
                    doc=(
                        "Whether the WCS for this detector is just a re-pointed raw WCS, "
                        "as opposed to something actually fit from stars."
                    )
                )
            self._delta1_key = schema.addField(
                "wcs_approx_delta1", type="Angle",
                doc=(
                    "Maximum of ``|target.pixelToSky(p) - approx.pixelToSky(p)|`` on a "
                    "grid offset from the one used to fit the SIP approximation."
                )
            )
            self._delta2_key = schema.addField(
                "wcs_approx_delta2", type="Angle",
                doc=(
                    "Maximum of "
                    "``|target.pixelToSky(p) - target.pixelToSky(approx.skyToPixel(target.pixelToSky(p)))|`` "
                    "on a grid offset from the one used to fit the SIP approximation."
                )
            )

    def run(self, *, catalog, camera):
        """Re-fit the pointing from the WCSs in a visit.

        Parameters
        ----------
        catalog : `lsst.afw.table.ExposureCatalog`
            A catalog of per-detector records for the visit.  WCS components
            are updated in-place: FITS-compatible SIP approximations are added,
            and re-pointed raw-like WCSs are added for detectors that had no
            fitted WCS if `~RefitPointingConfig.add_wcs_fallbacks` is `True`.
        camera : `lsst.afw.cameraGeom.Camera`
            Camera geometry.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            A struct with the following attributes:

            - boresight (`lsst.geom.SpherePoint`): new boresight location
            - orientation (`lsst.geom.Angle`): new orientation angle
            - catalog (`lsst.afw.table.ExposureCatalog`): the same catalog that
              was passed in, after modification in-place.
            - regions (`lsst.obs.base.VisitGeometry`): updated regions for the
              visit and all detectors.

        Raises
        ------
        NoVisitWcs
            Raised if ``catalog`` is empty or if there are no WCSs for any
            detectors.
        """
        if not catalog:
            raise NoVisitWcs("No detector rows in visit catalog.")
        boresight, orientation = self._fit_pointing(catalog, camera)
        if (visit_info := catalog[0].getVisitInfo()) is not None:
            old_boresight = visit_info.getBoresightRaDec()
            offset = old_boresight.separation(boresight)
            self.log.info(
                "Re-fit pointing is %s, orientation=%0.2f deg (%0.2g deg from the original boresight).",
                boresight,
                orientation.asDegrees(),
                offset.asDegrees(),
            )
        else:
            self.log.info("Re-fit pointing is %s, orientation=%0.2f deg.", boresight, orientation.asDegrees())
        max_delta1 = None
        max_delta2 = None
        for record in catalog:
            target_wcs = record.getWcs()
            detector = camera[record.getId()]
            fallback_wcs = createInitialSkyWcsFromBoresight(boresight, orientation, detector=detector)
            if target_wcs is None:
                if self.config.add_wcs_fallbacks:
                    self.log.info(
                        "Installing a re-pointed raw-like fallback WCS for detector %d.", record.getId()
                    )
                    target_wcs = fallback_wcs
                    if self._flag_key is not None:
                        record.set(self._flag_key, True)
                else:
                    continue
            sip_approx = self._fit_sip_approximation(target_wcs, fallback_wcs, detector=detector)
            fits_wcs = makeTanSipWcs(
                sip_approx.getPixelOrigin(),
                fallback_wcs.getSkyOrigin(),
                sip_approx.getCdMatrix(),
                sip_approx.getA(),
                sip_approx.getB(),
                sip_approx.getAP(),
                sip_approx.getBP(),
            )
            record.setWcs(target_wcs.copyWithFitsApproximation(fits_wcs))
            delta1, delta2 = self._validate_sip_approximation(record, detector)
            if max_delta1 is None or max_delta1 < delta1:
                max_delta1 = delta1
            if max_delta2 is None or max_delta2 < delta2:
                max_delta2 = delta2
        self.log.info(
            'SIP approximations have pixelToSky good to ≤ %0.2g mas and '
            'skyToPixel good to ≤ %0.2g mas.',
            max_delta1.asArcseconds() * 1000,
            max_delta2.asArcseconds() * 1000,
        )
        regions = self._make_visit_geometry(boresight, orientation, catalog, camera)
        return Struct(
            boresight=boresight,
            orientation=orientation,
            catalog=catalog,
            regions=regions,
        )

    def _fit_pointing(self, catalog, camera):
        """Fit the pointing for a visit from the detectors in that visit that
        have a fitted WCS.

        Parameters
        ----------
        catalog : `lsst.afw.table.ExposureCatalog`
            A catalog of per-detector records for the visit.
        camera : `lsst.afw.cameraGeom.Camera`
            Camera geometry.

        Returns
        -------
        boresight : `lsst.geom.SpherePoint`
            New boresight location.
        orientation : `lsst.geom.Angle`
            New orientation angle.
        """
        start_arrays = []
        target_arrays = []
        start_boresight: SpherePoint | None = None
        start_orientation = 0.0*degrees
        start_y_axis_point: SpherePoint | None = None
        for record in catalog:
            # We call the WCSs that were actually fit to the stars the "true"
            # WCSs.
            target_wcs = record.getWcs()
            if target_wcs is None:
                continue
            detector = camera[record.getId()]
            if start_boresight is None:
                # We just need some semi-arbitrary point on the sky that lets
                # extract the camera geometry part of a raw WCS.  Might be
                # helpful to have it in the right hemisphere, but otherwise it
                # shouldn't matter.
                start_boresight = target_wcs.pixelToSky(Point2D(0.0, 0.0))
            # Make a raw-like WCS at the arbitrary boresight and orientation.
            start_wcs = createInitialSkyWcsFromBoresight(
                start_boresight, start_orientation, detector=detector
            )
            # Make a grid of positions for the detector and map them to the sky
            # via both the true WCS and the arbitrary raw-like one, but in
            # xyz unit-vector form.
            pixel_x, pixel_y = self._make_grid(detector, self.config.pointing_grid_spacing)
            start_ra, start_dec = start_wcs.pixelToSkyArray(pixel_x, pixel_y)
            start_arrays.append(
                np.stack(
                    SpherePoint.toUnitXYZ(longitude=start_ra, latitude=start_dec, units=radians),
                    axis=1,
                )
            )
            target_ra, target_dec = target_wcs.pixelToSkyArray(pixel_x, pixel_y)
            target_arrays.append(
                np.stack(
                    SpherePoint.toUnitXYZ(longitude=target_ra, latitude=target_dec, units=radians),
                    axis=1,
                )
            )
            if start_y_axis_point is None:
                # Define a point at (0, 1deg) in the FIELD_ANGLE system,
                # according to the raw WCS, to let us fit the rotation angle.
                # We could do this with any detector and get the same answer.
                # We're going from field angle to pixels to sky even though the
                # raw WCS goes from pixels to field angle to sky, because that
                # minimizes how much this code knows about how the rotation
                # angle is defined.
                start_y_axis_point = start_wcs.pixelToSky(
                    detector.transform(Point2D(0.0, np.pi / 180.0), FIELD_ANGLE, PIXELS)
                )
        if not start_arrays or not target_arrays:
            raise NoVisitWcs("No fitted WCSs found for visit.")
        # Fit the spherical rotation that maps the points in the arbitrary
        # raw-like WCS to the true WCS.
        start_concat = np.concatenate(start_arrays)
        target_concat = np.concatenate(target_arrays)
        transform = SphereTransform.fit_unit_vectors(start_concat, target_concat)
        # If we apply that same rotation to our arbitrary start boresight, we
        # get the boresight predicted by the true WCSs.
        boresight = transform(start_boresight)
        # If we apply that rotation to our point on the FIELD_ANGLE y-axis, we
        # can similarly recover the orientation angle predicted by the true
        # WCSs.
        transformed_y_axis_point = transform(start_y_axis_point)
        orientation = Angle(90, degrees) - boresight.bearingTo(transformed_y_axis_point)
        if camera.getFocalPlaneParity():
            raise NotImplementedError("Cameras with focal plane parity flips are not yet supported.")
        return boresight, orientation

    def _fit_sip_approximation(self, target_wcs, fallback_wcs, detector):
        """Fit a SIP approximation to a WCS.

        Parameters
        ----------
        target_wcs : `lsst.afw.geom.SkyWcs`
            The WCS to approximate.  This is generally a WCS fit to stars, but
            it may just be the repointed raw WCS (in which case this method
            will effectively just fit a polynomial approximation to the camera
            geometry transforms).
        fallback_wcs : `lsst.afw.geom.SkyWcs`
            A raw-like WCS using the updated boresight and rotation angle.
        detector : `lsst.afw.cameraGeom.Detector`
            Camera geometry for this detector.

        Returns
        -------
        sip_wcs : `lsst.afw.geom.SkyWcs`
            Approximation to the given WCS.
        """
        pixels_to_iwc = target_wcs.getTransform().then(
            getIntermediateWorldCoordsToSky(fallback_wcs, simplify=True).inverted()
        )
        pixel_bbox = Box2D(detector.getBBox())
        grid_shape = Extent2I(
            # We add one to the dimensions since there's a point the min and
            # max in each dimension.
            np.ceil(pixel_bbox.width / self.config.sip_grid_spacing) + 1,
            np.ceil(pixel_bbox.height / self.config.sip_grid_spacing) + 1,
        )
        sip_approx = SipApproximation(
            pixels_to_iwc,
            fallback_wcs.getPixelOrigin(),
            fallback_wcs.getCdMatrix(),
            pixel_bbox,
            grid_shape,
            self.config.sip_order,
        )
        return sip_approx

    def _validate_sip_approximation(self, record, detector):
        """Validate the SIP approximation to a WCS.

        Parameters
        ----------
        record : `lsst.afw.table.ExposureRecord`
            Record for a detector in a visit, with a WCS that has a FITS
            approximation.
        detector : `lsst.afw.cameraGeom.Detector`
            Camera geometry for this detector.

        Returns
        -------
        delta1 : `lsst.geom.Angle`
            Maximum error in the pixelToSky approximation.
        delta2 : `lsst.geom.Angle`
            Maximum error in the skyToPixel approximation.
        """
        target_wcs = record.getWcs()
        approx_wcs = target_wcs.getFitsApproximation()
        x, y = self._make_grid(detector, self.config.sip_grid_spacing, offset=True)
        target_ra, target_dec = target_wcs.pixelToSkyArray(x, y)
        approx1_ra, approx1_dec = approx_wcs.pixelToSkyArray(x, y)
        delta1 = float(np.max(np.hypot(target_ra - approx1_ra, target_dec - approx1_dec))) * radians
        approx2_x, approx2_y = approx_wcs.skyToPixelArray(target_ra, target_dec)
        approx2_ra, approx2_dec = target_wcs.pixelToSkyArray(approx2_x, approx2_y)
        delta2 = float(np.max(np.hypot(target_ra - approx2_ra, target_dec - approx2_dec))) * radians
        if self._delta1_key is not None:
            record.set(self._delta1_key, delta1)
        if self._delta2_key is not None:
            record.set(self._delta2_key, delta2)
        self.log.verbose(
            'SIP approximation for detector %s has pixelToSky good to ≤ %0.2g mas and '
            'skyToPixel good to ≤ %0.2g mas.',
            record.getId(),
            delta1.asArcseconds() * 1000,
            delta2.asArcseconds() * 1000,
        )
        return delta1, delta2

    def _make_visit_geometry(self, boresight, orientation, catalog, camera):
        """Create new sky regions for the visit and its detectors.

        Parameters
        ----------
        boresight : `lsst.geom.SpherePoint`
            New boresight location.
        orientation : `lsst.geom.Angle`
            New orientation angle.
        catalog : `lsst.afw.table.ExposureCatalog`
            A catalog of per-detector records for the visit with WCSs
            A repointed raw-like WCS will be used for any detectors not in the
            catalog or for which the catalog record does not have a WCS.
        camera : `lsst.afw.cameraGeom.Camera`
            Camera geometry.

        Returns
        -------
        regions : `lsst.obs.base.visit_geometry.VisitGeometry`
            Updated regions for the visit and all detectors in the camera.
        """
        detector_regions: dict[int, ConvexPolygon] = {}
        all_vertices = []
        for detector in camera:
            wcs = None
            if (record := catalog.find(detector.getId())) is not None:
                wcs = record.getWcs()
            if wcs is None:
                wcs = createInitialSkyWcsFromBoresight(boresight, orientation, detector)
            corners = wcs.pixelToSky(detector.getCorners(PIXELS))
            vertices = [sp.getVector() for sp in corners]
            detector_regions[detector.getId()] = ConvexPolygon(vertices)
            all_vertices.extend(vertices)
        visit_region = ConvexPolygon.convexHull(all_vertices)
        return VisitGeometry(
            boresight_ra=boresight.getRa().asDegrees(),
            boresight_dec=boresight.getDec().asDegrees(),
            orientation=orientation.asDegrees(),
            visit_region=visit_region,
            detector_regions=detector_regions,
        )

    def _make_grid(self, detector, spacing, offset=False) -> tuple[np.ndarray, np.ndarray]:
        pixel_bbox = Box2D(detector.getBBox())
        n_x = math.ceil(pixel_bbox.width / spacing)
        n_y = math.ceil(pixel_bbox.height / spacing)
        # We add one to the dimensions since there's a point at the min and max
        # in each dimension.
        xs = np.linspace(pixel_bbox.x.min, pixel_bbox.x.max, n_x + 1)
        ys = np.linspace(pixel_bbox.y.min, pixel_bbox.y.max, n_y + 1)
        if offset:
            xs = 0.5*(xs[1:] + xs[:-1])
            ys = 0.5*(ys[1:] + ys[:-1])
        x, y = np.meshgrid(xs, ys)
        return x.ravel(), y.ravel()
