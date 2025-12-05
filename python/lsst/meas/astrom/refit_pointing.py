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

from lsst.geom import (
    Angle,
    Box2D,
    Point2D,
    SpherePoint,
    SphereTransform,
    arcseconds,
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
    grid_spacing = Field[float](
        doc=(
            "Spacing (in pixels) between grid points used to evaluate the WCS when fitting the pointing. "
            "This can be a very sparse grid (there are only three degrees of freedom).  "
            "If the spacing does not divide the detector bounding box evenly, it is decreased slightly."
        ),
        dtype=float,
        default=512.0,
    )
    rejection_threshold = Field(
        doc=(
            "If the distance between the target WCS position and the position predicted by the camera "
            "geometry after refitting the pointing using just one detector exceeds this value (in "
            "arcseconds) at any point on the pointing-fit grid, that detector is rejected from the pointing "
            "fit.  The quantity this threshold is applied to is saved in the wcs_detector_pointing_residual "
            "column."
        ),
        dtype=float,
        default=10.0,
    )
    nulling_threshold = Field(
        doc=(
            "If the distance between the target WCS position and the position predicted by the camera "
            "geometry after refitting the pointing using all detectors exceeds this value (in arcseconds) "
            "at any point on the pointing-fit grid, that detector's WCS is set to None in the catalog.  "
            "The quantity this threshold is applied to is saved in the wcs_visit_pointing_residual column."
        ),
        dtype=float,
        default=60.0,
    )
    schema_prefix = Field(
        doc="Prefix for all schema fields.",
        dtype=str,
        default="",
    )
    fallback_region_padding = Field(
        doc=(
            "Padding to add (in pixels) to the regions of detectors for which only a "
            "pointing + camera geometry WCS is available."
        ),
        dtype=int,
        default=50,
    )


class RefitPointingTask(Task):
    """A task that uses the available WCSs of the detectors in a visit to
    re-fit the pointing for that visit and compute new visit regions for the butler.
    """

    _DefaultName = "refitPointing"
    ConfigClass = RefitPointingConfig

    def __init__(self, config=None, *, schema, **kwargs):
        super().__init__(config, **kwargs)
        self._detector_pointing_residual_key = schema.addField(
            self.config.schema_prefix + "wcs_detector_pointing_residual",
            type="Angle",
            doc=(
                "Maximum difference (on the pointing-fit grid) between the target WCS position and "
                "the position predicted by camera geometry, after re-pointing using the target WCS "
                "for this detector only."
            ),
        )
        self._visit_pointing_residual_key = schema.addField(
            self.config.schema_prefix + "wcs_visit_pointing_residual",
            type="Angle",
            doc=(
                "Maximum difference (on the pointing-fit grid) between the target WCS position and "
                "the position predicted by camera geometry, after re-pointing using the target WCS "
                "of all non-rejected detectors in the visit."
            ),
        )
        self._rejected_key = schema.addField(
            self.config.schema_prefix + "wcs_detector_pointing_rejected",
            type="Flag",
            doc=(
                "Flag set if this detector was rejected from the pointing fit due to its "
                "wcs_detector_pointing_residual value."
            ),
        )
        self._rejection_threshold = self.config.rejection_threshold * arcseconds
        self._nulling_threshold = self.config.nulling_threshold * arcseconds

    def run(self, *, catalog, camera):
        """Re-fit the pointing from the WCSs in a visit.

        Parameters
        ----------
        catalog : `lsst.afw.table.ExposureCatalog`
            A catalog of per-detector records for the visit.  Columns with WCS
            diagnostics are updatd in-place, and WCSs may be set to `None` if
            they do not satisfy the `~RefitPointingConfig.nulling_threshold`.
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
        self._null_bad(catalog)
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
        start_boresight: SpherePoint | None = None
        start_orientation = 0.0 * degrees
        start_y_axis_point: SpherePoint | None = None
        detectors_kept: list[int] = []
        start_xyz: dict[int, np.ndarray] = {}
        target_xyz: dict[int, np.ndarray] = {}
        for record in catalog:
            detector_id = record.getId()
            # We call the WCSs that were actually fit to the stars the "true"
            # WCSs.
            target_wcs = record.getWcs()
            if target_wcs is None:
                continue
            try:
                detector = camera[detector_id]
            except LookupError:
                self.log.warning("Detector %d has no camera geometry; skipping it.", detector_id)
                continue
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
            pixel_x, pixel_y = self._make_grid(detector, self.config.grid_spacing)
            start_ra, start_dec = start_wcs.pixelToSkyArray(pixel_x, pixel_y)
            start_xyz[detector_id] = np.stack(
                SpherePoint.toUnitXYZ(longitude=start_ra, latitude=start_dec, units=radians),
                axis=1,
            )
            target_ra, target_dec = target_wcs.pixelToSkyArray(pixel_x, pixel_y)
            target_xyz[detector_id] = np.stack(
                SpherePoint.toUnitXYZ(longitude=target_ra, latitude=target_dec, units=radians),
                axis=1,
            )
            # Fit the pointing using just the grid for this detector to see if
            # the residuals are any good; they won't be if the target WCS is
            # bonkers and makes the detector non-rectangular on the sky.
            detector_transform = SphereTransform.fit_unit_vectors(
                start_xyz[detector_id],
                target_xyz[detector_id],
            )
            detector_pointing_residual = self._compute_pointing_residual(
                detector_transform, start_xyz[detector_id], target_xyz[detector_id]
            )
            record.set(self._detector_pointing_residual_key, detector_pointing_residual)
            if detector_pointing_residual > self._rejection_threshold:
                record.set(self._rejected_key, True)
                if not detectors_kept:
                    # This was the first detector we saw; need to reset.
                    start_boresight = None
                self.log.warning(
                    'Dropping detector %d with detector pointing residual %0.2g" from pointing fit.',
                    detector_id,
                    detector_pointing_residual.asArcseconds(),
                )
                continue
            detectors_kept.append(detector_id)
        if not detectors_kept:
            # Since we can't apply the nulling-threshold test, set all WCSs to
            # None.
            for record in catalog:
                record.setWcs(None)
            raise NoVisitWcs("No valid target WCSs were left after rejection.")
        # Fit the spherical rotation that maps the points in the arbitrary
        # start WCS to the target WCS, using all kept detectors.
        transform = SphereTransform.fit_unit_vectors(
            np.concatenate([start_xyz[i] for i in detectors_kept]),
            np.concatenate([target_xyz[i] for i in detectors_kept]),
        )
        # Compute and record the residuals for each detector with this
        # transform.
        for record in catalog:
            detector_id = record.getId()
            if detector_id not in start_xyz:
                # This detector already doesn't have a WCS.
                continue
            visit_pointing_residual = self._compute_pointing_residual(
                transform, start_xyz[detector_id], target_xyz[detector_id]
            )
            record.set(self._visit_pointing_residual_key, visit_pointing_residual)
        # If we apply that same rotation to our arbitrary start boresight, we
        # get the boresight predicted by the target WCSs.
        boresight = transform(start_boresight)
        # If we apply that rotation to a point on the FIELD_ANGLE y-axis, we
        # can similarly recover the orientation angle predicted by the target
        # WCSs.
        start_y_axis_point = start_boresight.offset(90 * degrees, 1.0 * degrees)
        transformed_y_axis_point = transform(start_y_axis_point)
        orientation = Angle(90, degrees) - boresight.bearingTo(transformed_y_axis_point)
        if camera.getFocalPlaneParity():
            raise NotImplementedError("Cameras with focal plane parity flips are not yet supported.")
        return boresight, orientation

    def _compute_pointing_residual(self, transform, from_xyz, to_xyz):
        # Apply the transform to the start positions and subtract the
        # target positions (all in 3-vector space) to get the residual
        # 3-vectors.
        residual_vecs = np.dot(transform.matrix, from_xyz.transpose()).transpose()
        residual_vecs -= to_xyz
        # Compute the squared chord length of the residual vectors, find
        # the maximum of that over the grid (since everything else we do
        # is monotonic), then translate that into an angle.
        return 2.0 * np.arcsin(0.5 * np.sum(residual_vecs**2, axis=1).max() ** 0.5) * radians

    def _null_bad(self, catalog):
        for record in catalog:
            visit_pointing_residual = record.get(self._visit_pointing_residual_key)
            if visit_pointing_residual > self._nulling_threshold:
                self.log.warning(
                    'Setting WCS to None for detector %d with visit pointing residual %0.2g".',
                    record.getId(),
                    visit_pointing_residual.asArcseconds(),
                )
                record.setWcs(None)

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
            pixel_bbox = Box2D(detector.getBBox())
            if wcs is None or record.get(self._rejected_key):
                wcs = createInitialSkyWcsFromBoresight(boresight, orientation, detector)
                pixel_bbox.grow(self.config.fallback_region_padding)
            corners = wcs.pixelToSky(pixel_bbox.getCorners())
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

    def _make_grid(self, detector, spacing) -> tuple[np.ndarray, np.ndarray]:
        pixel_bbox = Box2D(detector.getBBox())
        n_x = math.ceil(pixel_bbox.width / spacing)
        n_y = math.ceil(pixel_bbox.height / spacing)
        # We add one to the dimensions since there's a point at the min and max
        # in each dimension.
        xs = np.linspace(pixel_bbox.x.min, pixel_bbox.x.max, n_x + 1)
        ys = np.linspace(pixel_bbox.y.min, pixel_bbox.y.max, n_y + 1)
        x, y = np.meshgrid(xs, ys)
        return x.ravel(), y.ravel()
