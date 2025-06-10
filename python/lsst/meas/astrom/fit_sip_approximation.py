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

__all__ = ()

import math
from typing import ClassVar

import numpy as np

from lsst.afw.geom import SkyWcs, getIntermediateWorldCoordsToSky, SipApproximation, makeTanSipWcs
from lsst.afw.cameraGeom import Detector, FIELD_ANGLE, PIXELS
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
from lsst.pipe.base import Struct, Task
from lsst.obs.base import Instrument


class FitSipApproximationConfig(Config):
    pointing_grid_spacing = Field[float](
        doc=(
            "Spacing (in pixels) between grid points used to evaluate the WCS when fitting the pointing "
            "This can be a very sparse grid (there are only three degrees of freedom in this fit).  "
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


class FitSipApproximationTask(Task):
    _DefaultName: ClassVar[str] = "fitSipApproximation"
    ConfigClass: ClassVar[type[FitSipApproximationConfig]] = FitSipApproximationConfig
    config: FitSipApproximationConfig

    def run(self, *, true_wcs: SkyWcs, detector: Detector, instrument: Instrument) -> Struct:
        boresight, orientation = self.fit_pointing(
            true_wcs=true_wcs, detector=detector, instrument=instrument
        )
        repointed_raw_wcs, _ = self.make_raw_wcs(
            boresight, orientation, detector=detector, instrument=instrument
        )
        sip_approx = self.fit_sip_approximation(
            true_wcs=true_wcs, repointed_raw_wcs=repointed_raw_wcs, detector=detector
        )
        fits_wcs = makeTanSipWcs(
            sip_approx.getPixelOrigin(),
            repointed_raw_wcs.getSkyOrigin(),
            sip_approx.getCdMatrix(),
            sip_approx.getA(),
            sip_approx.getB(),
            sip_approx.getAP(),
            sip_approx.getBP(),
        )
        return Struct(
            wcs=true_wcs.withFitsApproximation(fits_wcs),
            boresight=boresight,
            orientation=orientation,
            repointed_raw_wcs=repointed_raw_wcs,
        )

    def fit_pointing(
        self, *, true_wcs: SkyWcs, detector: Detector, instrument: Instrument
    ) -> tuple[SpherePoint, Angle]:
        start = true_wcs.pixelToSky(Point2D(0.0, 0.0))
        unpointed_raw_wcs, flip_x = self.make_raw_wcs(
            start, 0.0 * degrees, detector=detector, instrument=instrument
        )
        pixel_x, pixel_y = self._make_pointing_grid(detector)
        start_ra, start_dec = unpointed_raw_wcs.pixelToSkyArray(pixel_x, pixel_y)
        start_uv = np.stack(
            SpherePoint.toUnitXYZ(longitude=start_ra, latitude=start_dec, units=radians),
            axis=1,
        )
        true_ra, true_dec = true_wcs.pixelToSkyArray(pixel_x, pixel_y)
        true_uv = np.stack(
            SpherePoint.toUnitXYZ(longitude=true_ra, latitude=true_dec, units=radians),
            axis=1,
        )
        unpointed_y_axis_point = unpointed_raw_wcs.pixelToSky(
            detector.transform(Point2D(0.0, np.pi / 180.0), FIELD_ANGLE, PIXELS)
        )
        transform = SphereTransform.fit_unit_vectors(start_uv, true_uv)
        boresight = transform(start)
        transformed_y_axis_point = transform(unpointed_y_axis_point)
        orientation = Angle(90, degrees) - boresight.bearingTo(transformed_y_axis_point)
        if flip_x:
            raise NotImplementedError("flip_x=True is not implemented yet")
        return boresight, orientation

    def fit_sip_approximation(
        self, *, true_wcs: SkyWcs, repointed_raw_wcs: SkyWcs, detector: Detector
    ) -> SkyWcs:
        pixels_to_iwc = true_wcs.getTransform().then(
            getIntermediateWorldCoordsToSky(repointed_raw_wcs, simplify=True).inverted()
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
            repointed_raw_wcs.getPixelOrigin(),
            repointed_raw_wcs.getCdMatrix(),
            pixel_bbox,
            grid_shape,
            self.config.sip_order,
        )
        return sip_approx

    def make_raw_wcs(
        self, boresight: SpherePoint, orientation: Angle, *, detector: Detector, instrument: Instrument
    ) -> tuple[SkyWcs, bool]:
        raw_formatter = instrument.getRawFormatter({"detector": detector.getId()})
        return (
            raw_formatter.makeRawSkyWcsFromBoresight(boresight, orientation, detector),
            raw_formatter.wcsFlipX,
        )

    def _make_pointing_grid(self, detector: Detector) -> tuple[np.ndarray, np.ndarray]:
        pixel_bbox = Box2D(detector.getBBox())
        n_x = math.ceil(pixel_bbox.width / self.config.pointing_grid_spacing)
        n_y = math.ceil(pixel_bbox.height / self.config.pointing_grid_spacing)
        x, y = np.meshgrid(
            # We add one to the dimensions since there's a point the min and
            # max in each dimension.
            np.linspace(pixel_bbox.x.min, pixel_bbox.x.max, n_x + 1),
            np.linspace(pixel_bbox.y.min, pixel_bbox.y.max, n_y + 1),
        )
        return x.ravel(), y.ravel()
