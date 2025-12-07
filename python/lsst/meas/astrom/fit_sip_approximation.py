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

__all__ = ("FitSipApproximationConfig", "FitSipApproximationTask",)

import math


from lsst.afw.geom import SipApproximation
from lsst.geom import Box2D, Extent2I
from lsst.pex.config import Config, Field
from lsst.pipe.base import Struct, Task


class FitSipApproximationConfig(Config):
    grid_spacing = Field[float](
        doc=(
            "Spacing (in pixels) between grid points used to evaluate the WCS when fitting the "
            "approximation. "
            "If the spacing does not divide the detector bounding box evenly, it is decreased slightly."
        ),
        dtype=float,
        default=32.0,
    )
    order = Field[int](
        doc="Polynomial order for the SIP approximation.",
        dtype=int,
        default=5,
    )


class FitSipApproximationTask(Task):
    """A simple convenience wrapper for `lsst.afw.geom.SipApproximation`.
    """

    _DefaultName = "fitSipApproximation"
    ConfigClass = FitSipApproximationConfig

    def run(self, *, wcs, bbox):
        """Re-fit the pointing from the WCSs in a visit.

        Parameters
        ----------
        wcs : `lsst.afw.geom.SkyWcs`
            Target WCS to approximate.
        bbox : `lsst.geom.Box2I`
            The region where the WCS and its approximation are expected to be
            valid.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            A struct with the following attributes:

            - ``wcs`` (`lsst.afw.geom.SkyWcs`): a copy of the input ``wcs``
              with a SIP approximation attached.
            - ``delta_sky`` (`lsst.geom.Angle`): maximum separation in
              ``pixelToSky`` values on a grid offset from the one used for the
              fit.
            - ``delta_pixel`` (`float`): maximum separation in ``skyToPixel``
              values on a grid offset from the one used for the fit.
        """
        grid_shape = Extent2I(
            math.ceil(bbox.width / self.config.grid_spacing) + 1,
            math.ceil(bbox.height / self.config.grid_spacing) + 1,
        )
        approx = SipApproximation(wcs, Box2D(bbox), grid_shape, self.config.order)
        delta_sky, delta_pixel = approx.computeDeltas()
        self.log.verbose('Fit TAN-SIP approximation good to %0.2g" and %0.2g pixels.',
                         delta_sky.asArcseconds(), delta_pixel)
        return Struct(
            wcs=wcs.copyWithFitsApproximation(approx.getWcs()),
            delta_sky=delta_sky,
            delta_pixel=delta_pixel,
        )
