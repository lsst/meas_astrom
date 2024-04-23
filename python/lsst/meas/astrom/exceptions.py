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

__all__ = ["AstrometryError", "AstrometryFitFailure", "BadAstrometryFit", "MatcherFailure"]

import lsst.pipe.base


class AstrometryError(lsst.pipe.base.AlgorithmError):
    """Parent class for failures in astrometric fitting.

    Parameters
    ----------
    msg : `str`
        Informative message about the nature of the error.
    **kwargs
        All other arguments are added to a ``_metadata`` attribute, which is
        used to generate the metadata property for Task annotation.
    """
    def __init__(self, msg, **kwargs):
        self.msg = msg
        self._metadata = kwargs
        super().__init__(msg, kwargs)

    def __str__(self):
        # Exception doesn't handle **kwargs, so we need a custom str.
        return f"{self.msg}: {self.metadata}"

    @property
    def metadata(self):
        for key, value in self._metadata.items():
            if not (isinstance(value, int) or isinstance(value, float) or isinstance(value, str)):
                raise TypeError(f"{key} is of type {type(value)}, but only (int, float, str) are allowed.")
        return self._metadata


class BadAstrometryFit(AstrometryError):
    """Raised if the quality of the astrometric fit is worse than some
    threshold.

    Parameters
    ----------
    distMean : `float`
        Mean on-sky separation of matched sources, in arcseconds.
    distMedian : `float`
        Median on-sky separation of matched sources, in arcseconds.
    """
    def __init__(self, distMean, maxMeanDist, distMedian, **kwargs):
        msg = f'Poor quality astrometric fit, {distMean}" > {maxMeanDist}"'
        super().__init__(msg, **kwargs)
        self._metadata["distMean"] = distMean
        self._metadata["maxMeanDist"] = distMean
        self._metadata["distMedian"] = distMedian


class AstrometryFitFailure(AstrometryError):
    """Raised if the astrometry fitter fails."""


class MatcherFailure(AstrometryError):
    """Raised if the matcher fails."""
