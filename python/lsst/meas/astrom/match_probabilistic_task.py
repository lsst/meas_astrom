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

from typing import Dict, List, Optional, Set, Tuple
import warnings

import astropy.table
import logging
import lsst.afw.geom as afwGeom
import lsst.geom as geom
import lsst.pipe.base as pipeBase
import lsst.utils as utils
import numpy as np
import pandas as pd

from .matcher_probabilistic import MatchProbabilisticConfig, MatcherProbabilistic


__all__ = ["MatchProbabilisticTask", "radec_to_xy"]


def radec_to_xy(ra_vec, dec_vec, factor, wcs: afwGeom.SkyWcs):
    radec_true = [
        geom.SpherePoint(ra*factor, dec*factor, geom.degrees)
        for ra, dec in zip(ra_vec, dec_vec)
    ]
    return wcs.skyToPixel(radec_true)


class MatchProbabilisticTask(pipeBase.Task):
    """Run MatchProbabilistic on a reference and target catalog covering the same tract."""

    ConfigClass = MatchProbabilisticConfig
    _DefaultName = "matchProbabilistic"

    @staticmethod
    def _apply_select_bool(
        catalog: astropy.table.Table | pd.DataFrame,
        columns_true: List[str],
        columns_false: List[str],
        selection: Optional[np.array],
    ) -> np.array:
        """Apply additional boolean selection columns.

        catalog : `pandas.DataFrame` | `astropy.table.Table`
            The catalog to select from.
        columns_true : `list` [`str`]
            Columns that must be True for selection.
        columns_false : `list` [`str`]
            Columns that must be False for selection.
        selection : `numpy.array`
            A prior selection array. Default all true.

        Returns
        -------
        selection : `numpy.array`
            The final selection array.

        """
        # TODO: Remove pandas support in DM-46523
        is_pd = isinstance(catalog, pd.DataFrame)
        if is_pd:
            warnings.warn("pandas usage in MatchProbabilisticTask is deprecated; it will be removed "
                          " in favour of astropy.table after release 28.0.0", category=FutureWarning)
        select_additional = (len(columns_true) + len(columns_false)) > 0
        if select_additional:
            if selection is None:
                selection = np.ones(len(catalog), dtype=bool)
            for column in columns_true:
                # This is intended for boolean columns, so the behaviour for non-boolean is not obvious
                # More config options and/or using a ConfigurableActionField might be best
                values = catalog[column] if not is_pd else catalog[column].values
                selection &= (np.isfinite(values) & (values != 0))
            for column in columns_false:
                values = catalog[column] if not is_pd else catalog[column].values
                selection &= (values == 0)
        return selection

    @property
    def columns_in_ref(self) -> Set[str]:
        return self.config.columns_in_ref

    @property
    def columns_in_target(self) -> Set[str]:
        return self.config.columns_in_target

    def match(
        self,
        catalog_ref: astropy.table.Table | pd.DataFrame,
        catalog_target: astropy.table.Table | pd.DataFrame,
        select_ref: np.array = None,
        select_target: np.array = None,
        wcs: afwGeom.SkyWcs = None,
        logger: logging.Logger = None,
        logging_n_rows: int = None,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[int, str]]:
        """Match sources in a reference tract catalog with a target catalog.

        Parameters
        ----------
        catalog_ref : `pandas.DataFrame` | `astropy.table.Table`
            A reference catalog to match objects/sources from.
        catalog_target : `pandas.DataFrame` | `astropy.table.Table`
            A target catalog to match reference objects/sources to.
        select_ref : `numpy.array`
            A boolean array of the same length as `catalog_ref` selecting the sources that can be matched.
        select_target : `numpy.array`
            A boolean array of the same length as `catalog_target` selecting the sources that can be matched.
        wcs : `lsst.afw.image.SkyWcs`
            A coordinate system to convert catalog positions to sky coordinates. Only used if
            `self.config.coords_ref_to_convert` is set.
        logger : `logging.Logger`
            A Logger for logging.
        logging_n_rows : `int`
            Number of matches to make before outputting incremental log message.

        Returns
        -------
        catalog_out_ref : `pandas.DataFrame` | `astropy.table.Table`
            Reference matched catalog with indices of target matches.
        catalog_out_target : `pandas.DataFrame` | `astropy.table.Table`
            Reference matched catalog with indices of target matches.
        """
        if logger is None:
            logger = self.log

        # TODO: Remove pandas support in DM-46523
        is_ref_pd = isinstance(catalog_ref, pd.DataFrame)
        is_target_pd = isinstance(catalog_target, pd.DataFrame)
        if is_ref_pd or is_target_pd:
            warnings.warn("pandas usage in MatchProbabilisticTask is deprecated; it will be removed "
                          " in favour of astropy.table after release 28.0.0", category=FutureWarning)

        config = self.config

        if config.column_ref_order is None:
            fluxes = (
                catalog_ref.loc[:, config.columns_ref_flux].values
                if is_ref_pd else
                [catalog_ref[key] for key in config.columns_ref_flux]
            )
            flux_tot = np.nansum(fluxes, axis=1 if is_ref_pd else 0)
            catalog_ref["flux_total"] = flux_tot
            if config.mag_brightest_ref != -np.inf or config.mag_faintest_ref != np.inf:
                mag_tot = (
                    -2.5 * np.log10(flux_tot) + config.coord_format.mag_zeropoint_ref
                )
                select_mag = (mag_tot >= config.mag_brightest_ref) & (
                    mag_tot <= config.mag_faintest_ref
                )
            else:
                select_mag = np.isfinite(flux_tot)
            if select_ref is None:
                select_ref = select_mag
            else:
                select_ref &= select_mag

        with warnings.catch_warnings():
            # We already issued a deprecation warning; no need to repeat it.
            warnings.filterwarnings(action="ignore", category=FutureWarning)
            select_ref = self._apply_select_bool(
                catalog=catalog_ref,
                columns_true=config.columns_ref_select_true,
                columns_false=config.columns_ref_select_false,
                selection=select_ref,
            )
            select_target = self._apply_select_bool(
                catalog=catalog_target,
                columns_true=config.columns_target_select_true,
                columns_false=config.columns_target_select_false,
                selection=select_target,
            )

        logger.info(
            "Beginning MatcherProbabilistic.match with %d/%d ref sources selected vs %d/%d target",
            len(catalog_ref) if select_ref is None else np.sum(select_ref),
            len(catalog_ref),
            len(catalog_target) if select_target is None else np.sum(select_target),
            len(catalog_target),
        )

        catalog_out_ref, catalog_out_target, exceptions = self.matcher.match(
            catalog_ref,
            catalog_target,
            select_ref=select_ref,
            select_target=select_target,
            logger=logger,
            logging_n_rows=logging_n_rows,
            wcs=wcs,
            radec_to_xy_func=radec_to_xy,
        )

        return catalog_out_ref, catalog_out_target, exceptions

    @utils.timer.timeMethod
    def run(
        self,
        catalog_ref: astropy.table.Table | pd.DataFrame,
        catalog_target: astropy.table.Table | pd.DataFrame,
        wcs: afwGeom.SkyWcs = None,
        **kwargs,
    ) -> pipeBase.Struct:
        """Match sources in a reference tract catalog with a target catalog.

        Parameters
        ----------
        catalog_ref : `pandas.DataFrame` | `astropy.table.Table`
            A reference catalog to match objects/sources from.
        catalog_target : `pandas.DataFrame` | `astropy.table.Table`
            A target catalog to match reference objects/sources to.
        wcs : `lsst.afw.image.SkyWcs`
            A coordinate system to convert catalog positions to sky coordinates.
            Only needed if `config.coords_ref_to_convert` is used to convert
            reference catalog sky coordinates to pixel positions.
        kwargs : Additional keyword arguments to pass to `match`.

        Returns
        -------
        retStruct : `lsst.pipe.base.Struct`
            A struct with output_ref and output_target attribute containing the
            output matched catalogs, as well as a dict
        """
        # TODO: Remove pandas support in DM-46523
        is_ref_pd = isinstance(catalog_ref, pd.DataFrame)
        is_target_pd = isinstance(catalog_target, pd.DataFrame)
        if is_ref_pd:
            catalog_ref.reset_index(inplace=True)
        if is_target_pd:
            catalog_target.reset_index(inplace=True)
        if is_ref_pd or is_target_pd:
            warnings.warn("pandas usage in MatchProbabilisticTask is deprecated; it will be removed "
                          " in favour of astropy.table after release 28.0.0", category=FutureWarning)
        with warnings.catch_warnings():
            # We already issued a deprecation warning; no need to repeat it.
            warnings.filterwarnings(action="ignore", category=FutureWarning)
            catalog_ref, catalog_target, exceptions = self.match(
                catalog_ref, catalog_target, wcs=wcs, **kwargs
            )

        return pipeBase.Struct(
            cat_output_ref=catalog_ref,
            cat_output_target=catalog_target,
            exceptions=exceptions,
        )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.matcher = MatcherProbabilistic(self.config)
