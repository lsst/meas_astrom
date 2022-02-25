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

__all__ = ['ConvertCatalogCoordinatesConfig', 'MatchProbabilisticConfig', 'MatcherProbabilistic']

import lsst.pex.config as pexConfig

from dataclasses import dataclass
import logging
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import time
from typing import Callable, Set

logger_default = logging.getLogger(__name__)


def _mul_column(column: np.array, value: float):
    if value is not None and value != 1:
        column *= value
    return column


def _radec_to_xyz(ra, dec):
    """Convert input ra/dec coordinates to spherical unit vectors.

    Parameters
    ----------
    ra, dec: `numpy.ndarray`
        Arrays of right ascension/declination in degrees.

    Returns
    -------
    vectors : `numpy.ndarray`, (N, 3)
        Output unit vectors.
    """
    if ra.size != dec.size:
        raise ValueError('ra and dec must be same size')
    ras = np.radians(ra)
    decs = np.radians(dec)
    vectors = np.empty((ras.size, 3))

    sin_dec = np.sin(np.pi / 2 - decs)
    vectors[:, 0] = sin_dec * np.cos(ras)
    vectors[:, 1] = sin_dec * np.sin(ras)
    vectors[:, 2] = np.cos(np.pi / 2 - decs)

    return vectors


@dataclass
class CatalogExtras:
    """Store frequently-reference (meta)data relevant for matching a catalog.

    Parameters
    ----------
    catalog : `pandas.DataFrame`
        A pandas catalog to store extra information for.
    select : `numpy.array`
        A numpy boolean array of the same length as catalog to be used for
        target selection.
    """
    n: int
    indices: np.array
    select: np.array

    coordinate_factor: float = None

    def __init__(self, catalog: pd.DataFrame, select: np.array = None, coordinate_factor: float = None):
        self.n = len(catalog)
        self.select = np.ones(self.n, dtype=bool) if select is None else select
        self.indices = np.flatnonzero(select) if select is not None else np.arange(self.n)
        self.coordinate_factor = coordinate_factor


@dataclass(frozen=True)
class ComparableCatalog:
    """A catalog with sources with coordinate columns in some standard format/units.

    catalog : `pandas.DataFrame`
        A catalog with comparable coordinate columns.
    column_coord1 : `str`
        The first spatial coordinate column name.
    column_coord2 : `str`
        The second spatial coordinate column name.
    coord1 : `numpy.array`
        The first spatial coordinate values.
    coord2 : `numpy.array`
        The second spatial coordinate values.
    extras : `CatalogExtras`
        Extra cached (meta)data for the `catalog`.
    """
    catalog: pd.DataFrame
    column_coord1: str
    column_coord2: str
    coord1: np.array
    coord2: np.array
    extras: CatalogExtras


class ConvertCatalogCoordinatesConfig(pexConfig.Config):
    """Configuration for the MatchProbabilistic matcher.
    """
    column_ref_coord1 = pexConfig.Field(
        dtype=str,
        default='ra',
        doc='The reference table column for the first spatial coordinate (usually x or ra).',
    )
    column_ref_coord2 = pexConfig.Field(
        dtype=str,
        default='dec',
        doc='The reference table column for the second spatial coordinate (usually y or dec).'
            'Units must match column_ref_coord1.',
    )
    column_target_coord1 = pexConfig.Field(
        dtype=str,
        default='coord_ra',
        doc='The target table column for the first spatial coordinate (usually x or ra).'
            'Units must match column_ref_coord1.',
    )
    column_target_coord2 = pexConfig.Field(
        dtype=str,
        default='coord_dec',
        doc='The target table column for the second spatial coordinate (usually y or dec).'
            'Units must match column_ref_coord2.',
    )
    coords_spherical = pexConfig.Field(
        dtype=bool,
        default=True,
        doc='Whether column_*_coord[12] are spherical coordinates (ra/dec) or not (pixel x/y)',
    )
    coords_ref_factor = pexConfig.Field(
        dtype=float,
        default=1.0,
        doc='Multiplicative factor for reference catalog coordinates.'
            'If coords_spherical is true, this must be the number of degrees per unit increment of '
            'column_ref_coord[12]. Otherwise, it must convert the coordinate to the same units'
            ' as the target coordinates.',
    )
    coords_target_factor = pexConfig.Field(
        dtype=float,
        default=1.0,
        doc='Multiplicative factor for target catalog coordinates.'
            'If coords_spherical is true, this must be the number of degrees per unit increment of '
            'column_target_coord[12]. Otherwise, it must convert the coordinate to the same units'
            ' as the reference coordinates.',
    )
    coords_ref_to_convert = pexConfig.DictField(
        default=None,
        keytype=str,
        itemtype=str,
        dictCheck=lambda x: len(x) == 2,
        doc='Dict mapping sky coordinate columns to be converted to pixel columns',
    )
    mag_zeropoint_ref = pexConfig.Field(
        dtype=float,
        default=31.4,
        doc='Magnitude zeropoint for reference catalog.',
    )

    def format_catalogs(
        self,
        catalog_ref: pd.DataFrame,
        catalog_target: pd.DataFrame,
        select_ref: np.array = None,
        select_target: np.array = None,
        radec_to_xy_func: Callable = None,
        return_converted_columns: bool = False,
        **kwargs,
    ):
        """Format matched catalogs that may require coordinate conversions.

        Parameters
        ----------
        catalog_ref : `pandas.DataFrame`
            A reference catalog for comparison to `catalog_target`.
        catalog_target : `pandas.DataFrame`
            A target catalog with measurements for comparison to `catalog_ref`.
        select_ref : `numpy.ndarray`, (Nref,)
            A boolean array of len `catalog_ref`, True for valid match candidates.
        select_target : `numpy.ndarray`, (Ntarget,)
            A boolean array of len `catalog_target`, True for valid match candidates.
        radec_to_xy_func : `typing.Callable`
            Function taking equal-length ra, dec arrays and returning an ndarray of
            - ``x``: current parameter (`float`).
            - ``extra_args``: additional arguments (`dict`).
        return_converted_columns : `bool`
            Whether to return converted columns in the `coord1` and `coord2`
            attributes, rather than keep the original values.
        kwargs

        Returns
        -------

        """
        convert_ref = self.coords_ref_to_convert
        if convert_ref and not callable(radec_to_xy_func):
            raise TypeError('radec_to_xy_func must be callable if converting ref coords')

        # Set up objects with frequently-used attributes like selection bool array
        extras_ref, extras_target = (
            CatalogExtras(catalog, select=select, coordinate_factor=coord_factor)
            for catalog, select, coord_factor in zip(
                (catalog_ref, catalog_target),
                (select_ref, select_target),
                (self.coords_ref_factor, self.coords_target_factor),
            )
        )

        compcats = []

        # Retrieve coordinates and multiply them by scaling factors
        for catalog, extras, (column1, column2), convert in (
            (catalog_ref, extras_ref, (self.column_ref_coord1, self.column_ref_coord2), convert_ref),
            (catalog_target, extras_target, (self.column_target_coord1, self.column_target_coord2), False),
        ):
            coord1, coord2 = (
                _mul_column(catalog.loc[:, column].values, extras.coordinate_factor)
                for column in (column1, column2)
            )
            if convert:
                xy_ref = radec_to_xy_func(coord1, coord2, self.coords_ref_factor, **kwargs)
                for idx_coord, column_out in enumerate(self.coords_ref_to_convert.values()):
                    coord = np.array([xy[idx_coord] for xy in xy_ref])
                    catalog[column_out] = coord
            if convert_ref and return_converted_columns:
                column1, column2 = self.coords_ref_to_convert.values()
                coord1, coord2 = catalog[column1], catalog[column2]
            if isinstance(coord1, pd.Series):
                coord1 = coord1.values
            if isinstance(coord2, pd.Series):
                coord2 = coord2.values

            compcats.append(ComparableCatalog(
                catalog=catalog, column_coord1=column1, column_coord2=column2,
                coord1=coord1, coord2=coord2, extras=extras,
            ))

        return tuple(compcats)


class MatchProbabilisticConfig(pexConfig.Config):
    """Configuration for the MatchProbabilistic matcher.
    """
    column_order = pexConfig.Field(
        dtype=str,
        default=None,
        optional=True,
        doc="Column to sort fit by. Derived from columns_ref_flux if not set.",
    )

    @property
    def columns_in_ref(self) -> Set[str]:
        columns_all = [x for x in self.columns_ref_flux]
        for columns in (
            [self.coord_format.column_ref_coord1, self.coord_format.column_ref_coord2],
            self.columns_ref_meas,
        ):
            columns_all.extend(columns)

        return set(columns_all)

    @property
    def columns_in_target(self) -> Set[str]:
        columns_all = [self.coord_format.column_target_coord1, self.coord_format.column_target_coord2]
        for columns in (
            self.columns_target_meas,
            self.columns_target_err,
            self.columns_target_select_false,
            self.columns_target_select_true,
        ):
            columns_all.extend(columns)
        return set(columns_all)

    columns_ref_flux = pexConfig.ListField(
        dtype=str,
        default=[],
        listCheck=lambda x: len(set(x)) == len(x),
        optional=True,
        doc="List of reference flux columns to nansum total magnitudes from if column_order is None",
    )
    columns_ref_meas = pexConfig.ListField(
        dtype=str,
        doc='The reference table columns to compute match likelihoods from '
            '(usually centroids and fluxes/magnitudes)',
    )
    columns_target_meas = pexConfig.ListField(
        dtype=str,
        doc='Target table columns with measurements corresponding to columns_ref_meas',
    )
    columns_target_err = pexConfig.ListField(
        dtype=str,
        doc='Target table columns with standard errors (sigma) corresponding to columns_ref_meas',
    )
    columns_target_select_true = pexConfig.ListField(
        dtype=str,
        default=('detect_isPrimary',),
        doc='Target table columns to require to be True for selecting sources',
    )
    columns_target_select_false = pexConfig.ListField(
        dtype=str,
        default=('merge_peak_sky',),
        doc='Target table columns to require to be False for selecting sources',
    )
    coord_format = pexConfig.ConfigField(
        dtype=ConvertCatalogCoordinatesConfig,
        doc="Configuration for coordinate conversion",
    )
    mag_brightest_ref = pexConfig.Field(
        dtype=float,
        default=-np.inf,
        doc='Bright magnitude cutoff for selecting reference sources to match.'
            ' Ignored if column_order is None.'
    )
    mag_faintest_ref = pexConfig.Field(
        dtype=float,
        default=np.Inf,
        doc='Faint magnitude cutoff for selecting reference sources to match.'
            ' Ignored if column_order is None.'
    )
    match_dist_max = pexConfig.Field(
        dtype=float,
        default=0.5,
        doc='Maximum match distance. Units must be arcseconds if coords_spherical, '
            'or else match those of column_*_coord[12] multiplied by coords_*_factor.',
    )
    match_n_max = pexConfig.Field(
        dtype=int,
        default=10,
        optional=True,
        doc='Maximum number of spatial matches to consider (in ascending distance order).',
    )
    match_n_finite_min = pexConfig.Field(
        dtype=int,
        default=3,
        optional=True,
        doc='Minimum number of columns with a finite value to measure match likelihood',
    )
    order_ascending = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc='Whether to order reference match candidates in ascending order of column_order '
            '(should be False if the column is a flux and True if it is a magnitude.',
    )


class MatcherProbabilistic:
    """A probabilistic, greedy catalog matcher.

    Parameters
    ----------
    config: `MatchProbabilisticConfig`
        A configuration instance.
    """
    config: MatchProbabilisticConfig

    def __init__(
            self,
            config: MatchProbabilisticConfig,
    ):
        self.config = config

    def match(
            self,
            catalog_ref: pd.DataFrame,
            catalog_target: pd.DataFrame,
            select_ref: np.array = None,
            select_target: np.array = None,
            logger: logging.Logger = None,
            logging_n_rows: int = None,
            **kwargs
    ):
        """Match catalogs.

        Parameters
        ----------
        catalog_ref : `pandas.DataFrame`
            A reference catalog to match in order of a given column (i.e. greedily).
        catalog_target : `pandas.DataFrame`
            A target catalog for matching sources from `catalog_ref`. Must contain measurements with errors.
        select_ref : `numpy.array`
            A boolean array of the same length as `catalog_ref` selecting the sources that can be matched.
        select_target : `numpy.array`
            A boolean array of the same length as `catalog_target` selecting the sources that can be matched.
        logger : `logging.Logger`
            A Logger for logging.
        logging_n_rows : `int`
            The number of sources to match before printing a log message.
        kwargs
            Additional keyword arguments to pass to `format_catalogs`.

        Returns
        -------
        catalog_out_ref : `pandas.DataFrame`
            A catalog of identical length to `catalog_ref`, containing match information for rows selected by
            `select_ref` (including the matching row index in `catalog_target`).
        catalog_out_target : `pandas.DataFrame`
            A catalog of identical length to `catalog_target`, containing the indices of matching rows in
            `catalog_ref`.
        exceptions : `dict` [`int`, `Exception`]
            A dictionary keyed by `catalog_target` row number of the first exception caught when matching.
        """
        if logger is None:
            logger = logger_default

        config = self.config

        # Transform any coordinates, if required
        # Note: The returned objects contain the original catalogs, as well as
        # transformed coordinates, and the selection of sources for matching.
        # These might be identical to the arrays passed as kwargs, but that
        # depends on config settings.
        # For the rest of this function, the selection arrays will be used,
        # but the indices of the original, unfiltered catalog will also be
        # output, so some further indexing steps are needed.
        ref, target = config.coord_format.format_catalogs(
            catalog_ref=catalog_ref, catalog_target=catalog_target,
            select_ref=select_ref, select_target=select_target, **kwargs
        )

        n_ref_select = len(ref.extras.indices)

        match_dist_max = config.match_dist_max
        coords_spherical = config.coord_format.coords_spherical
        if coords_spherical:
            match_dist_max = np.radians(match_dist_max / 3600.)

        # Convert ra/dec sky coordinates to spherical vectors for accurate distances
        func_convert = _radec_to_xyz if coords_spherical else np.vstack
        vec_ref, vec_target = (
            func_convert(cat.coord1[cat.extras.select], cat.coord2[cat.extras.select])
            for cat in (ref, target)
        )

        # Generate K-d tree to compute distances
        logger.info('Generating cKDTree with match_n_max=%d', config.match_n_max)
        tree_obj = cKDTree(vec_target)

        scores, idxs_target_select = tree_obj.query(
            vec_ref,
            distance_upper_bound=match_dist_max,
            k=config.match_n_max,
        )

        n_target_select = len(target.extras.indices)
        n_matches = np.sum(idxs_target_select != n_target_select, axis=1)
        n_matched_max = np.sum(n_matches == config.match_n_max)
        if n_matched_max > 0:
            logger.warning(
                '%d/%d (%.2f%%) selected true objects have n_matches=n_match_max(%d)',
                n_matched_max, n_ref_select, 100.*n_matched_max/n_ref_select, config.match_n_max
            )

        # Pre-allocate outputs
        target_row_match = np.full(target.extras.n, np.nan, dtype=np.int64)
        ref_candidate_match = np.zeros(ref.extras.n, dtype=bool)
        ref_row_match = np.full(ref.extras.n, np.nan, dtype=int)
        ref_match_count = np.zeros(ref.extras.n, dtype=int)
        ref_match_meas_finite = np.zeros(ref.extras.n, dtype=int)
        ref_chisq = np.full(ref.extras.n, np.nan, dtype=float)

        # If no order is specified, take nansum of all flux columns for a 'total flux'
        # Note: it won't actually be a total flux if bands overlap significantly
        # (or it might define a filter with >100% efficiency
        column_order = (
            catalog_ref.loc[ref.extras.select, config.column_order].values
            if config.column_order is not None else
            np.nansum(catalog_ref.loc[ref.extras.select, config.columns_ref_flux].values, axis=1)
        )
        order = np.argsort(column_order if config.order_ascending else -column_order)

        # Need the original reference row indices for output
        idx_orig_ref, idx_orig_target = (np.argwhere(cat.extras.select) for cat in (ref, target))

        # Retrieve required columns, including any converted ones (default to original column name)
        columns_convert = config.coord_format.coords_ref_to_convert
        if columns_convert is None:
            columns_convert = {}
        data_ref = catalog_ref.loc[
            ref.extras.indices[order],
            [columns_convert.get(column, column) for column in config.columns_ref_meas]
        ]
        data_target = catalog_target.loc[target.extras.select, config.columns_target_meas]
        errors_target = catalog_target.loc[target.extras.select, config.columns_target_err]

        exceptions = {}
        # The kdTree uses len(inputs) as a sentinel value for no match
        matched_target = {n_target_select, }

        t_begin = time.process_time()

        logger.info('Matching n_indices=%d/%d', len(order), len(catalog_ref))
        for index_n, index_row_select in enumerate(order):
            index_row = idx_orig_ref[index_row_select]
            ref_candidate_match[index_row] = True
            found = idxs_target_select[index_row_select, :]
            # Select match candidates from nearby sources not already matched
            # Note: set lookup is apparently fast enough that this is a few percent faster than:
            # found = [x for x in found[found != n_target_select] if x not in matched_target]
            # ... at least for ~1M sources
            found = [x for x in found if x not in matched_target]
            n_found = len(found)
            if n_found > 0:
                # This is an ndarray of n_found rows x len(data_ref/target) columns
                chi = (
                    (data_target.iloc[found].values - data_ref.iloc[index_n].values)
                    / errors_target.iloc[found].values
                )
                finite = np.isfinite(chi)
                n_finite = np.sum(finite, axis=1)
                # Require some number of finite chi_sq to match
                chisq_good = n_finite >= config.match_n_finite_min
                if np.any(chisq_good):
                    try:
                        chisq_sum = np.zeros(n_found, dtype=float)
                        chisq_sum[chisq_good] = np.nansum(chi[chisq_good, :] ** 2, axis=1)
                        idx_chisq_min = np.nanargmin(chisq_sum / n_finite)
                        ref_match_meas_finite[index_row] = n_finite[idx_chisq_min]
                        ref_match_count[index_row] = len(chisq_good)
                        ref_chisq[index_row] = chisq_sum[idx_chisq_min]
                        idx_match_select = found[idx_chisq_min]
                        row_target = target.extras.indices[idx_match_select]
                        ref_row_match[index_row] = row_target

                        target_row_match[row_target] = index_row
                        matched_target.add(idx_match_select)
                    except Exception as error:
                        # Can't foresee any exceptions, but they shouldn't prevent
                        # matching subsequent sources
                        exceptions[index_row] = error

            if logging_n_rows and ((index_n + 1) % logging_n_rows == 0):
                t_elapsed = time.process_time() - t_begin
                logger.info(
                    'Processed %d/%d in %.2fs at sort value=%.3f',
                    index_n + 1, n_ref_select, t_elapsed, column_order[order[index_n]],
                )

        catalog_out_ref = pd.DataFrame({
            'match_candidate': ref_candidate_match,
            'match_row': ref_row_match,
            'match_count': ref_match_count,
            'match_chisq': ref_chisq,
            'match_n_chisq_finite': ref_match_meas_finite,
        })

        catalog_out_target = pd.DataFrame({
            'match_candidate': target.extras.select if target.extras.select is not None else (
                np.ones(target.extras.n, dtype=bool)),
            'match_row': target_row_match,
        })

        return catalog_out_ref, catalog_out_target, exceptions
