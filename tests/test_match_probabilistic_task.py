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


import unittest
import lsst.utils.tests

import lsst.afw.geom as afwGeom
from lsst.meas.astrom import ConvertCatalogCoordinatesConfig, MatchProbabilisticConfig, MatchProbabilisticTask

import numpy as np
import pandas as pd


class MatchProbabilisticTaskTestCase(lsst.utils.tests.TestCase):
    """MatchProbabilisticTask test case."""
    def setUp(self):
        ra = np.array([-0.1, -0.2, 0., 0.1, 0.2])
        dec = np.array([-0.15, 0.15, 0, 0.15, -0.15])
        mag_g = np.array([23., 24., 25., 25.5, 26.])
        mag_r = mag_g + [0.5, -0.2, -0.8, -0.5, -1.5]
        coord_format = ConvertCatalogCoordinatesConfig
        zeropoint = coord_format.mag_zeropoint_ref.default
        fluxes = tuple(-0.4*10**(mag - zeropoint) for mag in (mag_g, mag_r))
        eps_coord = np.full_like(ra, lsst.geom.Angle(0.2, lsst.geom.arcseconds).asDegrees())
        eps_flux = np.full_like(eps_coord, 10)
        flags = np.ones_like(eps_coord, dtype=bool)
        name_index = 'index'

        columns_flux = ['flux_g', 'flux_r']
        columns_ref_meas = [
            coord_format.column_ref_coord1.default,
            coord_format.column_ref_coord2.default,
        ] + columns_flux

        data_ref = {
            name_index: np.arange(len(ra)),
            columns_ref_meas[0]: ra[::-1],
            columns_ref_meas[1]: dec[::-1],
            columns_flux[0]: fluxes[0][::-1],
            columns_flux[1]: fluxes[1][::-1],
        }
        self.catalog_ref = pd.DataFrame(data=data_ref)

        columns_target_meas = [
            coord_format.column_target_coord1.default,
            coord_format.column_target_coord2.default,
        ] + columns_flux
        columns_target_err = [f'{column}Err' for column in columns_target_meas]

        data_target = {
            name_index: np.arange(len(ra)),
            columns_target_meas[0]: ra + eps_coord,
            columns_target_meas[1]: dec + eps_coord,
            f'{columns_target_meas[0]}Err': eps_coord,
            f'{columns_target_meas[1]}Err': eps_coord,
            columns_flux[0]: fluxes[0] + eps_flux,
            columns_flux[1]: fluxes[1] - eps_flux,
            f'{columns_flux[0]}Err': eps_flux,
            f'{columns_flux[1]}Err': eps_flux,
            MatchProbabilisticConfig.columns_target_select_true.default[0]: flags,
            MatchProbabilisticConfig.columns_target_select_false.default[0]: ~flags,
        }
        self.catalog_target = pd.DataFrame(data=data_target)

        self.task = MatchProbabilisticTask(config=MatchProbabilisticConfig(
            columns_ref_flux=columns_flux,
            columns_ref_meas=columns_ref_meas,
            columns_ref_copy=[name_index],
            columns_target_meas=columns_target_meas,
            columns_target_err=columns_target_err,
            columns_target_copy=[name_index],
        ))
        self.wcs = afwGeom.makeSkyWcs(crpix=lsst.geom.Point2D(9000, 9000),
                                      crval=lsst.geom.SpherePoint(180., 0., lsst.geom.degrees),
                                      cdMatrix=afwGeom.makeCdMatrix(scale=0.2*lsst.geom.arcseconds))

    def tearDown(self):
        del self.catalog_ref
        del self.catalog_target
        del self.task
        del self.wcs

    def test_MatchProbabilisticTask(self):
        for (columns, catalog) in (
            (self.task.columns_in_ref, self.catalog_ref),
            (self.task.columns_in_target, self.catalog_target),
        ):
            self.assertTrue(all((column in catalog.columns for column in columns)))
        result = self.task.run(
            catalog_ref=self.catalog_ref,
            catalog_target=self.catalog_target,
            wcs=self.wcs,
            logging_n_rows=2,
        )
        indices_target = list(result.cat_output_target.match_row.values)
        self.assertEqual(indices_target, list(np.arange(len(indices_target))[::-1]))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
