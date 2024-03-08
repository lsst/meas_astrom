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

from lsst.meas.astrom import MatchProbabilisticConfig


class MatchProbabilisticConfigTestCase(lsst.utils.tests.TestCase):
    """MatchProbabilisticConfig test case."""
    def setUp(self):
        kwargs = dict(
            columns_ref_meas=["x", "y"],
            columns_target_meas=["x", "y"],
            columns_target_err=["xErr", "yErr"],
        )
        configs_bad = {"too_few_finite": MatchProbabilisticConfig(**kwargs)}
        self.config_good = MatchProbabilisticConfig(match_n_finite_min=2, **kwargs)
        kwargs["columns_target_meas"] = ["x"]
        configs_bad["too_few_target_meas"] = MatchProbabilisticConfig(**kwargs)
        kwargs["columns_target_meas"] = ["x", "y", "z"]
        configs_bad["too_many_target_meas"] = MatchProbabilisticConfig(**kwargs)
        kwargs["columns_target_meas"] = ["x", "y"]
        kwargs["columns_target_err"] = ["xErr"]
        configs_bad["too_few_target_err"] = MatchProbabilisticConfig(**kwargs)
        kwargs["columns_target_err"] = ["xErr", "yErr", "zErr"]
        configs_bad["too_many_target_err"] = MatchProbabilisticConfig(**kwargs)
        self.configs_bad = configs_bad

    def tearDown(self):
        del self.config_good
        del self.configs_bad

    def test_MatchProbabilisticTask(self):
        for name, config in self.configs_bad.items():
            with self.assertRaises(ValueError, msg=f"expected {name} failure"):
                config.validate()
        self.config_good.validate()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
