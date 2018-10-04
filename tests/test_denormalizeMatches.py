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

import sys
import unittest

import lsst.afw.geom
import lsst.afw.table

from lsst.meas.astrom import denormalizeMatches


class DenormalizeMatchesTestCase(unittest.TestCase):
    """Test the behaviour of the denormalizedMatches function"""

    def checkDenormalizeMatches(self, refType, srcType, MatchClass, num=10):
        """Check that denormalizeMatches works

        We create reference and source catalogs, generate matches,
        run denormalizeMatches and verify that the results are as expected.

        Parameters
        ----------
        refType : `str`
            Type of reference catalog/table; "Simple" or "Source".
        srcType : `str`
            Type of source catalog/table; "Simple" or "Source".
        MatchClass : `type`
            Class for match; should be suitable for the refType and srcType.
        """
        refSchema = getattr(lsst.afw.table, refType + "Table").makeMinimalSchema()
        refCat = getattr(lsst.afw.table, refType + "Catalog")(refSchema)
        for ii in range(num):
            ref = refCat.addNew()
            ref.set("id", ii)

        srcSchema = getattr(lsst.afw.table, srcType + "Table").makeMinimalSchema()
        srcCat = getattr(lsst.afw.table, srcType + "Catalog")(srcSchema)
        for ii in range(2*num, num, -1):
            src = srcCat.addNew()
            src.set("id", ii)

        matches = [MatchClass(ref, src, ref.get("id")) for ref, src in zip(refCat, srcCat)]
        catalog = denormalizeMatches(matches)
        for row, ref, src in zip(catalog, refCat, srcCat):
            self.assertEqual(row.get("ref_id"), ref.get("id"))
            self.assertEqual(row.get("src_id"), src.get("id"))
            self.assertEqual(row.get("distance"), ref.get("id"))

    def testDenormalizeMatches(self):
        """Test denormalizeMatches for various types"""
        for args in (("Simple", "Simple", lsst.afw.table.SimpleMatch),
                     ("Simple", "Source", lsst.afw.table.ReferenceMatch),
                     ("Source", "Source", lsst.afw.table.SourceMatch),
                     ):
            self.checkDenormalizeMatches(*args)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules[__name__])
    unittest.main()
