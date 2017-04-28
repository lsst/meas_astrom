from __future__ import absolute_import, division, print_function

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
