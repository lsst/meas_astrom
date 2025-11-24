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

from lsst.afw.geom import makeModifiedWcs, makeTransform
from lsst.afw.table import ExposureTable, ExposureCatalog
from lsst.geom import AffineTransform, LinearTransform, Box2D, SpherePoint, degrees
from lsst.meas.astrom.refit_pointing import RefitPointingConfig, RefitPointingTask
from lsst.obs.base.yamlCamera import makeCamera
from lsst.obs.base.utils import createInitialSkyWcsFromBoresight
from lsst.resources import ResourcePath


class RefitPointingTestCase(lsst.utils.tests.TestCase):
    """Tests for `RefitPointingTask`."""

    @classmethod
    def setUpClass(cls):
        with ResourcePath("resource://lsst.obs.base/test/dummycam.yaml").as_local() as camera_path:
            cls.camera = makeCamera(camera_path.ospath)
        assert len(cls.camera) == 2, "This is what the test expects."

    def test_config_validate(self):
        """Test that the default config options are valid."""
        config = RefitPointingConfig()
        config.validate()

    def test_init(self):
        """Test that the task can be constructed and that it adds the expected
        schema columns.
        """
        config = RefitPointingConfig()
        config.schema_prefix = "pre_"
        schema = ExposureTable.makeMinimalSchema()
        RefitPointingTask(config=config, schema=schema)
        self.assertIn("pre_wcs_detector_pointing_residual", schema)
        self.assertIn("pre_wcs_visit_pointing_residual", schema)
        self.assertIn("pre_wcs_detector_pointing_rejected", schema)

    def test_ideal_run(self):
        """Test running the task with only good detectors with perfect
        camera-geometry-only WCSs to recover the boresight and orientation with
        zero residuals.
        """
        schema = ExposureTable.makeMinimalSchema()
        task = RefitPointingTask(schema=schema)
        catalog = ExposureCatalog(schema)
        boresight = SpherePoint(45.0, 60.0, degrees)
        orientation = 30.0*degrees
        for detector in self.camera:
            record = catalog.addNew()
            record.setId(detector.getId())
            record.setBBox(detector.getBBox())
            record.setWcs(createInitialSkyWcsFromBoresight(boresight, orientation, detector))
        result = task.run(catalog=catalog, camera=self.camera)
        # Check that we recover the input pointing.
        self.assertAlmostEqual(result.regions.boresight_ra, boresight.getRa().asDegrees())
        self.assertAlmostEqual(result.regions.boresight_dec, boresight.getDec().asDegrees())
        self.assertAlmostEqual(result.regions.orientation, orientation.asDegrees())
        # Test that the residuals are zero, since the WCSs are idealized.
        self.assertAlmostEqual(catalog[0]["wcs_detector_pointing_residual"].asArcseconds(), 0.0)
        self.assertAlmostEqual(catalog[1]["wcs_detector_pointing_residual"].asArcseconds(), 0.0)
        self.assertAlmostEqual(catalog[0]["wcs_visit_pointing_residual"].asArcseconds(), 0.0)
        self.assertAlmostEqual(catalog[1]["wcs_visit_pointing_residual"].asArcseconds(), 0.0)
        self.assertFalse(catalog[0]["wcs_detector_pointing_rejected"])
        self.assertFalse(catalog[1]["wcs_detector_pointing_rejected"])
        # Test that the output regions include points inside the detectors.
        for record in catalog:
            pixel_bbox = Box2D(record.getBBox())
            pixel_bbox.grow(-1E-6)  # shrink to avoid floating-point equality issues
            for test_point in record.getWcs().pixelToSky(pixel_bbox.getCorners()):
                test_vec = test_point.getVector()
                self.assertTrue(result.regions.visit_region.contains(test_vec))
                self.assertTrue(result.regions.detector_regions[record.getId()].contains(test_vec))

    def test_rejection_run(self):
        """Test running the task with one good detector with a perfect
        camera-geometry-only WCS and one detector with a garbage WCS, to make
        sure the latter is rejected.
        """
        schema = ExposureTable.makeMinimalSchema()
        task = RefitPointingTask(schema=schema)
        catalog = ExposureCatalog(schema)
        boresight = SpherePoint(45.0, 60.0, degrees)
        orientation = 30.0*degrees
        good_id = None
        bad_id = None
        correct_wcs = None
        for detector in self.camera:
            record = catalog.addNew()
            record.setId(detector.getId())
            record.setBBox(detector.getBBox())
            if good_id is None:
                good_id = detector.getId()
                record.setWcs(createInitialSkyWcsFromBoresight(boresight, orientation, detector))
            else:
                bad_id = detector.getId()
                correct_wcs = createInitialSkyWcsFromBoresight(boresight, orientation, detector)
                bad_wcs = makeModifiedWcs(
                    makeTransform(AffineTransform(LinearTransform.makeScaling(1.5, 0.75))),
                    correct_wcs,
                    False
                )
                record.setWcs(bad_wcs)
        result = task.run(catalog=catalog, camera=self.camera)
        # Test that the residuals are large for the bad detector, and that we
        # rejected it.
        self.assertGreater(catalog[bad_id]["wcs_detector_pointing_residual"].asArcseconds(), 10.0)
        self.assertGreater(catalog[bad_id]["wcs_visit_pointing_residual"].asArcseconds(), 10.0)
        self.assertTrue(catalog[bad_id]["wcs_detector_pointing_rejected"])
        # Check that we recover the input pointing (one detector is
        # sufficient).
        self.assertAlmostEqual(result.regions.boresight_ra, boresight.getRa().asDegrees())
        self.assertAlmostEqual(result.regions.boresight_dec, boresight.getDec().asDegrees())
        self.assertAlmostEqual(result.regions.orientation, orientation.asDegrees())
        # Test that the residuals are zero for the good detector.
        self.assertAlmostEqual(catalog[good_id]["wcs_detector_pointing_residual"].asArcseconds(), 0.0)
        self.assertAlmostEqual(catalog[good_id]["wcs_visit_pointing_residual"].asArcseconds(), 0.0)
        self.assertFalse(catalog[good_id]["wcs_detector_pointing_rejected"])
        # Test that the output regions include points inside the detectors.
        for record in catalog:
            pixel_bbox = Box2D(record.getBBox())
            pixel_bbox.grow(-1E-6)  # shrink to avoid floating-point equality issues
            if record.getId() == good_id:
                wcs = record.getWcs()
            else:
                wcs = correct_wcs
                # Grow the box to test that we padded this less-trustworthy
                # region.
                pixel_bbox.grow(50.0)
            for test_point in wcs.pixelToSky(pixel_bbox.getCorners()):
                test_vec = test_point.getVector()
                self.assertTrue(result.regions.visit_region.contains(test_vec))
                self.assertTrue(result.regions.detector_regions[record.getId()].contains(test_vec))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
