# This file is part of obs_test.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This file is ported from obs_test and should be removed when
# gen2 testing of meas_astrom is moved to gen3.
__all__ = ["TestMapper"]

import os
import warnings
import numpy as np

import lsst.utils
import lsst.geom as geom
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.image.utils as afwImageUtils
import lsst.daf.persistence as dafPersist
from lsst.obs.base import CameraMapper
# from .testCamera import TestCamera
# from .makeTestRawVisitInfo import MakeTestRawVisitInfo


class TestMapper(CameraMapper):
    """Camera mapper for the Test camera.
    """
    packageName = 'meas_astrom'

    def __init__(self, inputPolicy=None, **kwargs):
        policyFilePath = dafPersist.Policy.defaultPolicyFile(self.packageName, "testMapper.yaml", "policy")
        policy = dafPersist.Policy(policyFilePath)

        self.doFootprints = False
        if inputPolicy is not None:
            for kw in inputPolicy.paramNames(True):
                if kw == "doFootprints":
                    self.doFootprints = True
                else:
                    kwargs[kw] = inputPolicy.get(kw)

        CameraMapper.__init__(self, policy, policyFilePath, **kwargs)
        self.filterIdMap = {
            'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4, 'y': 5, 'i2': 5}

        with warnings.catch_warnings():
            # surpress Filter warnings; we already know this is deprecated
            warnings.simplefilter('ignore', category=FutureWarning)

            # The LSST Filters from L. Jones 04/07/10
            afwImageUtils.defineFilter('u', 364.59)
            afwImageUtils.defineFilter('g', 476.31)
            afwImageUtils.defineFilter('r', 619.42)
            afwImageUtils.defineFilter('i', 752.06)
            afwImageUtils.defineFilter('z', 866.85)
            afwImageUtils.defineFilter('y', 971.68, alias=['y4'])  # official y filter

    def _extractDetectorName(self, dataId):
        return "0"

    def _defectLookup(self, dataId):
        """Find the defects for a given CCD.

        Parameters
        ----------
        dataId : `dict`
            Dataset identifier

        Returns
        -------
        result : `str`
            Path to the defects file.

        Raises
        ------
        RuntimeError
            If ``obs_test`` is not setup.
        """
        obsTestDir = lsst.utils.getPackageDir('obs_test')

        return os.path.join(obsTestDir, "data", "input", "defects", "defects.fits")

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        Parameters
        ----------
        dataId : `dict`
            Data identifier with visit
        """
        visit = dataId['visit']
        return int(visit)

    def bypass_ccdExposureId(self, datasetType, pythonType, location, dataId):
        return self._computeCcdExposureId(dataId)

    def bypass_ccdExposureId_bits(self, datasetType, pythonType, location, dataId):
        return 41

    def validate(self, dataId):
        visit = dataId.get("visit")
        if visit is not None and not isinstance(visit, int):
            dataId["visit"] = int(visit)
        return dataId

    def _setCcdExposureId(self, propertyList, dataId):
        propertyList.set("Computed_ccdExposureId", self._computeCcdExposureId(dataId))
        return propertyList

    def _makeCamera(self, policy, repositoryDir):
        """Make a camera describing the camera geometry.

        Returns
        -------
        testCamera : `TestCamera`
            Test camera.
        """
        return TestCamera()


class TestCamera(cameraGeom.Camera):
    """A simple test Camera.

    Notes
    -----
    The camera has one ccd with name "0".
    That CCD has four amplifiers with names "00", "01", "10", and "11".

    The camera is modeled after a small portion of the LSST sim
    Summer 2012 camera: a single detector with four amplifiers,
    consisting of raft 2,2 sensor 0,0, half of channels 0,0 0,1 1,0 and 1,1
    (the half closest to the Y centerline).

    Note that the Summer 2012 camera has one very weird feature:
    the bias region (rawHOverscanBbox) is actually a prescan
    (it appears before the data pixels).

    A raw image has the sky in the same orientation on all amplifier
    subregions, so no amplifier subregions need their pixels to be flipped.

    Standard keys are:

    * ``amp``: amplifier number: one of 00, 01, 10, 11
    * ``ccd``: ccd number: always 0
    * ``visit``: exposure number; test data includes one exposure
        with visit=1
    """
    def __new__(cls):
        plateScale = geom.Angle(20, geom.arcseconds)  # plate scale, in angle on sky/mm
        # Radial distortion is modeled as a radial polynomial that converts from focal plane (in mm)
        # to field angle (in radians). Thus the coefficients are:
        # C0: always 0, for continuity at the center of the focal plane; units are rad
        # C1: 1/plateScale; units are rad/mm
        # C2: usually 0; units are rad/mm^2
        # C3: radial distortion; units are rad/mm^3
        radialCoeff = np.array([0.0, 1.0, 0.0, 0.925]) / plateScale.asRadians()
        fieldAngleToFocalPlane = afwGeom.makeRadialTransform(radialCoeff)
        focalPlaneToFieldAngle = fieldAngleToFocalPlane.inverted()

        camera = cameraGeom.Camera.Builder("test")
        cls._makeDetectors(camera, focalPlaneToFieldAngle)
        camera.setTransformFromFocalPlaneTo(cameraGeom.FIELD_ANGLE, focalPlaneToFieldAngle)
        return camera.finish()

    def __init__(self):
        pass

    @classmethod
    def _makeDetectors(cls, camera, focalPlaneToFieldAngle):
        """Make a list of detectors

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.camera.Builder`
            Camera to append detectors to.
        focalPlaneToFieldAngle : `lsst.afw.geom.TransformPoint2ToPoint2`
            Transform from ``FOCAL_PLANE`` to ``FIELD_ANGLE`` coordinates
            in the forward direction.
        """
        detectorConfigList = cls._makeDetectorConfigList()
        for detectorConfig in detectorConfigList:
            amplifiers = cls._makeAmplifierCatalog()
            detBuilder = cameraGeom.addDetectorBuilderFromConfig(
                camera,
                detectorConfig,
                amplifiers,
                focalPlaneToFieldAngle,
            )
            if detBuilder is None:
                raise RuntimeError("Could not add detector!")

    @classmethod
    def _makeDetectorConfigList(cls):
        """Make a list of detector configs

        Returns
        -------
        detectorConfigList : `list` of `lsst.afw.cameraGeom.DetectorConfig`
            List of detector configs.
        """
        # this camera has one detector that corresponds to a subregion of lsstSim detector R:2,2 S:0,0
        # with lower left corner 0, 1000 and dimensions 1018 x 2000
        # i.e. half of each of the following channels: 0,0, 0,1, 1,0 and 1,1
        detector0Config = cameraGeom.DetectorConfig()
        detector0Config.name = '0'
        detector0Config.id = 0
        detector0Config.serial = '0000011'
        detector0Config.detectorType = 0
        detector0Config.bbox_x0 = 0
        detector0Config.bbox_x1 = 1017
        detector0Config.bbox_y0 = 0
        detector0Config.bbox_y1 = 1999
        detector0Config.pixelSize_x = 0.01
        detector0Config.pixelSize_y = 0.01
        detector0Config.transformDict.nativeSys = 'Pixels'
        detector0Config.transformDict.transforms = None
        detector0Config.refpos_x = 2035.5
        detector0Config.refpos_y = 999.5
        detector0Config.offset_x = -42.26073
        detector0Config.offset_y = -42.21914
        detector0Config.transposeDetector = False
        detector0Config.pitchDeg = 0.0
        detector0Config.yawDeg = 90.0
        detector0Config.rollDeg = 0.0
        return [detector0Config]

    @classmethod
    def _makeAmplifierCatalog(cls):
        """Construct an amplifier info catalog

        Returns
        -------
        ampCatalog : `List` of `lsst.afw.cameraGeom.Amplifier.Builder()
            Amplifier information catalog.

        Notes
        -----
        The LSSTSim S12 amplifiers are unusual in that they start with 4 pixels
        of usable bias region (which is used to set rawHOverscanBbox, despite the name),
        followed by the data. There is no other underscan or overscan.
        """
        xDataExtent = 509  # trimmed
        yDataExtent = 1000
        xBiasExtent = 4
        xRawExtent = xDataExtent + xBiasExtent
        yRawExtent = yDataExtent
        readNoise = 3.975  # amplifier read noise, in e-
        saturationLevel = 65535
        linearityType = cameraGeom.NullLinearityType
        linearityCoeffs = [0, 0, 0, 0]

        ampCatalog = []
        for ampX in (0, 1):
            for ampY in (0, 1):
                # amplifier gain (e-/ADU) and read noiuse (ADU/pixel) from lsstSim raw data
                # note that obs_test amp <ampX><ampY> = lsstSim amp C<ampY>,<ampX> (axes are swapped)
                gain = {
                    (0, 0): 1.7741,     # C0,0
                    (0, 1): 1.65881,    # C1,0
                    (1, 0): 1.74151,    # C0,1
                    (1, 1): 1.67073,    # C1,1
                }[(ampX, ampY)]
                readNoise = {
                    (0, 0): 3.97531706217237,   # C0,0
                    (0, 1): 4.08263755342685,   # C1,0
                    (1, 0): 4.02753931932633,   # C0,1
                    (1, 1): 4.1890610691135,    # C1,1
                }[(ampX, ampY)]
                amplifier = cameraGeom.Amplifier.Builder()
                amplifier.setName("%d%d" % (ampX, ampY))
                amplifier.setBBox(geom.Box2I(
                    geom.Point2I(ampX * xDataExtent, ampY * yDataExtent),
                    geom.Extent2I(xDataExtent, yDataExtent),
                ))

                x0Raw = ampX * xRawExtent
                y0Raw = ampY * yRawExtent

                # bias region (which is prescan, in this case) is before the data
                readCorner = cameraGeom.ReadoutCorner.LL
                x0Bias = x0Raw
                x0Data = x0Bias + xBiasExtent

                amplifier.setRawBBox(geom.Box2I(
                    geom.Point2I(x0Raw, y0Raw),
                    geom.Extent2I(xRawExtent, yRawExtent),
                ))
                amplifier.setRawDataBBox(geom.Box2I(
                    geom.Point2I(x0Data, y0Raw),
                    geom.Extent2I(xDataExtent, yDataExtent),
                ))
                amplifier.setRawHorizontalOverscanBBox(geom.Box2I(
                    geom.Point2I(x0Bias, y0Raw),
                    geom.Extent2I(xBiasExtent, yRawExtent),
                ))
                amplifier.setRawXYOffset(geom.Extent2I(x0Raw, y0Raw))
                amplifier.setReadoutCorner(readCorner)
                amplifier.setGain(gain)
                amplifier.setReadNoise(readNoise)
                amplifier.setSaturation(saturationLevel)
                amplifier.setSuspectLevel(float("nan"))
                amplifier.setLinearityCoeffs([float(val) for val in linearityCoeffs])
                amplifier.setLinearityType(linearityType)
                # amplifier.setHasRawInfo(True)
                amplifier.setRawFlipX(False)
                amplifier.setRawFlipY(False)
                amplifier.setRawVerticalOverscanBBox(geom.Box2I())  # no vertical overscan
                amplifier.setRawPrescanBBox(geom.Box2I())  # no horizontal prescan
                ampCatalog.append(amplifier)
        return ampCatalog
