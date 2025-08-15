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

__all__ = ()

import math
from typing import ClassVar

import numpy as np

from lsst.afw.geom import SkyWcs, getIntermediateWorldCoordsToSky, SipApproximation, makeTanSipWcs
from lsst.afw.cameraGeom import Detector, FIELD_ANGLE, PIXELS
from lsst.geom import (
    Angle,
    Box2D,
    Extent2I,
    Point2D,
    SpherePoint,
    SphereTransform,
    degrees,
    radians,
)
from lsst.pex.config import Config, Field
from lsst.pipe.base import Struct, Task
from lsst.obs.base import Instrument


class FitSipApproximationConfig(Config):
    pointing_grid_spacing = Field[float](
        doc=(
            "Spacing (in pixels) between grid points used to evaluate the WCS when fitting the pointing "
            "This can be a very sparse grid (there are only three degrees of freedom in this fit).  "
            "If the spacing does not divide the detector bounding box evenly, it is decreased slightly."
        ),
        dtype=float,
        default=512.0,
    )
    sip_grid_spacing = Field[float](
        doc=(
            "Spacing (in pixels) between grid points used to evaluate the WCS when fitting the "
            "SIP approximation.  This needs to be a fine grid to fit the number of degrees of freedom in "
            "the configured polynomial order. "
            "If the spacing does not divide the detector bounding box evenly, it is decreased slightly."
        ),
        dtype=float,
        default=32.0,
    )
    sip_order = Field[int](
        doc="Polynomial order for the SIP approximation.",
        dtype=int,
        default=5,
    )


class FitSipApproximationTask(Task):
    _DefaultName: ClassVar[str] = "fitSipApproximation"
    ConfigClass: ClassVar[type[FitSipApproximationConfig]] = FitSipApproximationConfig
    config: FitSipApproximationConfig

    def run(self, *, true_wcs: SkyWcs, detector: Detector, instrument: Instrument) -> Struct:
        boresight, orientation = self.fit_pointing(
            true_wcs=true_wcs, detector=detector, instrument=instrument
        )
        repointed_raw_wcs, _ = self.make_raw_wcs(
            boresight, orientation, detector=detector, instrument=instrument
        )
        sip_approx = self.fit_sip_approximation(
            true_wcs=true_wcs, repointed_raw_wcs=repointed_raw_wcs, detector=detector
        )
        fits_wcs = makeTanSipWcs(
            sip_approx.getPixelOrigin(),
            repointed_raw_wcs.getSkyOrigin(),
            sip_approx.getCdMatrix(),
            sip_approx.getA(),
            sip_approx.getB(),
            sip_approx.getAP(),
            sip_approx.getBP(),
        )
        return Struct(
            wcs=true_wcs.withFitsApproximation(fits_wcs),
            boresight=boresight,
            orientation=orientation,
            repointed_raw_wcs=repointed_raw_wcs,
        )

    def fit_pointing(
        self, *, true_wcs: SkyWcs, detector: Detector, instrument: Instrument
    ) -> tuple[SpherePoint, Angle]:
        start = true_wcs.pixelToSky(Point2D(0.0, 0.0))
        unpointed_raw_wcs, flip_x = self.make_raw_wcs(
            start, 0.0 * degrees, detector=detector, instrument=instrument
        )
        pixel_x, pixel_y = self._make_pointing_grid(detector)
        start_ra, start_dec = unpointed_raw_wcs.pixelToSkyArray(pixel_x, pixel_y)
        start_uv = np.stack(
            SpherePoint.toUnitXYZ(longitude=start_ra, latitude=start_dec, units=radians),
            axis=1,
        )
        true_ra, true_dec = true_wcs.pixelToSkyArray(pixel_x, pixel_y)
        true_uv = np.stack(
            SpherePoint.toUnitXYZ(longitude=true_ra, latitude=true_dec, units=radians),
            axis=1,
        )
        unpointed_y_axis_point = unpointed_raw_wcs.pixelToSky(
            detector.transform(Point2D(0.0, np.pi / 180.0), FIELD_ANGLE, PIXELS)
        )
        transform = SphereTransform.fit_unit_vectors(start_uv, true_uv)
        boresight = transform(start)
        transformed_y_axis_point = transform(unpointed_y_axis_point)
        orientation = Angle(90, degrees) - boresight.bearingTo(transformed_y_axis_point)
        if flip_x:
            raise NotImplementedError("flip_x=True is not implemented yet")
        return boresight, orientation

    def fit_sip_approximation(
        self, *, true_wcs: SkyWcs, repointed_raw_wcs: SkyWcs, detector: Detector
    ) -> SkyWcs:
        pixels_to_iwc = true_wcs.getTransform().then(
            getIntermediateWorldCoordsToSky(repointed_raw_wcs, simplify=True).inverted()
        )
        pixel_bbox = Box2D(detector.getBBox())
        grid_shape = Extent2I(
            # We add one to the dimensions since there's a point the min and
            # max in each dimension.
            np.ceil(pixel_bbox.width / self.config.sip_grid_spacing) + 1,
            np.ceil(pixel_bbox.height / self.config.sip_grid_spacing) + 1,
        )
        sip_approx = SipApproximation(
            pixels_to_iwc,
            repointed_raw_wcs.getPixelOrigin(),
            repointed_raw_wcs.getCdMatrix(),
            pixel_bbox,
            grid_shape,
            self.config.sip_order,
        )
        return sip_approx

    def make_raw_wcs(
        self, boresight: SpherePoint, orientation: Angle, *, detector: Detector, instrument: Instrument
    ) -> tuple[SkyWcs, bool]:
        raw_formatter = instrument.getRawFormatter({"detector": detector.getId()})
        return (
            raw_formatter.makeRawSkyWcsFromBoresight(boresight, orientation, detector),
            raw_formatter.wcsFlipX,
        )

    def _make_pointing_grid(self, detector: Detector) -> tuple[np.ndarray, np.ndarray]:
        pixel_bbox = Box2D(detector.getBBox())
        n_x = math.ceil(pixel_bbox.width / self.config.pointing_grid_spacing)
        n_y = math.ceil(pixel_bbox.height / self.config.pointing_grid_spacing)
        x, y = np.meshgrid(
            # We add one to the dimensions since there's a point the min and
            # max in each dimension.
            np.linspace(pixel_bbox.x.min, pixel_bbox.x.max, n_x + 1),
            np.linspace(pixel_bbox.y.min, pixel_bbox.y.max, n_y + 1),
        )
        return x.ravel(), y.ravel()


import astropy.table
import click
import asyncio
import tqdm
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from lsst.daf.butler import Butler, DatasetRef, DatasetProvenance, DataCoordinate
from lsst.obs.base import Instrument


class Runner:

    def __init__(self, root: str, collections: str, run: str):
        self.butler = Butler.from_config(root, collections=collections, run=run)
        self.task = FitSipApproximationTask()
        self.instrument = Instrument.fromName("LSSTComCam", self.butler.registry)
        self.camera = self.butler.get("camera", instrument=self.instrument.getName())

    instance: ClassVar[Runner]

    @staticmethod
    def initialize(root: str, collections: str, run: str) -> None:
        Runner.instance = Runner(root, collections, run)

    def run(self, *refs: DatasetRef) -> tuple[DataCoordinate, float, float]:
        approx_wcs = None
        for ref in refs:
            exposure = self.butler.get(ref)
            exposure.metadata["BUNIT"] = "nJy"
            if approx_wcs is None:
                detector_id = ref.dataId["detector"]
                results = self.task.run(true_wcs=exposure.wcs, detector=self.camera[detector_id], instrument=self.instrument)
                bbox = self.camera[detector_id].getBBox()
                x, y = np.meshgrid(np.linspace(bbox.x.min, bbox.x.max, 20), np.linspace(bbox.y.min, bbox.y.max, 20))
                true_ra, true_dec = results.wcs.pixelToSkyArray(x.ravel(), y.ravel())
                approx_ra, approx_dec = results.wcs.getFitsApproximation().pixelToSkyArray(x.ravel(), y.ravel())
                x1, y1 = results.wcs.getFitsApproximation().skyToPixelArray(true_ra, true_dec)
                x2, y2 = results.wcs.skyToPixelArray(approx_ra, approx_dec)
                delta1 = ((x.ravel() - x1)**2 + (y.ravel() - y1)**2)**0.5
                delta2 = ((x.ravel() - x2)**2 + (y.ravel() - y2)**2)**0.5
                approx_wcs = results.wcs
            exposure.setWcs(approx_wcs)
            provenance, _ = DatasetProvenance.from_flat_dict(exposure.metadata, self.butler)
            self.butler.put(exposure, ref.datasetType, ref.dataId, provenance=provenance)
        return ref.dataId, float(delta1.max()), float(delta2.max())

    @staticmethod
    def run_in_pool(refs: list[DatasetRef]) -> tuple[DataCoordinate, float, float]:
        return Runner.instance.run(*refs)

    @staticmethod
    async def query_and_fix(*, root: str, collections: str, run: str, jobs: int) -> None:
        butler = Butler.from_config(root, collections=collections, run=run)
        visit_image_refs = {ref.dataId: ref for ref in butler.query_datasets("visit_image", limit=None)}
        for ref in butler.query_datasets("visit_image", limit=None, collections=run):
            del visit_image_refs[ref.dataId]
        difference_image_refs = {ref.dataId: ref for ref in butler.query_datasets("difference_image", limit=None)}
        for ref in butler.query_datasets("difference_image", limit=None, collections=run):
            del difference_image_refs[ref.dataId]
        pairs = {
            data_id: ([visit_image_refs[data_id]] if data_id in visit_image_refs else [])
                + ([difference_image_refs[data_id]] if data_id in difference_image_refs else [])
            for data_id in visit_image_refs.keys() | difference_image_refs.keys()
        }
        with ProcessPoolExecutor(
            max_workers=jobs,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=Runner.initialize,
            initargs=(root, collections, run),
        ) as executor:
            loop = asyncio.get_running_loop()
            work: set[asyncio.Future[tuple[DataCoordinate, float, float]]] = set()
            for pair in pairs.values():
                work.add(loop.run_in_executor(executor, Runner.run_in_pool, pair))
            stats = []
            with tqdm.tqdm(asyncio.as_completed(work), total=len(pairs)) as pbar:
                for f in pbar:
                    data_id, delta1, delta2 = await f
                    if not delta1 < 0.01 or not delta2 < 0.01:
                        pbar.write(f"{data_id}: {delta1:0.5f} {delta2:0.5f}")
                    stats.append(data_id.required_values + (delta1, delta2))
        t = astropy.table.Table(rows=stats, names=list(data_id.dimensions.required) + ["delta1", "delta2"])
        t.write("stats.ecsv", overwrite=True)

@click.group()
def main() -> None:
    pass

@main.command()
@click.argument("butler")
@click.argument("collections")
@click.argument("run")
@click.option("-j", "--jobs", default=1, type=int)
def fix_visit_images(
    *,
    butler: str,
    collections: str,
    run: str,
    jobs: int = 1,
) -> None:
    asyncio.run(Runner.query_and_fix(root=butler, collections=collections, run=run, jobs=jobs))


if __name__ == "__main__":
    main()
