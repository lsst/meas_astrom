#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
__all__ = ["FitSipDistortionTask", "FitSipDistortionConfig"]


import lsst.sphgeom
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.geom
import lsst.afw.display

from .scaledPolynomialTransformFitter import ScaledPolynomialTransformFitter, OutlierRejectionControl
from .sipTransform import SipForwardTransform, SipReverseTransform, makeWcs
from .makeMatchStatistics import makeMatchStatisticsInRadians

from .setMatchDistance import setMatchDistance


class FitSipDistortionConfig(lsst.pex.config.Config):
    """""Config for FitSipDistortionTask"""
    order = lsst.pex.config.RangeField(
        doc="Order of SIP polynomial",
        dtype=int,
        default=4,
        min=0,
    )
    numRejIter = lsst.pex.config.RangeField(
        doc="Number of rejection iterations",
        dtype=int,
        default=3,
        min=0,
    )
    rejSigma = lsst.pex.config.RangeField(
        doc="Number of standard deviations for clipping level",
        dtype=float,
        default=3.0,
        min=0.0,
    )
    nClipMin = lsst.pex.config.Field(
        doc="Minimum number of matches to reject when sigma-clipping",
        dtype=int,
        default=0
    )
    nClipMax = lsst.pex.config.Field(
        doc="Maximum number of matches to reject when sigma-clipping",
        dtype=int,
        default=1
    )
    maxScatterArcsec = lsst.pex.config.RangeField(
        doc="Maximum median scatter of a WCS fit beyond which the fit fails (arcsec); "
            "be generous, as this is only intended to catch catastrophic failures",
        dtype=float,
        default=10,
        min=0,
    )
    refUncertainty = lsst.pex.config.Field(
        doc="RMS uncertainty in reference catalog positions, in pixels.  Will be added "
            "in quadrature with measured uncertainties in the fit.",
        dtype=float,
        default=0.25,
    )
    nGridX = lsst.pex.config.Field(
        doc="Number of X grid points used to invert the SIP reverse transform.",
        dtype=int,
        default=100,
    )
    nGridY = lsst.pex.config.Field(
        doc="Number of Y grid points used to invert the SIP reverse transform.",
        dtype=int,
        default=100,
    )
    gridBorder = lsst.pex.config.Field(
        doc="When setting the gird region, how much to extend the image "
            "bounding box (in pixels) before transforming it to intermediate "
            "world coordinates using the initial WCS.",
        dtype=float,
        default=50.0,
    )


class FitSipDistortionTask(lsst.pipe.base.Task):
    """Fit a TAN-SIP WCS given a list of reference object/source matches.
    """
    ConfigClass = FitSipDistortionConfig
    _DefaultName = "fitWcs"

    def __init__(self, **kwargs):
        lsst.pipe.base.Task.__init__(self, **kwargs)
        self.outlierRejectionCtrl = OutlierRejectionControl()
        self.outlierRejectionCtrl.nClipMin = self.config.nClipMin
        self.outlierRejectionCtrl.nClipMax = self.config.nClipMax
        self.outlierRejectionCtrl.nSigma = self.config.rejSigma

    @lsst.pipe.base.timeMethod
    def fitWcs(self, matches, initWcs, bbox=None, refCat=None, sourceCat=None, exposure=None):
        """Fit a TAN-SIP WCS from a list of reference object/source matches.

        Parameters
        ----------
        matches : `list` of :cpp:class:`lsst::afw::table::ReferenceMatch`
            A sequence of reference object/source matches.
            The following fields are read:
            - match.first (reference object) coord
            - match.second (source) centroid

            The following fields are written:
            - match.first (reference object) centroid
            - match.second (source) centroid
            - match.distance (on sky separation, in radians)

        initWcs : :cpp:class:`lsst::afw::geom::SkyWcs`
            An initial WCS whose CD matrix is used as the final CD matrix.
        bbox : :cpp:class:`lsst::afw::geom::Box2I`
            The region over which the WCS will be valid (PARENT pixel coordinates);
            if `None` or an empty box then computed from matches
        refCat : :cpp:class:`lsst::afw::table::SimpleCatalog`
            Reference object catalog, or `None`.
            If provided then all centroids are updated with the new WCS,
            otherwise only the centroids for ref objects in matches are updated.
            Required fields are "centroid_x", "centroid_y", "coord_ra", and "coord_dec".
        sourceCat : :cpp:class:`lsst::afw::table::SourceCatalog`
            Source catalog, or `None`.
            If provided then coords are updated with the new WCS;
            otherwise only the coords for sources in matches are updated.
            Required input fields are "slot_Centroid_x", "slot_Centroid_y",
            "slot_Centroid_xErr", "slot_Centroid_yErr", and optionally
            "slot_Centroid_x_y_Cov".  The "coord_ra" and "coord_dec" fields
            will be updated but are not used as input.
        exposure : :cpp:class:`lsst::afw::image::Exposure`
            An Exposure or other displayable image on which matches can be
            overplotted.  Ignored (and may be `None`) if display-based debugging
            is not enabled via lsstDebug.

        Returns
        -------
        An lsst.pipe.base.Struct with the following fields:
            - wcs : :cpp:class:`lsst::afw::geom::SkyWcs`
                The best-fit WCS.
            - scatterOnSky : :cpp:class:`lsst::afw::geom::Angle`
                The median on-sky separation between reference objects and
                sources in "matches", as an lsst.afw.geom.Angle
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayFrame = lsstDebug.Info(__name__).frame
        displayPause = lsstDebug.Info(__name__).pause

        if bbox is None:
            bbox = lsst.afw.geom.Box2D()
            for match in matches:
                bbox.include(match.second.getCentroid())
            bbox = lsst.afw.geom.Box2I(bbox)

        wcs = self.makeInitialWcs(matches, initWcs)
        cdMatrix = lsst.afw.geom.LinearTransform(wcs.getCdMatrix())

        # Fit the "reverse" mapping from intermediate world coordinates to
        # pixels, rejecting outliers. Fitting in this direction first makes it
        # easier to handle the case where we have uncertainty on source
        # positions but not reference positions.  That's the case we have
        # right now for purely bookeeeping reasons, and it may be the case we
        # have in the future when we us Gaia as the reference catalog.
        revFitter = ScaledPolynomialTransformFitter.fromMatches(self.config.order, matches, wcs,
                                                                self.config.refUncertainty)
        revFitter.fit()
        for nIter in range(self.config.numRejIter):
            revFitter.updateModel()
            intrinsicScatter = revFitter.updateIntrinsicScatter()
            clippedSigma, nRejected = revFitter.rejectOutliers(self.outlierRejectionCtrl)
            self.log.debug(
                "Iteration {0}: intrinsic scatter is {1:4.3f} pixels, "
                "rejected {2} outliers at {3:3.2f} sigma.".format(
                    nIter+1, intrinsicScatter, nRejected, clippedSigma
                )
            )
            if display:
                displayFrame = self.display(revFitter, exposure=exposure, bbox=bbox,
                                            frame=displayFrame, displayPause=displayPause)
            revFitter.fit()
        revScaledPoly = revFitter.getTransform()
        # Convert the generic ScaledPolynomialTransform result to SIP form
        # with given CRPIX and CD (this is an exact conversion, up to
        # floating-point round-off error)
        sipReverse = SipReverseTransform.convert(revScaledPoly, wcs.getPixelOrigin(), cdMatrix)

        # Fit the forward mapping to a grid of points created from the reverse
        # transform.  Because that grid needs to be defined in intermediate
        # world coordinates, and we don't have a good way to get from pixels to
        # intermediate world coordinates yet (that's what we're fitting), we'll
        # first grow the box to make it conservatively large...
        gridBBoxPix = lsst.afw.geom.Box2D(bbox)
        gridBBoxPix.grow(self.config.gridBorder)
        # ...and then we'll transform using just the CRPIX offset and CD matrix
        # linear transform, which is the TAN-only (no SIP distortion, and
        # hence approximate) mapping from pixels to intermediate world
        # coordinates.
        gridBBoxIwc = lsst.afw.geom.Box2D()
        for point in gridBBoxPix.getCorners():
            point -= lsst.afw.geom.Extent2D(wcs.getPixelOrigin())
            gridBBoxIwc.include(cdMatrix(point))
        fwdFitter = ScaledPolynomialTransformFitter.fromGrid(self.config.order, gridBBoxIwc,
                                                             self.config.nGridX, self.config.nGridY,
                                                             revScaledPoly)
        fwdFitter.fit()
        # Convert to SIP forward form.
        fwdScaledPoly = fwdFitter.getTransform()
        sipForward = SipForwardTransform.convert(fwdScaledPoly, wcs.getPixelOrigin(), cdMatrix)

        # Make a new WCS from the SIP transform objects and the CRVAL in the
        # initial WCS.
        wcs = makeWcs(sipForward, sipReverse, wcs.getSkyOrigin())

        if refCat is not None:
            self.log.debug("Updating centroids in refCat")
            lsst.afw.table.updateRefCentroids(wcs, refList=refCat)
        else:
            self.log.warn("Updating reference object centroids in match list; refCat is None")
            lsst.afw.table.updateRefCentroids(wcs, refList=[match.first for match in matches])

        if sourceCat is not None:
            self.log.debug("Updating coords in sourceCat")
            lsst.afw.table.updateSourceCoords(wcs, sourceList=sourceCat)
        else:
            self.log.warn("Updating source coords in match list; sourceCat is None")
            lsst.afw.table.updateSourceCoords(wcs, sourceList=[match.second for match in matches])

        self.log.debug("Updating distance in match list")
        setMatchDistance(matches)

        stats = makeMatchStatisticsInRadians(wcs, matches, lsst.afw.math.MEDIAN)
        scatterOnSky = stats.getValue()*lsst.afw.geom.radians

        if scatterOnSky.asArcseconds() > self.config.maxScatterArcsec:
            raise lsst.pipe.base.TaskError(
                "Fit failed: median scatter on sky = %0.3f arcsec > %0.3f config.maxScatterArcsec" %
                (scatterOnSky.asArcseconds(), self.config.maxScatterArcsec))

        return lsst.pipe.base.Struct(
            wcs=wcs,
            scatterOnSky=scatterOnSky,
        )

    def display(self, revFitter, exposure=None, bbox=None, frame=0, pause=True):
        """Display positions and outlier status overlaid on an image.

        This method is called by fitWcs when display debugging is enabled.  It
        always drops into pdb before returning to allow interactive inspection,
        and hence it should never be called in non-interactive contexts.

        Parameters
        ----------
        revFitter : :cpp:class:`lsst::meas::astrom::ScaledPolynomialTransformFitter`
            Fitter object initialized with `fromMatches` for fitting a "reverse"
            distortion: the mapping from intermediate world coordinates to
            pixels.
        exposure : :cpp:class:`lsst::afw::image::Exposure`
            An Exposure or other displayable image on which matches can be
            overplotted.
        bbox : :cpp:class:`lsst::afw::geom::Box2I`
            Bounding box of the region on which matches should be plotted.
        """
        data = revFitter.getData()
        disp = lsst.afw.display.getDisplay(frame=frame)
        if exposure is not None:
            disp.mtv(exposure)
        elif bbox is not None:
            disp.mtv(exposure=lsst.afw.image.ExposureF(bbox))
        else:
            raise TypeError("At least one of 'exposure' and 'bbox' must be provided.")
        data = revFitter.getData()
        srcKey = lsst.afw.table.Point2DKey(data.schema["src"])
        srcErrKey = lsst.afw.table.CovarianceMatrix2fKey(data.schema["src"], ["x", "y"])
        refKey = lsst.afw.table.Point2DKey(data.schema["initial"])
        modelKey = lsst.afw.table.Point2DKey(data.schema["model"])
        rejectedKey = data.schema.find("rejected").key
        with disp.Buffering():
            for record in data:
                colors = ((lsst.afw.display.RED, lsst.afw.display.GREEN)
                          if not record.get(rejectedKey) else
                          (lsst.afw.display.MAGENTA, lsst.afw.display.CYAN))
                rx, ry = record.get(refKey)
                disp.dot("x", rx, ry, size=10, ctype=colors[0])
                mx, my = record.get(modelKey)
                disp.dot("o", mx, my, size=10, ctype=colors[0])
                disp.line([(rx, ry), (mx, my)], ctype=colors[0])
                sx, sy = record.get(srcKey)
                sErr = record.get(srcErrKey)
                sEllipse = lsst.afw.geom.Quadrupole(sErr[0, 0], sErr[1, 1], sErr[0, 1])
                disp.dot(sEllipse, sx, sy, ctype=colors[1])
        if pause or pause is None:  # default is to pause
            print("Dropping into debugger to allow inspection of display. Type 'continue' when done.")
            import pdb
            pdb.set_trace()
            return frame
        else:
            return frame + 1    # increment and return the frame for the next iteration.

    def makeInitialWcs(self, matches, wcs):
        """Generate a guess Wcs from the astrometric matches

        We create a Wcs anchored at the center of the matches, with the scale
        of the input Wcs.  This is necessary because the Wcs may have a very
        approximation position (as is common with telescoped-generated Wcs).
        We're using the best of each: positions from the matches, and scale
        from the input Wcs.

        Parameters
        ----------
        matches : list of :cpp:class:`lsst::afw::table::ReferenceMatch`
            A sequence of reference object/source matches.
            The following fields are read:

            - match.first (reference object) coord
            - match.second (source) centroid

        wcs : :cpp:class:`lsst::afw::geom::SkyWcs`
            An initial WCS whose CD matrix is used as the CD matrix of the
            result.

        Returns
        -------
        newWcsL : :cpp:class:`lsst::afw::geom::SkyWcs`
            A new WCS guess.
        """
        crpix = lsst.afw.geom.Extent2D(0, 0)
        crval = lsst.sphgeom.Vector3d(0, 0, 0)
        for mm in matches:
            crpix += lsst.afw.geom.Extent2D(mm.second.getCentroid())
            crval += mm.first.getCoord().getVector()
        crpix /= len(matches)
        crval /= len(matches)
        cd = wcs.getCdMatrix()
        newWcs = lsst.afw.geom.makeSkyWcs(crpix=lsst.afw.geom.Point2D(crpix),
                                          crval=lsst.afw.geom.SpherePoint(crval),
                                          cdMatrix=cd)
        return newWcs
