// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2016 LSST/AURA
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include <sstream>

#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/astrom/SipTransform.h"

namespace lsst { namespace meas { namespace astrom {

void SipTransformBase::transformPixelsInPlace(afw::geom::AffineTransform const & s) {
    // The implementation for transformPixels is identical for
    // SipForwardTransform and SipReverseTransform.  That's pretty obvious for
    // the pixel origin and CD matrix, which are the same in both cases, but
    // it wasn't obvious to me until I did the math that the polynomial
    // transforms are composed with the affine transform the same way.
    auto sInv = s.invert();
    _pixelOrigin = s.getLinear()(_pixelOrigin - sInv.getTranslation());
    _cdMatrix = _cdMatrix * sInv.getLinear();
    _poly = compose(s.getLinear(), compose(getPoly(), sInv.getLinear()));
}

SipForwardTransform SipForwardTransform::convert(
    PolynomialTransform const & poly,
    afw::geom::Point2D const & pixelOrigin,
    afw::geom::LinearTransform const & cdMatrix
) {
    auto forwardSipPoly = compose(
        afw::geom::AffineTransform(cdMatrix.invert()),
        compose(
            poly,
            afw::geom::AffineTransform(afw::geom::Extent2D(pixelOrigin))
        )
    );
    // Subtracting 1 here accounts for the extra terms outside the sum in the
    // transform definition (see class docs) - note that you can fold those
    // terms into the sum by adding 1 from the A_10 and B_01 terms.
    forwardSipPoly._xCoeffs(1, 0) -= 1;
    forwardSipPoly._yCoeffs(0, 1) -= 1;
    return SipForwardTransform(pixelOrigin, cdMatrix, forwardSipPoly);
}

SipForwardTransform SipForwardTransform::convert(
    ScaledPolynomialTransform const & scaled,
    afw::geom::Point2D const & pixelOrigin,
    afw::geom::LinearTransform const & cdMatrix
) {
    auto forwardSipPoly = compose(
        afw::geom::AffineTransform(cdMatrix.invert())*scaled.getOutputScalingInverse(),
        compose(
            scaled.getPoly(),
            scaled.getInputScaling()*afw::geom::AffineTransform(afw::geom::Extent2D(pixelOrigin))
        )
    );
    // Account for the terms outside the sum in the definition (see comment
    // earlier in the file for more explanation).
    forwardSipPoly._xCoeffs(1, 0) -= 1;
    forwardSipPoly._yCoeffs(0, 1) -= 1;
    return SipForwardTransform(pixelOrigin, cdMatrix, forwardSipPoly);
}

SipForwardTransform SipForwardTransform::convert(ScaledPolynomialTransform const & scaled) {
    afw::geom::Point2D pixelOrigin(-scaled.getOutputScalingInverse().getTranslation());
    afw::geom::LinearTransform cdMatrix(scaled.getInputScaling().getLinear().invert());
    return convert(scaled, pixelOrigin, cdMatrix);
}

SipForwardTransform SipForwardTransform::extract(afw::image::TanWcs const & wcs) {
    if (!wcs.hasDistortion()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "Constructing a SipForwardTransform from a TanWcs with no distortions is not implemented."
        );
    }
    ndarray::Array<double,2,2> xCoeffs = ndarray::allocate(wcs.getSipA().rows(), wcs.getSipA().cols());
    ndarray::Array<double,2,2> yCoeffs = ndarray::allocate(wcs.getSipB().rows(), wcs.getSipB().cols());
    xCoeffs.asEigen() = wcs.getSipA();
    yCoeffs.asEigen() = wcs.getSipB();
    return SipForwardTransform(
        wcs.getPixelOrigin(),
        afw::geom::LinearTransform(wcs.getCDMatrix()),
        PolynomialTransform(xCoeffs, yCoeffs)
    );
}

afw::geom::AffineTransform SipForwardTransform::linearize(afw::geom::Point2D const & in) const {
    afw::geom::AffineTransform tail(-afw::geom::Extent2D(getPixelOrigin()));
    return afw::geom::AffineTransform(_cdMatrix)
        * (afw::geom::AffineTransform() + _poly.linearize(tail(in)))
        * tail;
}

afw::geom::Point2D SipForwardTransform::operator()(afw::geom::Point2D const & uv) const {
    afw::geom::Point2D duv(uv - afw::geom::Extent2D(getPixelOrigin()));
    return getCDMatrix()(afw::geom::Extent2D(duv) + getPoly()(duv));
}

SipForwardTransform SipForwardTransform::transformPixels(afw::geom::AffineTransform const & s) const {
    SipForwardTransform result(*this);
    result.transformPixelsInPlace(s);
    return result;
}

SipReverseTransform SipReverseTransform::convert(
    PolynomialTransform const & poly,
    afw::geom::Point2D const & pixelOrigin,
    afw::geom::LinearTransform const & cdMatrix
) {
    auto reverseSipPoly = compose(
        afw::geom::AffineTransform(-afw::geom::Extent2D(pixelOrigin)),
        compose(
            poly,
            afw::geom::AffineTransform(cdMatrix)
        )
    );
    // Account for the terms outside the sum in the definition (see comment
    // earlier in the file for more explanation).
    reverseSipPoly._xCoeffs(1, 0) -= 1;
    reverseSipPoly._yCoeffs(0, 1) -= 1;
    return SipReverseTransform(pixelOrigin, cdMatrix, reverseSipPoly);
}

SipReverseTransform SipReverseTransform::convert(
    ScaledPolynomialTransform const & scaled,
    afw::geom::Point2D const & pixelOrigin,
    afw::geom::LinearTransform const & cdMatrix
) {
    auto reverseSipPoly = compose(
        afw::geom::AffineTransform(-afw::geom::Extent2D(pixelOrigin))
        *scaled.getOutputScalingInverse(),
        compose(
            scaled.getPoly(),
            scaled.getInputScaling()*afw::geom::AffineTransform(cdMatrix)
        )
    );
    // Account for the terms outside the sum in the definition (see comment
    // earlier in the file for more explanation).
    reverseSipPoly._xCoeffs(1, 0) -= 1;
    reverseSipPoly._yCoeffs(0, 1) -= 1;
    return SipReverseTransform(pixelOrigin, cdMatrix, reverseSipPoly);
}

SipReverseTransform SipReverseTransform::convert(ScaledPolynomialTransform const & scaled) {
    return convert(
        scaled,
        afw::geom::Point2D(scaled.getOutputScalingInverse().getTranslation()),
        scaled.getInputScaling().getLinear()
    );
}

SipReverseTransform SipReverseTransform::extract(afw::image::TanWcs const & wcs) {
    if (!wcs.hasDistortion()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "Constructing a SipReverseTransform from a TanWcs with no distortions is not implemented."
        );
    }
    ndarray::Array<double,2,2> xCoeffs = ndarray::allocate(wcs.getSipAp().rows(), wcs.getSipAp().cols());
    ndarray::Array<double,2,2> yCoeffs = ndarray::allocate(wcs.getSipBp().rows(), wcs.getSipBp().cols());
    xCoeffs.asEigen() = wcs.getSipAp();
    yCoeffs.asEigen() = wcs.getSipBp();
    return SipReverseTransform(
        wcs.getPixelOrigin(),
        afw::geom::LinearTransform(wcs.getCDMatrix()),
        PolynomialTransform(xCoeffs, yCoeffs)
    );
}

SipReverseTransform SipReverseTransform::transformPixels(afw::geom::AffineTransform const & s) const {
    SipReverseTransform result(*this);
    result.transformPixelsInPlace(s);
    result._cdInverse = result._cdMatrix.invert();
    return result;
}

afw::geom::AffineTransform SipReverseTransform::linearize(afw::geom::Point2D const & in) const {
    return afw::geom::AffineTransform(afw::geom::Extent2D(getPixelOrigin()))
        * (afw::geom::AffineTransform() + _poly.linearize(_cdInverse(in)))
        * _cdInverse;
}

afw::geom::Point2D SipReverseTransform::operator()(afw::geom::Point2D const & xy) const {
    afw::geom::Point2D UV = _cdInverse(xy);
    return afw::geom::Extent2D(UV) + afw::geom::Extent2D(getPixelOrigin()) + getPoly()(UV);
}


PTR(afw::image::TanWcs) makeWcs(
    SipForwardTransform const & sipForward,
    SipReverseTransform const & sipReverse,
    afw::coord::Coord const & skyOrigin
) {
    if (!sipForward.getPixelOrigin().asEigen().isApprox(sipReverse.getPixelOrigin().asEigen())) {
        std::ostringstream oss;
        oss << "SIP forward and reverse transforms have inconsistent CRPIX: "
            << sipForward.getPixelOrigin() << " != " << sipReverse.getPixelOrigin();
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            oss.str()
        );
    }
    if (!sipForward.getCDMatrix().getMatrix().isApprox(sipReverse.getCDMatrix().getMatrix())) {
        std::ostringstream oss;
        oss << "SIP forward and reverse transforms have inconsistent CD matrix: "
            << sipForward.getCDMatrix() << "\n!=\n" << sipReverse.getCDMatrix();
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            oss.str()
        );
    }
    Eigen::MatrixXd sipA(sipForward.getPoly().getXCoeffs().asEigen());
    Eigen::MatrixXd sipB(sipForward.getPoly().getYCoeffs().asEigen());
    Eigen::MatrixXd sipAP(sipReverse.getPoly().getXCoeffs().asEigen());
    Eigen::MatrixXd sipBP(sipReverse.getPoly().getYCoeffs().asEigen());
    // TanWcs uses strings for coordinate systems, while Coord uses an enum.
    // Frustratingly, there's no way to convert from the enum to the string.
    std::string coordSys;
    switch (skyOrigin.getCoordSystem()) {
    case afw::coord::ICRS:
        coordSys = "ICRS";
        break;
    case afw::coord::FK5:
        coordSys = "FK5";
        break;
    default:
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "Coordinate system not supported"
        );
    }
    return std::make_shared<afw::image::TanWcs>(
        skyOrigin.getPosition(afw::geom::degrees),
        sipForward.getPixelOrigin(),
        sipForward.getCDMatrix().getMatrix(),
        sipA, sipB, sipAP, sipBP,
        skyOrigin.getEpoch(),
        coordSys
    );
}

std::shared_ptr<afw::image::TanWcs> transformWcsPixels(
    afw::image::TanWcs const & wcs,
    afw::geom::AffineTransform const & s
) {
    if (wcs.hasDistortion()) {
        auto fwd = SipForwardTransform::extract(wcs).transformPixels(s);
        auto rev = SipReverseTransform::extract(wcs).transformPixels(s);
        return makeWcs(fwd, rev, *wcs.getSkyOrigin());
    } else {
        auto sInv = s.invert();
        auto pixelOrigin = s.getLinear()(wcs.getPixelOrigin() - sInv.getTranslation());
        Eigen::Matrix2d cdMatrix = wcs.getCDMatrix() * sInv.getLinear().getMatrix();
        return std::make_shared<afw::image::TanWcs>(
            wcs.getSkyOrigin()->toIcrs().getPosition(afw::geom::degrees),
            pixelOrigin,
            cdMatrix,
            wcs.getEquinox(),
            "ICRS"
        );
    }
}

std::shared_ptr<afw::image::TanWcs> rotateWcsPixelsBy90(
    afw::image::TanWcs const & wcs,
    int nQuarter,
    afw::geom::Extent2I const & dimensions
) {
    afw::geom::Extent2D offset;
    switch(nQuarter % 4) {
    case 0:
        offset = afw::geom::Extent2D(0, 0);
        break;
    case 1:
        offset = afw::geom::Extent2D(dimensions.getY() - 1, 0);
        break;
    case 2:
        offset = afw::geom::Extent2D(dimensions - afw::geom::Extent2I(1, 1));
        break;
    case 3:
        offset = afw::geom::Extent2D(0, dimensions.getX() - 1);
        break;
    }
    auto rot = afw::geom::LinearTransform::makeRotation(nQuarter*90.0*afw::geom::degrees);
    return transformWcsPixels(
        wcs,
        afw::geom::AffineTransform(rot, offset)
    );
}

}}} // namespace lsst::meas::astrom
