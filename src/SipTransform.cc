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

#include "lsst/geom/Point.h"
#include "lsst/geom/Angle.h"
#include "lsst/geom/SpherePoint.h"
#include "lsst/geom/AffineTransform.h"
#include "lsst/geom/LinearTransform.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/geom/transformFactory.h"
#include "lsst/meas/astrom/SipTransform.h"

namespace lsst {
namespace meas {
namespace astrom {

void SipTransformBase::transformPixelsInPlace(geom::AffineTransform const& s) {
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

SipForwardTransform SipForwardTransform::convert(PolynomialTransform const& poly,
                                                 geom::Point2D const& pixelOrigin,
                                                 geom::LinearTransform const& cdMatrix) {
    auto forwardSipPoly = compose(geom::AffineTransform(cdMatrix.invert()),
                                  compose(poly, geom::AffineTransform(geom::Extent2D(pixelOrigin))));
    // Subtracting 1 here accounts for the extra terms outside the sum in the
    // transform definition (see class docs) - note that you can fold those
    // terms into the sum by adding 1 from the A_10 and B_01 terms.
    forwardSipPoly._xCoeffs(1, 0) -= 1;
    forwardSipPoly._yCoeffs(0, 1) -= 1;
    return SipForwardTransform(pixelOrigin, cdMatrix, forwardSipPoly);
}

SipForwardTransform SipForwardTransform::convert(ScaledPolynomialTransform const& scaled,
                                                 geom::Point2D const& pixelOrigin,
                                                 geom::LinearTransform const& cdMatrix) {
    auto forwardSipPoly =
            compose(geom::AffineTransform(cdMatrix.invert()) * scaled.getOutputScalingInverse(),
                    compose(scaled.getPoly(),
                            scaled.getInputScaling() * geom::AffineTransform(geom::Extent2D(pixelOrigin))));
    // Account for the terms outside the sum in the definition (see comment
    // earlier in the file for more explanation).
    forwardSipPoly._xCoeffs(1, 0) -= 1;
    forwardSipPoly._yCoeffs(0, 1) -= 1;
    return SipForwardTransform(pixelOrigin, cdMatrix, forwardSipPoly);
}

SipForwardTransform SipForwardTransform::convert(ScaledPolynomialTransform const& scaled) {
    geom::Point2D pixelOrigin(-scaled.getOutputScalingInverse().getTranslation());
    geom::LinearTransform cdMatrix(scaled.getInputScaling().getLinear().invert());
    return convert(scaled, pixelOrigin, cdMatrix);
}

geom::AffineTransform SipForwardTransform::linearize(geom::Point2D const& in) const {
    geom::AffineTransform tail(-geom::Extent2D(getPixelOrigin()));
    return geom::AffineTransform(_cdMatrix) * (geom::AffineTransform() + _poly.linearize(tail(in))) * tail;
}

geom::Point2D SipForwardTransform::operator()(geom::Point2D const& uv) const {
    geom::Point2D duv(uv - geom::Extent2D(getPixelOrigin()));
    return getCdMatrix()(geom::Extent2D(duv) + getPoly()(duv));
}

SipForwardTransform SipForwardTransform::transformPixels(geom::AffineTransform const& s) const {
    SipForwardTransform result(*this);
    result.transformPixelsInPlace(s);
    return result;
}

SipReverseTransform SipReverseTransform::convert(PolynomialTransform const& poly,
                                                 geom::Point2D const& pixelOrigin,
                                                 geom::LinearTransform const& cdMatrix) {
    auto reverseSipPoly = compose(geom::AffineTransform(-geom::Extent2D(pixelOrigin)),
                                  compose(poly, geom::AffineTransform(cdMatrix)));
    // Account for the terms outside the sum in the definition (see comment
    // earlier in the file for more explanation).
    reverseSipPoly._xCoeffs(1, 0) -= 1;
    reverseSipPoly._yCoeffs(0, 1) -= 1;
    return SipReverseTransform(pixelOrigin, cdMatrix, reverseSipPoly);
}

SipReverseTransform SipReverseTransform::convert(ScaledPolynomialTransform const& scaled,
                                                 geom::Point2D const& pixelOrigin,
                                                 geom::LinearTransform const& cdMatrix) {
    auto reverseSipPoly =
            compose(geom::AffineTransform(-geom::Extent2D(pixelOrigin)) * scaled.getOutputScalingInverse(),
                    compose(scaled.getPoly(), scaled.getInputScaling() * geom::AffineTransform(cdMatrix)));
    // Account for the terms outside the sum in the definition (see comment
    // earlier in the file for more explanation).
    reverseSipPoly._xCoeffs(1, 0) -= 1;
    reverseSipPoly._yCoeffs(0, 1) -= 1;
    return SipReverseTransform(pixelOrigin, cdMatrix, reverseSipPoly);
}

SipReverseTransform SipReverseTransform::convert(ScaledPolynomialTransform const& scaled) {
    return convert(scaled, geom::Point2D(scaled.getOutputScalingInverse().getTranslation()),
                   scaled.getInputScaling().getLinear());
}

SipReverseTransform SipReverseTransform::transformPixels(geom::AffineTransform const& s) const {
    SipReverseTransform result(*this);
    result.transformPixelsInPlace(s);
    result._cdInverse = result._cdMatrix.invert();
    return result;
}

geom::AffineTransform SipReverseTransform::linearize(geom::Point2D const& in) const {
    return geom::AffineTransform(geom::Extent2D(getPixelOrigin())) *
           (geom::AffineTransform() + _poly.linearize(_cdInverse(in))) * _cdInverse;
}

geom::Point2D SipReverseTransform::operator()(geom::Point2D const& xy) const {
    geom::Point2D UV = _cdInverse(xy);
    return geom::Extent2D(UV) + geom::Extent2D(getPixelOrigin()) + getPoly()(UV);
}

std::shared_ptr<afw::geom::SkyWcs> makeWcs(SipForwardTransform const& sipForward,
                                           SipReverseTransform const& sipReverse,
                                           geom::SpherePoint const& skyOrigin) {
    if (!sipForward.getPixelOrigin().asEigen().isApprox(sipReverse.getPixelOrigin().asEigen())) {
        std::ostringstream oss;
        oss << "SIP forward and reverse transforms have inconsistent CRPIX: " << sipForward.getPixelOrigin()
            << " != " << sipReverse.getPixelOrigin();
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, oss.str());
    }
    if (!sipForward.getCdMatrix().getMatrix().isApprox(sipReverse.getCdMatrix().getMatrix())) {
        std::ostringstream oss;
        oss << "SIP forward and reverse transforms have inconsistent CD matrix: " << sipForward.getCdMatrix()
            << "\n!=\n"
            << sipReverse.getCdMatrix();
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, oss.str());
    }
    Eigen::MatrixXd sipA(ndarray::asEigenMatrix(sipForward.getPoly().getXCoeffs()));
    Eigen::MatrixXd sipB(ndarray::asEigenMatrix(sipForward.getPoly().getYCoeffs()));
    Eigen::MatrixXd sipAP(ndarray::asEigenMatrix(sipReverse.getPoly().getXCoeffs()));
    Eigen::MatrixXd sipBP(ndarray::asEigenMatrix(sipReverse.getPoly().getYCoeffs()));

    return afw::geom::makeTanSipWcs(sipForward.getPixelOrigin(), skyOrigin,
                                    sipForward.getCdMatrix().getMatrix(), sipA, sipB, sipAP, sipBP);
}

std::shared_ptr<afw::geom::SkyWcs> transformWcsPixels(afw::geom::SkyWcs const& wcs,
                                                      geom::AffineTransform const& s) {
    auto affineTransform22 = afw::geom::makeTransform(s);
    return afw::geom::makeModifiedWcs(*affineTransform22->getInverse(), wcs, true);
}

std::shared_ptr<afw::geom::SkyWcs> rotateWcsPixelsBy90(afw::geom::SkyWcs const& wcs, int nQuarter,
                                                       geom::Extent2I const& dimensions) {
    geom::Extent2D offset;
    switch (nQuarter % 4) {
        case 0:
            offset = geom::Extent2D(0, 0);
            break;
        case 1:
            offset = geom::Extent2D(dimensions.getY() - 1, 0);
            break;
        case 2:
            offset = geom::Extent2D(dimensions - geom::Extent2I(1, 1));
            break;
        case 3:
            offset = geom::Extent2D(0, dimensions.getX() - 1);
            break;
    }
    auto rot = geom::LinearTransform::makeRotation(nQuarter * 90.0 * geom::degrees);
    return transformWcsPixels(wcs, geom::AffineTransform(rot, offset));
}

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
