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

#include "lsst/geom/Extent.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/AffineTransform.h"
#include "lsst/meas/astrom/PolynomialTransform.h"
#include "lsst/meas/astrom/SipTransform.h"
#include "lsst/meas/astrom/detail/polynomialUtils.h"

namespace lsst {
namespace meas {
namespace astrom {

PolynomialTransform PolynomialTransform::convert(ScaledPolynomialTransform const& scaled) {
    return compose(scaled.getOutputScalingInverse(), compose(scaled.getPoly(), scaled.getInputScaling()));
}

PolynomialTransform PolynomialTransform::convert(SipForwardTransform const& other) {
    PolynomialTransform poly = other.getPoly();
    // Adding 1 here accounts for the extra terms outside the sum in the SIP
    // transform definition (see SipForwardTransform docs) - note that you can
    // fold those terms into the sum by adding 1 from the A_10 and B_01 terms.
    poly._xCoeffs(1, 0) += 1;
    poly._yCoeffs(0, 1) += 1;
    return compose(other.getCdMatrix(),
                   compose(poly, geom::AffineTransform(geom::Point2D() - other.getPixelOrigin())));
}

PolynomialTransform PolynomialTransform::convert(SipReverseTransform const& other) {
    PolynomialTransform poly = other.getPoly();
    // Account for the terms outside the sum in the SIP definition (see comment
    // earlier in the file for more explanation).
    poly._xCoeffs(1, 0) += 1;
    poly._yCoeffs(0, 1) += 1;
    return compose(geom::AffineTransform(geom::Extent2D(other.getPixelOrigin())),
                   compose(poly, other._cdInverse));
}

PolynomialTransform::PolynomialTransform(int order) : _xCoeffs(), _yCoeffs(), _u(), _v() {
    if (order < 0) {
        throw LSST_EXCEPT(pex::exceptions::LengthError, "PolynomialTransform order must be >= 0");
    }
    // Delay allocation until after error checking.
    _xCoeffs.reset(ndarray::Array<double, 2, 2>(ndarray::allocate(order + 1, order + 1)));
    _yCoeffs.reset(ndarray::Array<double, 2, 2>(ndarray::allocate(order + 1, order + 1)));
    _xCoeffs.setZero();
    _yCoeffs.setZero();
    _u = Eigen::VectorXd(order + 1);
    _v = Eigen::VectorXd(order + 1);
}

PolynomialTransform::PolynomialTransform(ndarray::Array<double const, 2, 0> const& xCoeffs,
                                         ndarray::Array<double const, 2, 0> const& yCoeffs)
        : _xCoeffs(ndarray::copy(xCoeffs)),
          _yCoeffs(ndarray::copy(yCoeffs)),
          _u(_xCoeffs.rows()),
          _v(_xCoeffs.rows()) {
    if (xCoeffs.getShape() != yCoeffs.getShape()) {
        throw LSST_EXCEPT(
                pex::exceptions::LengthError,
                (boost::format("X and Y coefficient matrices must have the same shape: "
                               " (%d,%d) != (%d,%d)") %
                 xCoeffs.getSize<0>() % xCoeffs.getSize<1>() % yCoeffs.getSize<0>() % yCoeffs.getSize<1>())
                        .str());
    }
    if (_xCoeffs.cols() != _xCoeffs.rows()) {
        throw LSST_EXCEPT(pex::exceptions::LengthError,
                          (boost::format("Coefficient matrices must be triangular, not trapezoidal: "
                                         " %d != %d ") %
                           _xCoeffs.rows() % _xCoeffs.cols())
                                  .str());
    }
}

PolynomialTransform::PolynomialTransform(PolynomialTransform const& other)
        : _xCoeffs(ndarray::copy(other.getXCoeffs())),
          _yCoeffs(ndarray::copy(other.getYCoeffs())),
          _u(other._u.size()),
          _v(other._v.size()) {}

PolynomialTransform::PolynomialTransform(PolynomialTransform&& other) : _xCoeffs(), _yCoeffs(), _u(), _v() {
    this->swap(other);
}

PolynomialTransform& PolynomialTransform::operator=(PolynomialTransform const& other) {
    if (&other != this) {
        PolynomialTransform tmp(other);
        tmp.swap(*this);
    }
    return *this;
}

PolynomialTransform& PolynomialTransform::operator=(PolynomialTransform&& other) {
    if (&other != this) {
        other.swap(*this);
    }
    return *this;
}

void PolynomialTransform::swap(PolynomialTransform& other) {
    _xCoeffs.swap(other._xCoeffs);
    _yCoeffs.swap(other._yCoeffs);
    _u.swap(other._u);
    _v.swap(other._v);
}

geom::AffineTransform PolynomialTransform::linearize(geom::Point2D const& in) const {
    double xu = 0.0, xv = 0.0, yu = 0.0, yv = 0.0, x = 0.0, y = 0.0;
    int const order = getOrder();
    detail::computePowers(_u, in.getX());
    detail::computePowers(_v, in.getY());
    for (int p = 0; p <= order; ++p) {
        for (int q = 0; q <= order; ++q) {
            if (p > 0) {
                xu += _xCoeffs(p, q) * p * _u[p - 1] * _v[q];
                yu += _yCoeffs(p, q) * p * _u[p - 1] * _v[q];
            }
            if (q > 0) {
                xv += _xCoeffs(p, q) * q * _u[p] * _v[q - 1];
                yv += _yCoeffs(p, q) * q * _u[p] * _v[q - 1];
            }
            x += _xCoeffs(p, q) * _u[p] * _v[q];
            y += _yCoeffs(p, q) * _u[p] * _v[q];
        }
    }
    geom::LinearTransform linear;
    linear.getMatrix()(0, 0) = xu;
    linear.getMatrix()(0, 1) = xv;
    linear.getMatrix()(1, 0) = yu;
    linear.getMatrix()(1, 1) = yv;
    geom::Point2D origin(x, y);
    return geom::AffineTransform(linear, origin - linear(in));
}

geom::Point2D PolynomialTransform::operator()(geom::Point2D const& in) const {
    int const order = getOrder();
    detail::computePowers(_u, in.getX());
    detail::computePowers(_v, in.getY());
    double x = 0;
    double y = 0;
    for (int p = 0; p <= order; ++p) {
        for (int q = 0; q <= order; ++q) {
            x += _xCoeffs(p, q) * _u[p] * _v[q];
            y += _yCoeffs(p, q) * _u[p] * _v[q];
        }
    }
    return geom::Point2D(x, y);
}

ScaledPolynomialTransform ScaledPolynomialTransform::convert(PolynomialTransform const& poly) {
    return ScaledPolynomialTransform(poly, geom::AffineTransform(), geom::AffineTransform());
}

ScaledPolynomialTransform ScaledPolynomialTransform::convert(SipForwardTransform const& sipForward) {
    ScaledPolynomialTransform result(sipForward.getPoly(),
                                     geom::AffineTransform(geom::Point2D(0, 0) - sipForward.getPixelOrigin()),
                                     geom::AffineTransform(sipForward.getCdMatrix()));
    // Account for the terms outside the sum in the SIP definition (see comment
    // earlier in the file for more explanation).
    result._poly._xCoeffs(1, 0) += 1;
    result._poly._yCoeffs(0, 1) += 1;
    return result;
}

ScaledPolynomialTransform ScaledPolynomialTransform::convert(SipReverseTransform const& sipReverse) {
    ScaledPolynomialTransform result(sipReverse.getPoly(), geom::AffineTransform(sipReverse._cdInverse),
                                     geom::AffineTransform(geom::Extent2D(sipReverse.getPixelOrigin())));
    result._poly._xCoeffs(1, 0) += 1;
    result._poly._yCoeffs(0, 1) += 1;
    return result;
}

ScaledPolynomialTransform::ScaledPolynomialTransform(PolynomialTransform const& poly,
                                                     geom::AffineTransform const& inputScaling,
                                                     geom::AffineTransform const& outputScalingInverse)
        : _poly(poly), _inputScaling(inputScaling), _outputScalingInverse(outputScalingInverse) {}

void ScaledPolynomialTransform::swap(ScaledPolynomialTransform& other) {
    _poly.swap(other._poly);
    std::swap(_inputScaling, other._inputScaling);
    std::swap(_outputScalingInverse, other._outputScalingInverse);
}

geom::AffineTransform ScaledPolynomialTransform::linearize(geom::Point2D const& in) const {
    return _outputScalingInverse * _poly.linearize(_inputScaling(in)) * _inputScaling;
}

geom::Point2D ScaledPolynomialTransform::operator()(geom::Point2D const& in) const {
    return _outputScalingInverse(_poly(_inputScaling(in)));
}

PolynomialTransform compose(geom::AffineTransform const& t1, PolynomialTransform const& t2) {
    typedef geom::AffineTransform AT;
    PolynomialTransform result(t2.getOrder());
    result._xCoeffs = t2._xCoeffs * t1[AT::XX] + t2._yCoeffs * t1[AT::XY];
    result._yCoeffs = t2._xCoeffs * t1[AT::YX] + t2._yCoeffs * t1[AT::YY];
    result._xCoeffs(0, 0) += t1[AT::X];
    result._yCoeffs(0, 0) += t1[AT::Y];
    return result;
}

PolynomialTransform compose(PolynomialTransform const& t1, geom::AffineTransform const& t2) {
    typedef geom::AffineTransform AT;
    int const order = t1.getOrder();
    if (order < 1) {
        PolynomialTransform t1a(1);
        t1a._xCoeffs(0, 0) = t1._xCoeffs(0, 0);
        t1a._yCoeffs(0, 0) = t1._yCoeffs(0, 0);
        return compose(t1a, t2);
    }
    detail::BinomialMatrix binomial(order);
    // For each of these, (e.g.) a[n] == pow(a, n)
    auto const t2u = detail::computePowers(t2[AT::X], order);
    auto const t2v = detail::computePowers(t2[AT::Y], order);
    auto const t2uu = detail::computePowers(t2[AT::XX], order);
    auto const t2uv = detail::computePowers(t2[AT::XY], order);
    auto const t2vu = detail::computePowers(t2[AT::YX], order);
    auto const t2vv = detail::computePowers(t2[AT::YY], order);
    PolynomialTransform result(order);
    for (int p = 0; p <= order; ++p) {
        for (int m = 0; m <= p; ++m) {
            for (int j = 0; j <= m; ++j) {
                for (int q = 0; p + q <= order; ++q) {
                    for (int n = 0; n <= q; ++n) {
                        for (int k = 0; k <= n; ++k) {
                            double z = binomial(p, m) * t2u[p - m] * binomial(m, j) * t2uu[j] * t2uv[m - j] *
                                       binomial(q, n) * t2v[q - n] * binomial(n, k) * t2vu[k] * t2vv[n - k];
                            result._xCoeffs(j + k, m + n - j - k) += t1._xCoeffs(p, q) * z;
                            result._yCoeffs(j + k, m + n - j - k) += t1._yCoeffs(p, q) * z;
                        }  // k
                    }      // n
                }          // q
            }              // j
        }                  // m
    }                      // p
    return result;
}

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
