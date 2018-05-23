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
#ifndef LSST_MEAS_ASTROM_PolynomialTransform_INCLUDED
#define LSST_MEAS_ASTROM_PolynomialTransform_INCLUDED

#include "ndarray/eigen.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/AffineTransform.h"

namespace lsst { namespace meas { namespace astrom {

class SipForwardTransform;
class SipReverseTransform;
class ScaledPolynomialTransform;

/**
 *  A 2-d coordinate transform represented by a pair of standard polynomials
 *  (one for each coordinate).
 *
 *  PolynomialTransform instances should be confined to a single thread.
 */
class PolynomialTransform {
public:

    /**
     *  Convert a ScaledPolynomialTransform to an equivalent PolynomialTransform.
     */
    static PolynomialTransform convert(ScaledPolynomialTransform const & other);

    /**
     *  Convert a SipForwardTransform to an equivalent PolynomialTransform.
     */
    static PolynomialTransform convert(SipForwardTransform const & other);

    /**
     *  Convert a SipReverseTransform to an equivalent PolynomialTransform.
     */
    static PolynomialTransform convert(SipReverseTransform const & other);

    /**
     *  Construct a new transform from existing coefficient arrays.
     *
     *  For both input arguments, the array element at [p, q] corresponds
     *  to the polynomial term x^p y^q.
     *
     *  Both arrays are expected be square and triangular; if N is the
     *  order of the transform, both arrays should be (N+1)x(N+1), and
     *  elements with p + q > N should be zero.
     */
    PolynomialTransform(
        ndarray::Array<double const,2,0> const & xCoeffs,
        ndarray::Array<double const,2,0> const & yCoeffs
    );

    /**
     *  Copy constructor.
     *
     *  Coefficient arrays are deep-copied.
     */
    PolynomialTransform(PolynomialTransform const & other);

    /**
     *  Move constructor.
     *
     *  Coefficient arrays are moved.
     */
    PolynomialTransform(PolynomialTransform && other);

    /**
     *  Copy assignment.
     *
     *  Coefficient arrays are deep-copied.
     */
    PolynomialTransform & operator=(PolynomialTransform const & other);

    /**
     *  Move constructor.
     *
     *  Coefficient arrays are moved.
     */
    PolynomialTransform & operator=(PolynomialTransform && other);

    /// Lightweight swap.
    void swap(PolynomialTransform & other);

    /// Return the order of the polynomials.
    int getOrder() const { return _xCoeffs.rows() - 1; }

    /**
     * 2-D polynomial coefficients that compute the output x coordinate.
     *
     * Indexing the result by [p][q] gives the coefficient of
     * @f$x_{\mathrm{in}}^p\,y_{\mathrm{in}}^q@f$.
     */
    ndarray::Array<double const,2,2> getXCoeffs() const { return _xCoeffs.shallow(); }

    /**
     * 2-D polynomial coefficients that compute the output x coordinate.
     *
     * Indexing the result by [p][q] gives the coefficient of
     * @f$x_{\mathrm{in}}^p\,y_{\mathrm{in}}^q@f$.
     */
    ndarray::Array<double const,2,2> getYCoeffs() const { return _yCoeffs.shallow(); }

    /**
     * Return an approximate affine transform at the given point.
     */
    geom::AffineTransform linearize(geom::Point2D const & in) const;

    /**
     * Apply the transform to a point.
     */
    geom::Point2D operator()(geom::Point2D const & in) const;

private:

    PolynomialTransform(int order);

    friend PolynomialTransform compose(geom::AffineTransform const & t1, PolynomialTransform const & t2);
    friend PolynomialTransform compose(PolynomialTransform const & t1, geom::AffineTransform const & t2);
    friend class ScaledPolynomialTransformFitter;
    friend class SipForwardTransform;
    friend class SipReverseTransform;
    friend class ScaledPolynomialTransform;

    ndarray::EigenView<double,2,2> _xCoeffs;
    ndarray::EigenView<double,2,2> _yCoeffs;
    mutable Eigen::VectorXd _u; // workspace for operator() and linearize
    mutable Eigen::VectorXd _v;
};

/**
 *  A 2-d coordinate transform represented by a lazy composition of an AffineTransform,
 *  a PolynomialTransform, and another AffineTransform.
 *
 *  ScaledPolynomialTransform instances should be confined to a single thread.
 */
class ScaledPolynomialTransform {
public:

    /**
     *  Convert a PolynomialTransform to an equivalent ScaledPolynomialTransform.
     *
     *  This simply inserts identity AffineTransforms before and after applying
     *  the given PolynomialTransform.
     */
    static ScaledPolynomialTransform convert(PolynomialTransform const & poly);

    /**
     *  Convert a SipForwardTransform to an equivalent ScaledPolynomialTransform.
     *
     *  The input transform's CRPIX offset and CD matrix scaling are used to define
     *  the input and output affine transforms, respectively, leaving the polynomial
     *  coefficients unmodified.
     */
    static ScaledPolynomialTransform convert(SipForwardTransform const & sipForward);

    /**
     *  Convert a SipForwardTransform to an equivalent ScaledPolynomialTransform.
     *
     *  The input transform's CD matrix scaling and CRPIX offset are used to define
     *  the input and output affine transforms, respectively, leaving the polynomial
     *  coefficients unmodified.
     */
    static ScaledPolynomialTransform convert(SipReverseTransform const & sipReverse);

    /**
     *  Construct a new ScaledPolynomialTransform from its constituents.
     *
     *  @param[in]  poly           A PolynomialTransform to be applied to points.
     *                             after the inputScaling transform.
     *  @param[in]  inputScaling   An AffineTransform to be applied immediately
     *                             to input points.
     *  @param[in]  outputScalingInverse  An AffineTransform to be applied to points
     *                                    after the PolynomialTransform.
     */
    ScaledPolynomialTransform(
        PolynomialTransform const & poly,
        geom::AffineTransform const & inputScaling,
        geom::AffineTransform const & outputScalingInverse
    );

    ScaledPolynomialTransform(ScaledPolynomialTransform const & other) = default;

    ScaledPolynomialTransform(ScaledPolynomialTransform && other) = default;

    ScaledPolynomialTransform & operator=(ScaledPolynomialTransform const & other) = default;

    ScaledPolynomialTransform & operator=(ScaledPolynomialTransform && other) = default;

    void swap(ScaledPolynomialTransform & other);

    /// Return the polynomial transform applied after the input scaling.
    PolynomialTransform const & getPoly() const { return _poly; }

    /// Return the first affine transform applied to input points.
    geom::AffineTransform const & getInputScaling() const { return _inputScaling; }

    /// Return the affine transform applied to points after the polynomial transform.
    geom::AffineTransform const & getOutputScalingInverse() const { return _outputScalingInverse; }

    /**
     * Return an approximate affine transform at the given point.
     */
    geom::AffineTransform linearize(geom::Point2D const & in) const;

    /**
     * Apply the transform to a point.
     */
    geom::Point2D operator()(geom::Point2D const & in) const;

private:
    friend class ScaledPolynomialTransformFitter;
    PolynomialTransform _poly;
    geom::AffineTransform _inputScaling;
    geom::AffineTransform _outputScalingInverse;
};

/**
 *  Return a PolynomialTransform that is equivalent to the composition t1(t2())
 *
 *  The returned composition would be exact in ideal arithmetic, but may suffer from
 *  significant round-off error for high-order polynomials.
 */
PolynomialTransform compose(geom::AffineTransform const & t1, PolynomialTransform const & t2);

/**
 *  Return a PolynomialTransform that is equivalent to the composition t1(t2())
 *
 *  The returned composition would be exact in ideal arithmetic, but may suffer from
 *  significant round-off error for high-order polynomials.
 */
PolynomialTransform compose(PolynomialTransform const & t1, geom::AffineTransform const & t2);

}}} // namespace lsst::meas::astrom

#endif // !LSST_MEAS_ASTROM_PolynomialTransform_INCLUDED