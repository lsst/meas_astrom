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
#ifndef LSST_MEAS_ASTROM_SipTransform_INCLUDED
#define LSST_MEAS_ASTROM_SipTransform_INCLUDED

#include "lsst/afw/geom/LinearTransform.h"
#include "lsst/meas/astrom/PolynomialTransform.h"


namespace lsst { namespace afw { namespace image {

class TanWcs;

} // namespace image

namespace coord {

class Coord;

} // namespace coord

} // namespace afw

namespace meas { namespace astrom {

/**
 *  Base class for SIP transform objects.
 *
 *  This class simply provides some getters for its derived classes.
 *  It should not be used directly, and does not define a polymorphic
 *  interface.
 */
class SipTransformBase {
public:

    /**
     *  Return the pixel origin (CRPIX) of the transform.
     */
    afw::geom::Point2D const & getPixelOrigin() const { return _pixelOrigin; }

    /**
     *  Return the CD matrix of the transform.
     */
    afw::geom::LinearTransform const & getCDMatrix() const { return _cdMatrix; }

    /**
     *  Return the polynomial component of the transform (A,B) or (AP,BP).
     */
    PolynomialTransform const & getPoly() const { return _poly; }

protected:

    /**
     *  Construct a SipTransformBase from its components.
     *
     *  See SipForwardTransform and SipReverseTransform for more extensive
     *  definitions.
     *
     *  @param[in]   pixelOrigin    CRPIX @f$(u_0,v_0)@f$ (zero-indexed)
     *  @param[in]   cdMatrix       CD matrix @f$Z@f$
     *  @param[in]   poly           Either the forward or reverse
     *                              SIP polynomial (depending on the derived
     *                              class).
     */
    SipTransformBase(
        afw::geom::Point2D const & pixelOrigin,
        afw::geom::LinearTransform const & cdMatrix,
        PolynomialTransform const & poly
    ) : _pixelOrigin(pixelOrigin),
        _cdMatrix(cdMatrix),
        _poly(poly)
    {}

    SipTransformBase(SipTransformBase const & other) = default;
    SipTransformBase(SipTransformBase && other) = default;
    SipTransformBase & operator=(SipTransformBase const & other) = default;
    SipTransformBase & operator=(SipTransformBase && other) = default;

    void swap(SipTransformBase & other) {
        std::swap(_pixelOrigin, other._pixelOrigin);
        std::swap(_cdMatrix, other._cdMatrix);
        _poly.swap(other._poly);
    }

    void transformPixelsInPlace(afw::geom::AffineTransform const & s);

    afw::geom::Point2D _pixelOrigin;
    afw::geom::LinearTransform _cdMatrix;
    PolynomialTransform _poly;
};


/**
 *  A transform that maps pixel coordinates to intermediate world coordinates
 *  according to the SIP convention.
 *
 *  The SIP forward transform is defined as
 *  @f[
 *     \left[\begin{array}{ c }
 *       x \\
 *       y
 *     \end{array}\right]
 *     = \mathbf{Z}
 *     \left[\begin{array}{ c }
 *       (u - u_0) + {\displaystyle\sum_{p,q}^{2 \le p + q \le N}} \mathrm{A}_{p,q} (u-u_0)^p (v-v_0)^q \\
 *       (v - v_0) + {\displaystyle\sum_{p,q}^{2 \le p + q \le N}} \mathrm{B}_{p,q} (u-u_0)^p (v-v_0)^q
 *     \end{array}\right]
 *  @f]
 *  where
 *   - @f$(u,v)@f$ are pixel coordinates (zero-indexed).
 *   - @f$(x,y)@f$ are "intermediate world coordinates" -- the result of
 *     applying the gnomonic (TAN) projection at sky origin CRVAL to sky
 *     coordinates).
 *   - @f$\mathbf{Z}@f$ is the @f$2 \times 2@f$ linear transform (@f$\mathrm{CD}@f$) matrix.
 *   - @f$(u_0,v_0)@f$ is the pixel origin @f$\mathrm{CRPIX}@f$ (but zero-indexed;
 *     the FITS standard is 1-indexed).
 *   - @f$\mathrm{A}@f$, @f$\mathrm{B}@f$ are the polynomial coefficients of
 *     the forward transform.
 *
 *  The SIP convention encourages (but does not require) nulling the zeroth-
 *  and first-order elements of
 *  @f$\mathrm{A}@f$ and @f$\mathrm{B}@f$, which ensures the representation of
 *  a given transform is unique.  This also makes fitting a SIP transform to
 *  data a nonlinear operation, as well as making the conversion from standard
 *  polynomial transforms to SIP form impossible in general.  Accordingly,
 *  this class does not attempt to null low-order polynomial terms at all
 *  when converting from other transforms.
 *
 *  SipForwardTransform instances should be confined to a single thread.
 */
class SipForwardTransform : public SipTransformBase {
public:

    /**
     *  Convert a PolynomialTransform to an equivalent SipForwardTransform.
     *
     *  @param[in]   poly           PolynomialTransform to convert.
     *  @param[in]   pixelOrigin    CRPIX @f$(u_0,v_0)@f$ (zero-indexed)
     *  @param[in]   cdMatrix       CD matrix @f$Z@f$
     */
    static SipForwardTransform convert(
        PolynomialTransform const & poly,
        afw::geom::Point2D const & pixelOrigin,
        afw::geom::LinearTransform const & cdMatrix
    );

    /**
     *  Convert a ScaledPolynomialTransform to an equivalent SipForwardTransform.
     *
     *  @param[in]   scaled         ScaledPolynomialTransform to convert.
     *  @param[in]   pixelOrigin    CRPIX @f$(u_0,v_0)@f$ (zero-indexed)
     *  @param[in]   cdMatrix       CD matrix @f$Z@f$
     */
    static SipForwardTransform convert(
        ScaledPolynomialTransform const & scaled,
        afw::geom::Point2D const & pixelOrigin,
        afw::geom::LinearTransform const & cdMatrix
    );

    /**
     *  Convert a ScaledPolynomialTransform to an equivalent SipForwardTransform.
     *
     *  The pixel origin CRPIX and CD matrix are defined to reproduce the translation
     *  and linear transformation in the ScaledPolynomialTransform's input and
     *  output scalings (respectively).
     */
    static SipForwardTransform convert(ScaledPolynomialTransform const & scaled);

    /**
     *  Construct a SipForwardTransform from its components.
     *
     *  @param[in]   pixelOrigin    CRPIX @f$(u_0,v_0)@f$ (zero-indexed).
     *  @param[in]   cdMatrix       CD matrix @f$Z@f$
     *  @param[in]   forwardSipPoly Polynomial transform @f$(A,B)@f$
     */
    SipForwardTransform(
        afw::geom::Point2D const & pixelOrigin,
        afw::geom::LinearTransform const & cdMatrix,
        PolynomialTransform const & forwardSipPoly
    ) :
        SipTransformBase(pixelOrigin, cdMatrix, forwardSipPoly)
    {}

    SipForwardTransform(SipForwardTransform const & other) = default;

    SipForwardTransform(SipForwardTransform && other) = default;

    SipForwardTransform & operator=(SipForwardTransform const & other) = default;

    SipForwardTransform & operator=(SipForwardTransform && other) = default;

    void swap(SipForwardTransform & other) {
        SipTransformBase::swap(other);
    }

    /**
     * Return an approximate affine transform at the given point.
     */
    afw::geom::AffineTransform linearize(afw::geom::Point2D const & in) const;

    /**
     * Apply the transform to a point.
     */
    afw::geom::Point2D operator()(afw::geom::Point2D const & uv) const;

    /**
     * Return a new forward SIP transform that includes a transformation of
     * the pixel coordinate system by the given affine transform.
     */
    SipForwardTransform transformPixels(afw::geom::AffineTransform const & s) const;

};


/**
 *  A transform that maps intermediate world coordinates to pixel coordinates
 *  according to the SIP convention.
 *
 *  The SIP reverse transform is defined as
 *  @f[
 *     \left[\begin{array}{ c }
 *       u \\
 *       v
 *     \end{array}\right]
 *     = \left[\begin{array}{ c }
 *       u_0 + U + {\displaystyle\sum_{p,q}^{0 \le p + q \le N}} \mathrm{AP}_{p,q} U^p V^q \\
 *       v_0 + V + {\displaystyle\sum_{p,q}^{0 \le p + q \le N}} \mathrm{BP}_{p,q} U^p V^q \\
 *     \end{array}\right]
 *  @f]
 *  with
 *  @f[
 *     \left[\begin{array}{ c }
 *       U \\
 *       V
 *     \end{array}\right]
 *     = \mathbf{Z}^{-1}
 *     \left[\begin{array}{ c }
 *       x \\
 *       y
 *     \end{array}\right]
 *  @f]
 *  and
 *   - @f$(u,v)@f$ are pixel coordinates.
 *   - @f$(x,y)@f$ are "intermediate world coordinates" -- the result of
 *     applying the gnomonic (TAN) projection at sky origin CRVAL to sky
 *     coordinates).
 *   - @f$\mathbf{Z}@f$ is the @f$2 \times 2@f$ linear transform (@f$\mathrm{CD}@f$) matrix.
 *   - @f$(u_0,v_0)@f$ is the pixel origin @f$\mathrm{CRPIX}@f$ (but zero-indexed;
 *     the FITS standard is 1-indexed).
 *   - @f$\mathrm{AP}@f$, @f$\mathrm{BP}@f$ are the polynomial coefficients of
 *     the reverse transform.
 *
 *  SipForwardTransform instances should be confined to a single thread.
 */
class SipReverseTransform : public SipTransformBase {
public:

    /**
     *  Convert a PolynomialTransform to an equivalent SipReverseTransform.
     *
     *  @param[in]   poly           PolynomialTransform to convert.
     *  @param[in]   pixelOrigin    CRPIX @f$(u_0,v_0)@f$ (zero-indexed)
     *  @param[in]   cdMatrix       CD matrix @f$Z@f$
     */
    static SipReverseTransform convert(
        PolynomialTransform const & poly,
        afw::geom::Point2D const & pixelOrigin,
        afw::geom::LinearTransform const & cdMatrix
    );

    /**
     *  Convert a ScaledPolynomialTransform to an equivalent SipReverseTransform.
     *
     *  @param[in]   scaled         ScaledPolynomialTransform to convert.
     *  @param[in]   pixelOrigin    CRPIX @f$(u_0,v_0)@f$ (zero-indexed)
     *  @param[in]   cdMatrix       CD matrix @f$Z@f$
     */
    static SipReverseTransform convert(
        ScaledPolynomialTransform const & scaled,
        afw::geom::Point2D const & pixelOrigin,
        afw::geom::LinearTransform const & cdMatrix
    );

    /**
     *  Convert a ScaledPolynomialTransform to an equivalent SipReverseTransform.
     *
     *  The pixel origin CRPIX and CD matrix are defined to reproduce the translation
     *  and linear transformation in the ScaledPolynomialTransforms output and
     *  input scalings (respectively).
     */
    static SipReverseTransform convert(ScaledPolynomialTransform const & scaled);

    /**
     *  Construct a SipReverseTransform from its components.
     *
     *  @param[in]   pixelOrigin    CRPIX @f$(u_0,v_0)@f$ (zero-indexed)
     *  @param[in]   cdMatrix       CD matrix @f$Z@f$
     *  @param[in]   reverseSipPoly Polynomial transform @f$(AP,BP)@f$
     */
    SipReverseTransform(
        afw::geom::Point2D const & pixelOrigin,
        afw::geom::LinearTransform const & cdMatrix,
        PolynomialTransform const & reverseSipPoly
    ) : SipTransformBase(pixelOrigin, cdMatrix, reverseSipPoly),
        _cdInverse(cdMatrix.invert())
    {}

    SipReverseTransform(SipReverseTransform const & other) = default;

    SipReverseTransform(SipReverseTransform && other) = default;

    SipReverseTransform & operator=(SipReverseTransform const & other) = default;

    SipReverseTransform & operator=(SipReverseTransform && other) = default;

    void swap(SipReverseTransform & other) {
        SipTransformBase::swap(other);
        std::swap(_cdInverse, other._cdInverse);
    }

    /**
     * Return an approximate affine transform at the given point.
     */
    afw::geom::AffineTransform linearize(afw::geom::Point2D const & in) const;

    /**
     * Apply the transform to a point.
     */
    afw::geom::Point2D operator()(afw::geom::Point2D const & xy) const;

    /**
     * Return a new reverse SIP transform that includes a transformation of
     * the pixel coordinate system by the given affine transform.
     */
    SipReverseTransform transformPixels(afw::geom::AffineTransform const & s) const;

private:
    friend class PolynomialTransform;
    friend class ScaledPolynomialTransform;
    afw::geom::LinearTransform _cdInverse;
};

/**
 *  Create a new TAN SIP Wcs from a pair of SIP transforms and the sky origin.
 *
 *  @param[in]   sipForward    Mapping from pixel coordinates to intermediate
 *                             world coordinates.
 *  @param[in]   sipReverse    Mapping from intermediate world coordinates to
 *                             pixel coordinates.
 *  @param[in]   skyOrigin     Position of the gnomonic projection that maps
 *                             sky coordinates to intermediate world coordinates
 *                             (CRVAL).
 *
 *  @throw pex::exceptions::InvalidParameterError if the forward and reverse
 *         SIP transforms have different CRPIX values or CD matrices.
 */
std::shared_ptr<afw::image::TanWcs> makeWcs(
    SipForwardTransform const & sipForward,
    SipReverseTransform const & sipReverse,
    afw::coord::Coord const & skyOrigin
);

}}} // namespace lsst::meas::astrom

#endif // !LSST_MEAS_ASTROM_SipTransform_INCLUDED