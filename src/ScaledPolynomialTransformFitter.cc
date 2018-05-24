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

#include <map>

#include "Eigen/LU"  // for determinant, even though it's a 2x2 that doesn't use actual LU implementation

#include "boost/math/tools/minima.hpp"

#include "lsst/geom/Box.h"
#include "lsst/geom/AffineTransform.h"
#include "lsst/meas/astrom/ScaledPolynomialTransformFitter.h"
#include "lsst/meas/astrom/detail/polynomialUtils.h"
#include "lsst/afw/table/aggregates.h"
#include "lsst/afw/table/Match.h"
#include "lsst/afw/math/LeastSquares.h"

// When developing this code, it was used to add an in-line check for the
// correctness of a particularly complex calculation that was difficult to
// factor out into something unit-testable (the larger context that code is in
// *is* unit tested, which should be guard against regressions).  The test
// code is still in place, guarded by a check against this preprocessor macro;
// set it to 1 to re-enable that test during development, but be sure to
// set it back to zero before committing.
#define LSST_ScaledPolynomialTransformFitter_TEST_IN_PLACE 0

namespace lsst {
namespace meas {
namespace astrom {

// A singleton struct that manages the schema and keys for polynomial fitting catalogs.
class ScaledPolynomialTransformFitter::Keys {
public:
    afw::table::Schema schema;
    afw::table::Key<afw::table::RecordId> refId;
    afw::table::Key<afw::table::RecordId> srcId;
    afw::table::Point2DKey output;
    afw::table::Point2DKey input;
    afw::table::Point2DKey initial;
    afw::table::Point2DKey model;
    afw::table::CovarianceMatrixKey<float, 2> outputErr;
    // We use uint16 instead of Flag since it's the only bool we have here, we
    // may want NumPy views, and afw::table doesn't support [u]int8 fields.
    afw::table::Key<std::uint16_t> rejected;

    Keys(Keys const &) = delete;
    Keys(Keys &&) = delete;
    Keys &operator=(Keys const &) = delete;
    Keys &operator=(Keys &&) = delete;

    static Keys const &forMatches() {
        static Keys const it(0);
        return it;
    }

    static Keys const &forGrid() {
        static Keys const it;
        return it;
    }

private:
    Keys(int)
            : schema(),
              refId(schema.addField<afw::table::RecordId>("ref_id", "ID of reference object in this match.")),
              srcId(schema.addField<afw::table::RecordId>("src_id", "ID of source object in this match.")),
              output(afw::table::Point2DKey::addFields(schema, "src",
                                                       "source positions in pixel coordinates.", "pix")),
              input(afw::table::Point2DKey::addFields(schema, "intermediate",
                                                      "reference positions in intermediate world coordinates",
                                                      "deg")),
              initial(afw::table::Point2DKey::addFields(
                      schema, "initial", "reference positions transformed by initial WCS", "pix")),
              model(afw::table::Point2DKey::addFields(
                      schema, "model", "result of applying transform to reference positions", "pix")),
              outputErr(
                      afw::table::CovarianceMatrixKey<float, 2>::addFields(schema, "src", {"x", "y"}, "pix")),
              rejected(schema.addField<std::uint16_t>("rejected",
                                                      "True if the match should be rejected from the fit.")) {
        schema.getCitizen().markPersistent();
    }

    Keys()
            : schema(),
              output(afw::table::Point2DKey::addFields(
                      schema, "output", "grid output positions in intermediate world coordinates", "deg")),
              input(afw::table::Point2DKey::addFields(schema, "input",
                                                      "grid input positions in pixel coordinates.", "pix")),
              model(afw::table::Point2DKey::addFields(
                      schema, "model", "result of applying transform to input positions", "deg")) {
        schema.getCitizen().markPersistent();
    }
};

namespace {

// Return the AffineTransforms that maps the given (x,y) coordinates to lie within (-1, 1)x(-1, 1)
geom::AffineTransform computeScaling(afw::table::BaseCatalog const &data, afw::table::Point2DKey const &key) {
    geom::Box2D bbox;
    for (auto const &record : data) {
        bbox.include(geom::Point2D(record.get(key)));
    };
    return geom::AffineTransform(
                   geom::LinearTransform::makeScaling(0.5 * bbox.getWidth(), 0.5 * bbox.getHeight()))
                   .invert() *
           geom::AffineTransform(-geom::Extent2D(bbox.getCenter()));
}

}  // namespace

ScaledPolynomialTransformFitter ScaledPolynomialTransformFitter::fromMatches(
        int maxOrder, afw::table::ReferenceMatchVector const &matches, afw::geom::SkyWcs const &initialWcs,
        double intrinsicScatter) {
    Keys const &keys = Keys::forMatches();
    afw::table::BaseCatalog catalog(keys.schema);
    catalog.reserve(matches.size());
    float var2 = intrinsicScatter * intrinsicScatter;
    auto initialIwcToSky = getIntermediateWorldCoordsToSky(initialWcs);
    for (auto const &match : matches) {
        auto record = catalog.addNew();
        record->set(keys.refId, match.first->getId());
        record->set(keys.srcId, match.second->getId());
        record->set(keys.input, initialIwcToSky->applyInverse(match.first->getCoord()));
        record->set(keys.initial, initialWcs.skyToPixel(match.first->getCoord()));
        record->set(keys.output, match.second->getCentroid());
        record->set(keys.outputErr, match.second->getCentroidErr() + var2 * Eigen::Matrix2f::Identity());
        record->set(keys.rejected, false);
    }
    return ScaledPolynomialTransformFitter(catalog, keys, maxOrder, intrinsicScatter,
                                           computeScaling(catalog, keys.input),
                                           computeScaling(catalog, keys.output));
}

ScaledPolynomialTransformFitter ScaledPolynomialTransformFitter::fromGrid(
        int maxOrder, geom::Box2D const &bbox, int nGridX, int nGridY,
        ScaledPolynomialTransform const &toInvert) {
    Keys const &keys = Keys::forGrid();
    afw::table::BaseCatalog catalog(keys.schema);
    catalog.reserve(nGridX * nGridY);
    geom::Extent2D dx(bbox.getWidth() / nGridX, 0.0);
    geom::Extent2D dy(0.0, bbox.getHeight() / nGridY);
    for (int iy = 0; iy < nGridY; ++iy) {
        for (int ix = 0; ix < nGridX; ++ix) {
            geom::Point2D point = bbox.getMin() + dx * ix + dy * iy;
            auto record = catalog.addNew();
            record->set(keys.output, point);
            record->set(keys.input, toInvert(point));
        }
    }
    return ScaledPolynomialTransformFitter(catalog, keys, maxOrder, 0.0, computeScaling(catalog, keys.input),
                                           computeScaling(catalog, keys.output));
}

ScaledPolynomialTransformFitter::ScaledPolynomialTransformFitter(afw::table::BaseCatalog const &data,
                                                                 Keys const &keys, int maxOrder,
                                                                 double intrinsicScatter,
                                                                 geom::AffineTransform const &inputScaling,
                                                                 geom::AffineTransform const &outputScaling)
        : _keys(keys),
          _intrinsicScatter(intrinsicScatter),
          _data(data),
          _outputScaling(outputScaling),
          _transform(PolynomialTransform(maxOrder), inputScaling, outputScaling.invert()),
          _vandermonde(data.size(), detail::computePackedSize(maxOrder)) {
    // Create a matrix that evaluates the max-order polynomials of all the (scaled) input positions;
    // we'll extract subsets of this later when fitting to a subset of the matches and a lower order.
    for (std::size_t i = 0; i < data.size(); ++i) {
        geom::Point2D input = getInputScaling()(_data[i].get(_keys.input));
        // x[k] == pow(x, k), y[k] == pow(y, k)
        detail::computePowers(_transform._poly._u, input.getX());
        detail::computePowers(_transform._poly._v, input.getY());
        // We pack coefficients in the following order:
        // (0,0), (0,1), (1,0), (0,2), (1,1), (2,0)
        // Note that this lets us choose the just first N(N+1)/2 columns to
        // evaluate an Nth order polynomial, even if N < maxOrder.
        for (int n = 0, j = 0; n <= maxOrder; ++n) {
            for (int p = 0, q = n; p <= n; ++p, --q, ++j) {
                _vandermonde(i, j) = _transform._poly._u[p] * _transform._poly._v[q];
            }
        }
    }
}

void ScaledPolynomialTransformFitter::fit(int order) {
    int maxOrder = _transform.getPoly().getOrder();
    if (order < 0) {
        order = maxOrder;
    }
    if (order > maxOrder) {
        throw LSST_EXCEPT(
                pex::exceptions::LengthError,
                (boost::format("Order (%d) exceeded maximum order for the fitter (%d)") % order % maxOrder)
                        .str());
    }

    int const packedSize = detail::computePackedSize(order);
    std::size_t nGood = 0;
    if (_keys.rejected.isValid()) {
        for (auto const &record : _data) {
            if (!record.get(_keys.rejected)) {
                ++nGood;
            }
        }
    } else {
        nGood = _data.size();
    }
    // One block of the block-diagonal (2x2) unweighted design matrix M;
    // m[i,j] = u_i^{p(j)} v_i^{q(j)}.  The two nonzero blocks are the same,
    // because we're using the same polynomial basis for x and y.
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nGood, packedSize);
    // vx, vy: (2x1) blocks of the unweighted data vector v
    Eigen::VectorXd vx = Eigen::VectorXd::Zero(nGood);
    Eigen::VectorXd vy = Eigen::VectorXd::Zero(nGood);
    // sxx, syy, sxy: (2x2) blocks of the covariance matrix S, each of which is
    // individually diagonal.
    Eigen::ArrayXd sxx(nGood);
    Eigen::ArrayXd syy(nGood);
    Eigen::ArrayXd sxy(nGood);
    Eigen::Matrix2d outS = _outputScaling.getLinear().getMatrix();
    for (std::size_t i1 = 0, i2 = 0; i1 < _data.size(); ++i1) {
        // check that the 'rejected' field (== 'not outlier-rejected') is both
        // present in the schema and not rejected.
        if (!_keys.rejected.isValid() || !_data[i1].get(_keys.rejected)) {
            geom::Point2D output = _outputScaling(_data[i1].get(_keys.output));
            vx[i2] = output.getX();
            vy[i2] = output.getY();
            m.row(i2) = _vandermonde.row(i1).head(packedSize);
            if (_keys.outputErr.isValid()) {
                Eigen::Matrix2d modelErr =
                        outS * _data[i1].get(_keys.outputErr).cast<double>() * outS.adjoint();
                sxx[i2] = modelErr(0, 0);
                sxy[i2] = modelErr(0, 1);
                syy[i2] = modelErr(1, 1);
            } else {
                sxx[i2] = 1.0;
                sxy[i2] = 0.0;
                syy[i2] = 1.0;
            }
            ++i2;
        }
    }
    // Do a blockwise inverse of S.  Note that the result F is still symmetric
    Eigen::ArrayXd fxx = 1.0 / (sxx - sxy.square() / syy);
    Eigen::ArrayXd fyy = 1.0 / (syy - sxy.square() / sxx);
    Eigen::ArrayXd fxy = -(sxy / sxx) * fyy;
#ifdef LSST_ScaledPolynomialTransformFitter_TEST_IN_PLACE
    assert((sxx * fxx + sxy * fxy).isApproxToConstant(1.0));
    assert((syy * fyy + sxy * fxy).isApproxToConstant(1.0));
    assert((sxx * fxy).isApprox(-sxy * fyy));
    assert((sxy * fxx).isApprox(-syy * fxy));
#endif
    // Now that we've got all the block quantities, we'll form the full normal equations matrix.
    // That's H = M^T F M:
    Eigen::MatrixXd h(2 * packedSize, 2 * packedSize);
    h.topLeftCorner(packedSize, packedSize) = m.adjoint() * fxx.matrix().asDiagonal() * m;
    h.topRightCorner(packedSize, packedSize) = m.adjoint() * fxy.matrix().asDiagonal() * m;
    h.bottomLeftCorner(packedSize, packedSize) = h.topRightCorner(packedSize, packedSize).adjoint();
    h.bottomRightCorner(packedSize, packedSize) = m.adjoint() * fyy.matrix().asDiagonal() * m;
    // And here's the corresponding RHS vector, g = M^T F v
    Eigen::VectorXd g(2 * packedSize);
    g.head(packedSize) = m.adjoint() * (fxx.matrix().asDiagonal() * vx + fxy.matrix().asDiagonal() * vy);
    g.tail(packedSize) = m.adjoint() * (fxy.matrix().asDiagonal() * vx + fyy.matrix().asDiagonal() * vy);
    // Solve the normal equations.
    auto lstsq = afw::math::LeastSquares::fromNormalEquations(h, g);
    auto solution = lstsq.getSolution();
    // Unpack the solution vector back into the polynomial coefficient matrices.
    for (int n = 0, j = 0; n <= order; ++n) {
        for (int p = 0, q = n; p <= n; ++p, --q, ++j) {
            _transform._poly._xCoeffs(p, q) = solution[j];
            _transform._poly._yCoeffs(p, q) = solution[j + packedSize];
        }
    }
}

void ScaledPolynomialTransformFitter::updateModel() {
    for (auto &record : _data) {
        record.set(_keys.model, _transform(record.get(_keys.input)));
    }
}

double ScaledPolynomialTransformFitter::updateIntrinsicScatter() {
    if (!_keys.rejected.isValid()) {
        throw LSST_EXCEPT(pex::exceptions::LogicError,
                          "Cannot compute intrinsic scatter on fitter initialized with fromGrid.");
    }
    double newIntrinsicScatter = computeIntrinsicScatter();
    float varDiff = newIntrinsicScatter * newIntrinsicScatter - _intrinsicScatter * _intrinsicScatter;
    for (auto &record : _data) {
        record.set(_keys.outputErr, record.get(_keys.outputErr) + varDiff * Eigen::Matrix2f::Identity());
    }
    _intrinsicScatter = newIntrinsicScatter;
    return _intrinsicScatter;
}

double ScaledPolynomialTransformFitter::computeIntrinsicScatter() const {
    // We model the variance matrix of each match as the sum of the intrinsic variance (which we're
    // trying to fit) and the per-source measurement uncertainty.
    // We start by computing the variance directly, which yields the sum.
    // At the same time, we find the maximum per-source variance.  Since the per-source uncertainties
    // are actually 2x2 matrices, we use the square of the major axis of that ellipse.
    double directVariance = 0.0;          // direct estimate of total scatter (includes measurement errors)
    double maxMeasurementVariance = 0.0;  // maximum of the per-match measurement uncertainties
    double oldIntrinsicVariance = _intrinsicScatter * _intrinsicScatter;
    std::size_t nGood = 0;
    for (auto const &record : _data) {
        if (!_keys.rejected.isValid() || !record.get(_keys.rejected)) {
            auto delta = record.get(_keys.output) - record.get(_keys.model);
            directVariance += 0.5 * delta.computeSquaredNorm();
            double cxx = _keys.outputErr.getElement(record, 0, 0) - oldIntrinsicVariance;
            double cyy = _keys.outputErr.getElement(record, 1, 1) - oldIntrinsicVariance;
            double cxy = _keys.outputErr.getElement(record, 0, 1);
            // square of semimajor axis of uncertainty error ellipse
            double ca2 = 0.5 * (cxx + cyy + std::sqrt(cxx * cxx + cyy * cyy + 4 * cxy * cxy - 2 * cxx * cyy));
            maxMeasurementVariance = std::max(maxMeasurementVariance, ca2);
            ++nGood;
        }
    }
    directVariance /= nGood;

    // Function that computes the -log likelihood of the current deltas with
    // the variance modeled as described above.
    auto logLikelihood = [this](double intrinsicVariance) {
        // Uncertainties in the table right now include the old intrinsic scatter; need to
        // subtract it off as we add the new one in.
        double varDiff = intrinsicVariance - this->_intrinsicScatter * this->_intrinsicScatter;
        double q = 0.0;
        for (auto &record : this->_data) {
            double x = record.get(this->_keys.output.getX()) - record.get(_keys.model.getX());
            double y = record.get(this->_keys.output.getY()) - record.get(_keys.model.getY());
            double cxx = this->_keys.outputErr.getElement(record, 0, 0) + varDiff;
            double cyy = this->_keys.outputErr.getElement(record, 1, 1) + varDiff;
            double cxy = this->_keys.outputErr.getElement(record, 0, 1);
            double det = cxx * cyy - cxy * cxy;
            q += (x * x * cyy - 2 * x * y * cxy + y * y * cxx) / det + std::log(det);
        }
        return q;
    };

    // directVariance brackets the intrinsic variance from above, and this quantity
    // brackets it from below:
    double minIntrinsicVariance = std::max(0.0, directVariance - maxMeasurementVariance);

    // Minimize the negative log likelihood to find the best-fit intrinsic variance.
    static constexpr int BITS_REQUIRED = 16;  // solution good to ~1E-4
    boost::uintmax_t maxIterations = 20;
    auto result = boost::math::tools::brent_find_minima(logLikelihood, minIntrinsicVariance, directVariance,
                                                        BITS_REQUIRED, maxIterations);
    return std::sqrt(result.first);  // return RMS instead of variance
}

std::pair<double, std::size_t> ScaledPolynomialTransformFitter::rejectOutliers(
        OutlierRejectionControl const &ctrl) {
    // If the 'rejected' field isn't present in the schema (because the fitter
    // was constructed with fromGrid), we can't do outlier rejection.
    if (!_keys.rejected.isValid()) {
        throw LSST_EXCEPT(pex::exceptions::LogicError,
                          "Cannot reject outliers on fitter initialized with fromGrid.");
    }
    if (static_cast<std::size_t>(ctrl.nClipMin) >= _data.size()) {
        throw LSST_EXCEPT(
                pex::exceptions::LogicError,
                (boost::format("Not enough values (%d) to clip %d.") % _data.size() % ctrl.nClipMin).str());
    }
    std::map<double, afw::table::BaseRecord *> rankings;
    for (auto &record : _data) {
        Eigen::Matrix2d cov = record.get(_keys.outputErr).cast<double>();
        Eigen::Vector2d d = (record.get(_keys.output) - record.get(_keys.model)).asEigen();
        double r2 = d.dot(cov.inverse() * d);
        rankings.insert(std::make_pair(r2, &record));
    }
    auto cutoff = rankings.upper_bound(ctrl.nSigma * ctrl.nSigma);
    int nClip = 0, nGood = 0;
    for (auto iter = rankings.begin(); iter != cutoff; ++iter) {
        iter->second->set(_keys.rejected, false);
        ++nGood;
    }
    for (auto iter = cutoff; iter != rankings.end(); ++iter) {
        iter->second->set(_keys.rejected, true);
        ++nClip;
    }
    assert(static_cast<std::size_t>(nGood + nClip) == _data.size());
    while (nClip < ctrl.nClipMin) {
        --cutoff;
        cutoff->second->set(_keys.rejected, true);
        ++nClip;
    }
    while (nClip > ctrl.nClipMax && cutoff != rankings.end()) {
        cutoff->second->set(_keys.rejected, false);
        ++cutoff;
        --nClip;
    }
    std::pair<double, std::size_t> result(ctrl.nSigma, nClip);
    if (cutoff != rankings.end()) {
        result.first = std::sqrt(cutoff->first);
    }
    return result;
}

}  // namespace astrom
}  // namespace meas
}  // namespace lsst
