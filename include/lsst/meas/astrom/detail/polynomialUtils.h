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

#ifndef LSST_MEAS_ASTROM_DETAIL_polynomialUtils_h_INCLUDED
#define LSST_MEAS_ASTROM_DETAIL_polynomialUtils_h_INCLUDED

#include "Eigen/Core"

namespace lsst {
namespace meas {
namespace astrom {
namespace detail {

/**
 *  Compute the index of the first coefficient with the given order in
 *  a packed 2-d polynomial coefficient array.
 *
 *  This defines the ordering as
 *  @code
 *  [(0,0), (0,1), (1,0), (0,2), (1,1), (2,0), ...]
 *  @endcode
 *  (or the same with indices swapped).
 */
inline int computePackedOffset(int order) { return (order * (order + 1)) / 2; }

/**
 *  Compute this size of a packed 2-d polynomial coefficient array.
 */
inline int computePackedSize(int order) { return computePackedOffset(order + 1); }

/**
 *  Fill an array with integer powers of x, so @f$$r[n] == r^n@f$.
 *
 *  When multiple powers are needed, this should be signficantly faster than
 *  repeated calls to std::pow().
 */
void computePowers(Eigen::VectorXd& r, double x);

/**
 *  Return an array with integer powers of x, so @f$$r[n] == r^n@f$.
 *
 *  When multiple powers are needed, this should be signficantly faster than
 *  repeated calls to std::pow().
 */
Eigen::VectorXd computePowers(double x, int n);

/**
 *  A class that computes binomial coefficients up to a certain power.
 *
 *  The binomial coefficient is defined as:
 *  @f[
 *     \left(\begin{array}{ c }
 *       n
 *       k
 *     \end{array}right\)
 *     = \frac{n!}{k!(n-k)!}
 *  @f]
 *  with both @f$n@f$ and @f$k@f$ nonnegative integers and @f$k \le n@f$
 *
 *  This class uses recurrence relations to avoid computing factorials directly,
 *  making it both more efficient and numerically stable.
 */
class BinomialMatrix {
public:
    /**
     *  Construct an object that can compute binomial coefficients with @f$n@f$
     *  up to and including the given value.
     */
    explicit BinomialMatrix(int const nMax) { extend(nMax); }

    /**
     *  Return the binomial coefficient.
     *
     *  No error checking is performed; the behavior of this method is is
     *  undefined if the given values do not satisfy
     *  @code
     *  n <= nMax && k <= n && n >=0 && k >= 0
     *  @endcode
     */
    double operator()(int n, int k) const { return getMatrix()(n, k); }

private:
    static void extend(int const n);

    static Eigen::MatrixXd& getMatrix();
};

}  // namespace detail
}  // namespace astrom
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_ASTROM_DETAIL_polynomialUtils_h_INCLUDED