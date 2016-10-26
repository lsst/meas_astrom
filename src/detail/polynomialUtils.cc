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

#include <mutex>
#include "lsst/meas/astrom/detail/polynomialUtils.h"

namespace lsst { namespace meas { namespace astrom { namespace detail {

void computePowers(Eigen::VectorXd & r, double x) {
    r[0] = 1.0;
    for (int i = 1; i < r.size(); ++i) {
        r[i] = r[i - 1] * x;
    }
}

Eigen::VectorXd computePowers(double x, int n) {
    Eigen::VectorXd r(n + 1);
    computePowers(r, x);
    return r;
}

void BinomialMatrix::extend(int const n) {
    static std::mutex mutex;
    auto & old = getMatrix();
    int const m = old.rows() - 1;
    if (n <= m) return;
    Eigen::MatrixXd updated = Eigen::MatrixXd::Zero(n + 1, n + 1);
    updated.topLeftCorner(old.rows(), old.cols()) = old;
    for (int i = m + 1; i <= n; ++i) {
        updated(i, 0) = 1.0;
        updated(i, i) = 1.0;
        for (int j = 1; j < i; ++j) {
            updated(i, j) = updated(i - 1, j - 1) *
                (static_cast<double>(i) / static_cast<double>(j));
        }
    }
    std::unique_lock<std::mutex> lock(mutex);
    old.swap(updated);
}

Eigen::MatrixXd & BinomialMatrix::getMatrix() {
    static Eigen::MatrixXd it = Eigen::MatrixXd::Constant(2, 2, 1.0);
    return it;
}

}}}} // namespace lsst::meas::astrom::detail
