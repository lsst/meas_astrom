/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testlsf1d

// The boost unit test header
#include "boost/test/unit_test.hpp"

using namespace std;

#include <vector>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <memory>
#include "Eigen/Core"
#include "lsst/afw/math/FunctionLibrary.h"

#include "lsst/meas/astrom/sip/LeastSqFitter1d.h"

namespace sip = lsst::meas::astrom::sip;
namespace math = lsst::afw::math;

BOOST_AUTO_TEST_CASE(fitLine) {
    vector<double> x;
    vector<double> y;
    vector<double> s;

    for (int i = 0; i < 7; ++i) {
        x.push_back((double)i);
        y.push_back((double)3 + 2 * i);
        s.push_back((double)1);
    }

    int order = 2;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);

    Eigen::VectorXd par = lsf.getParams();

    BOOST_CHECK(par.size() == order);
    BOOST_CHECK_CLOSE(par(0), 3., .001);
    BOOST_CHECK_CLOSE(par(1), 2., .001);
}

BOOST_AUTO_TEST_CASE(fitQuadratic) {
    vector<double> x;
    vector<double> y;
    vector<double> s;

    for (int i = 0; i < 7; ++i) {
        x.push_back((double)i);
        y.push_back((double)4 + i * (3 + i * 2));
        s.push_back((double)1);
    }

    int order = 3;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);

    math::PolynomialFunction1<double> f = lsf.getBestFitFunction();

    BOOST_CHECK((int)f.getNParameters() == order);
    BOOST_CHECK_CLOSE(f.getParameter(0), 4., .001);
    BOOST_CHECK_CLOSE(f.getParameter(1), 3., .001);
    BOOST_CHECK_CLOSE(f.getParameter(2), 2., .001);

    BOOST_CHECK_CLOSE(lsf.valueAt(1.), 9., .001);
}

// Check a more extreme case where there's a big difference between the x and y values
BOOST_AUTO_TEST_CASE(fitLinear2) {
    vector<double> x;
    vector<double> y;
    vector<double> s;

    x.push_back(631.920343622);
    x.push_back(1000.02009955);
    x.push_back(1205.63487689);
    x.push_back(1427.86231185);

    y.push_back(0.622249643993);
    y.push_back(1.01887634344);
    y.push_back(1.1982950985);
    y.push_back(1.4294170431);

    s.push_back(1);
    s.push_back(1);
    s.push_back(1);
    s.push_back(1);

    int order = 2;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);
    Eigen::VectorXd par = lsf.getParams();
    assert(par[0] == par[0]);  // stop compiler whining about par not being used
}

BOOST_AUTO_TEST_CASE(fitLinear3) {
    vector<double> x;
    vector<double> y;
    vector<double> s;

    x.push_back(689.301136505);
    x.push_back(1112.8573687);
    x.push_back(1386.67168477);

    y.push_back(0.66911456573);
    y.push_back(1.1147439759);
    y.push_back(1.39597284177);

    s.push_back(1);
    s.push_back(1);
    s.push_back(1);

    int order = 2;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);
    Eigen::VectorXd par = lsf.getParams();

    BOOST_CHECK_CLOSE(par(0), -0.0488399, .001);
    BOOST_CHECK_CLOSE(par(1), 0.00104313, .001);
}

BOOST_AUTO_TEST_CASE(fitQuadratic2) {
    vector<double> x;
    vector<double> y;
    vector<double> s;

    x.push_back(628.857680996);
    x.push_back(995.008255088);
    x.push_back(1203.39412154);
    x.push_back(1425.1404727);

    y.push_back(0.61672822987);
    y.push_back(1.01887634344);
    y.push_back(1.19679830465);
    y.push_back(1.42873084062);

    s.push_back(1);
    s.push_back(1);
    s.push_back(1);
    s.push_back(1);

    int order = 3;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);
    Eigen::VectorXd par = lsf.getParams();
    assert(par[0] == par[0]);  // stop compiler whining about par not being used

    for (unsigned int i = 0; i != x.size(); ++i) {
        // printf("%.3f %.3f %.3f \n", x[i], y[i], lsf.valueAt(x[i]));
        BOOST_CHECK_CLOSE(y[i], lsf.valueAt(x[i]), 10);
    }
}

BOOST_AUTO_TEST_CASE(errorbars) {
    vector<double> x;
    vector<double> y;
    vector<double> s;

    x.push_back(1);
    x.push_back(2);
    x.push_back(3);

    y.push_back(1);
    y.push_back(2);
    y.push_back(3);

    s.push_back(1);
    s.push_back(1);
    s.push_back(1);

    int order = 2;
    sip::LeastSqFitter1d<math::PolynomialFunction1<double> > lsf(x, y, s, order);
    Eigen::VectorXd par = lsf.getParams();
    Eigen::VectorXd err = lsf.getErrors();

    // Boost doesn't do well at checking for zero
    BOOST_CHECK_CLOSE(par[0] + 1, 1., 1e-6);
    BOOST_CHECK_CLOSE(par[1], 1., 1e-6);

    // Calculated by hand
    BOOST_CHECK_CLOSE(err[0], 1.52752523165195, 1e-6);
    BOOST_CHECK_CLOSE(err[1], 0.7071067811865481, 1e-6);
}
