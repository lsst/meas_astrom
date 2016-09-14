from __future__ import absolute_import
from builtins import range
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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


import os
import numpy as np

import lsst.afw.detection as det

from . import sipLib as sip


def clean(srcMatch, wcs, order=3, nsigma=3):
    """Remove bad points from srcMatch

    Input:
    srcMatch : std::vector<det::SourceMatch>
    order:      Order of polynomial to use in robust fitting
    nsigma:    Sources more than this far away from the robust best fit
                polynomial are removed

    Return:
    std::vector<det::SourceMatch> of the good data points
    """

    N = len(srcMatch)
    catX = np.zeros(N)
    #catY = np.zeros(N)
    for i in range(N):
        x, y = wcs.skyToPixel(srcMatch[i].first.getCoord())
        catX[i] = x
        #catY[i] = y

    # FIXME -- why does this only use X?

    x = np.array([s.second.getX() for s in srcMatch])
    dx = x - catX
    sigma = np.zeros_like(dx) + 0.1

    idx = indicesOfGoodPoints(x, dx, sigma, order=order, nsigma=nsigma)

    clean = []
    for i in idx:
        clean.append(srcMatch[i])
    return clean


def indicesOfGoodPoints(x, y, s, order=1, nsigma=3, maxiter=100):
    """Return a list of indices in the range [0, len(x)]
    of points that lie less than nsigma away from the robust
    best fit polynomial
    """

    # Indices of elements of x sorted in order of increasing value
    idx = x.argsort()
    newidx = []
    for niter in range(maxiter):
        rx = chooseRx(x, idx, order)
        ry = chooseRy(y, idx, order)
        rs = np.ones((len(rx)))

        lsf = sip.LeastSqFitter1dPoly(list(rx), list(ry), list(rs), order)
        fit = [lsf.valueAt(value) for value in rx]
        f = [lsf.valueAt(value) for value in x]

        sigma = (y-f).std()
        deviance = np.fabs((y - f) / sigma)
        newidx = np.flatnonzero(deviance < nsigma)

        if False:
            import matplotlib.pyplot as mpl
            mpl.plot(x, y, 'ks')
            mpl.plot(rx, ry, 'b-')
            mpl.plot(rx, ry, 'bs')
            mpl.plot(rx, fit, 'ms')
            mpl.plot(rx, fit, 'm-')
            #mpl.plot(x[newidx], y[newidx], 'rs')
            mpl.show()

        # If we haven't culled any points we're finished cleaning
        if len(newidx) == len(idx):
            break

    # We get here because we either a) stopped finding bad points
    # or b) ran out of iterations. Either way, we just return our
    # list of indices of good points.
    return newidx


def chooseRx(x, idx, order):
    """Create order+1 values of the ordinate based on the median of groups of elements of x"""
    rSize = len(idx)/float(order+1)  # Note, a floating point number
    rx = np.zeros((order+1))

    for i in range(order+1):
        rng = list(range(int(rSize*i), int(rSize*(i+1))))
        rx[i] = np.mean(x[idx[rng]])
    return rx


def chooseRy(y, idx, order):
    """Create order+1 values of the ordinate based on the median of groups of elements of y"""
    rSize = len(idx)/float(order+1)  # Note, a floating point number
    ry = np.zeros((order+1))

    for i in range(order+1):
        rng = list(range(int(rSize*i), int(rSize*(i+1))))
        ry[i] = np.median(y[idx[rng]])
    return ry
