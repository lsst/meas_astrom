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

__all__ = ["sourceMatchStatistics"]

import numpy as np


def sourceMatchStatistics(matchList, log=None):
    """Compute statistics on the accuracy of a wcs solution, using a precomputed list
    of matches between an image and a catalogue

    Input:
    matchList is a lsst::afw::detection::SourceMatch object

    Output:
    A dictionary storing the following quanities
    meanOfDiffInPixels          Average distance between image and catalogue position (in pixels)
    rmsOfDiffInPixels           Root mean square of distribution of distances
    quartilesOfDiffInPixels     An array of 5 values giving the boundaries of the quartiles of the
                                distribution.
    """

    size = len(matchList)
    if size == 0:
        raise ValueError("matchList contains no elements")

    dist = np.zeros(size)
    i = 0
    for match in matchList:
        catObj = match.first
        srcObj = match.second

        cx = catObj.getXAstrom()
        cy = catObj.getYAstrom()

        sx = srcObj.getXAstrom()
        sy = srcObj.getYAstrom()

        dist[i] = np.hypot(cx-sx, cy-sy)
        i = i+1

    dist.sort()

    quartiles = []
    for f in (0.25, 0.50, 0.75):
        i = int(f*size + 0.5)
        if i >= size:
            i = size - 1
        quartiles.append(dist[i])

    values = {}
    values['diffInPixels_Q25'] = quartiles[0]
    values['diffInPixels_Q50'] = quartiles[1]
    values['diffInPixels_Q75'] = quartiles[2]
    values['diffInPixels_mean'] = dist.mean()
    values['diffInPixels_std'] = dist.std()

    return values
