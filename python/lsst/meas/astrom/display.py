#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
from __future__ import absolute_import, division, print_function
import math

import numpy as np

import lsst.afw.display as afwDisplay
from lsst.afw.image import ExposureF
from lsst.afw.table import Point2DKey

__all__ = ["displayAstrometry", "plotAstrometry"]

def displayAstrometry(refCat=None, sourceCat=None, distortedCentroidKey=None, bbox=None, exposure=None,
                      matches=None, frame=1, title="", pause=True):
    """Show an astrometry debug image

    - reference objects in refCat are shown as red X
    - sources in sourceCat are shown as green +
    - distorted sources in sourceCat (position given by distortedCentroidKey) are shown as green o
    - matches are shown as a yellow circle around the source and a yellow line
        connecting the reference object and source
    - if both exposure and bbox are None, no image is displayed

    @param[in] refCat  reference object catalog; must have fields "centroid_x" and "centroid_y"
    @param[in] sourceCat  source catalog; must have field "slot_Centroid_x" and "slot_Centroid_y"
    @param[in] distortedCentroidKey  key for sourceCat with field to use for distorted positions, or None
    @param[in] exposure  exposure to display, or None for a blank exposure
    @param[in] bbox  bounding box of exposure; used if exposure is None for a blank image
    @param[in] matches  list of matches (an lsst.afw.table.ReferenceMatchVector), or None
    @param[in] frame  frame number for ds9 display
    @param[in] title  title for ds9 display
    @param[in] pause  pause for inspection of display? This is done by dropping into pdb.
    """
    disp = afwDisplay.getDisplay(frame)

    if exposure is not None:
        disp.mtv(exposure, title=title)
    elif bbox is not None:
        disp.mtv(exposure=ExposureF(bbox), title=title)

    with disp.Buffering():
        if refCat is not None:
            refCentroidKey = Point2DKey(refCat.schema["centroid"])
            for refObj in refCat:
                rx, ry = refObj.get(refCentroidKey)
                disp.dot("x", rx, ry, size=10, ctype=afwDisplay.RED)

        if sourceCat is not None:
            sourceCentroidKey = Point2DKey(sourceCat.schema["slot_Centroid"])
            for source in sourceCat:
                sx, sy = source.get(sourceCentroidKey)
                disp.dot("+", sx,  sy, size=10, ctype=afwDisplay.GREEN)
                if distortedCentroidKey is not None:
                    dx, dy = source.get(distortedCentroidKey)
                    disp.dot("o", dx, dy, size=10, ctype=afwDisplay.GREEN)
                    disp.line([(sx, sy), (dx, dy)], ctype=afwDisplay.GREEN)

        if matches is not None:
            refCentroidKey = Point2DKey(matches[0].first.schema["centroid"])
            sourceCentroidKey = Point2DKey(matches[0].second.schema["slot_Centroid"])
            radArr = np.ndarray(len(matches))

            for i, m in enumerate(matches):
                refCentroid = m.first.get(refCentroidKey)
                sourceCentroid = m.second.get(sourceCentroidKey)
                radArr[i] = math.hypot(*(refCentroid - sourceCentroid))
                sx, sy = sourceCentroid
                disp.dot("o", sx, sy, size=10, ctype=afwDisplay.YELLOW)
                disp.line([refCentroid, sourceCentroid], ctype=afwDisplay.YELLOW)

            print("<match radius> = %.4g +- %.4g [%d matches]" %
                (radArr.mean(), radArr.std(), len(matches)))

    if pause:
        print("Dropping into debugger to allow inspection of display. Type 'continue' when done.")
        import pdb
        pdb.set_trace()

def plotAstrometry(
    matches,
    refCat=None,
    sourceCat=None,
    refMarker="x",
    refColor="r",
    sourceMarker="+",
    sourceColor="g",
    matchColor="y"):
    """Plot reference objects, sources and matches

    By default:
    - reference objects in refCat are shown as red X
    - sources in sourceCat are shown as green +
    - matches are shown as a yellow circle around the source and a line
        connecting the reference object to the source

    @param[in] matches  list of matches
    @param[in] refCat  reference object catalog, or None to not plot reference objects
    @param[in] sourceCat  source catalog, or None to not plot sources
    @param[in] refMarker  pyplot marker for reference objects
    @param[in] refColor  pyplot color for reference objects
    @param[in] sourceMarker  pyplot marker for sources
    @param[in] sourceColor  pyplot color for sources
    @param[in] matchColor  color for matches; can be a constant
        or a function taking one match and returning a string
    """
    # delay importing plt to give users a chance to set the backend before calling this function
    import matplotlib.pyplot as plt
    refSchema = matches[0][0].schema
    refCentroidKey = Point2DKey(refSchema["centroid"])
    srcSchema = matches[0][1].schema
    srcCentroidKey = Point2DKey(srcSchema["slot_Centroid"])

    if refCat is not None:
        refXArr, refYArr = zip(*[ref.get(refCentroidKey) for ref in refCat])
        plt.plot(refXArr, refYArr, marker=refMarker, color=refColor, linestyle="")

    if sourceCat is not None:
        srcXArr, srcYArr = zip(*[src.get(srcCentroidKey) for src in sourceCat])
        plt.plot(srcXArr, srcYArr, marker=sourceMarker, color=sourceColor, linestyle="")

    def makeLineSegmentData(refXYArr, srcXYArr, colorArr):
        """Make a list of line segement data

        This is used to draw line segements between ref and src positions in the specified color

        The returned data has the format:
         [(refX0, srcX0), (refY0, srcY0), color0, (refX1, srcX1), (refY1, srcY1), color1,...]
        """
        if len(refXYArr) != len(srcXYArr):
            raise RuntimeError("len(refXYArr) = %d != %d = len(srcXYArr)" %
                (len(refXYArr), len(srcXYArr)))
        if len(refXYArr) != len(colorArr):
            raise RuntimeError("len(refXYArr) = %d != %d = len(colorArr)" %
                (len(refXYArr), len(colorArr)))

        refXArr, refYArr = zip(*refXYArr)
        srcXArr, srcYArr = zip(*srcXYArr)
        refSrcXArr = zip(refXArr, srcXArr)
        refSrcYArr = zip(refYArr, srcYArr)
        dataList = []
        for xycolor in zip(refSrcXArr, refSrcYArr, colorArr):
            for val in xycolor:
                dataList.append(val)
        return dataList

    refXYArr, srcXYArr = \
        zip(*[(match[0].get(refCentroidKey), match[1].get(srcCentroidKey)) for match in matches])

    def plotSourceCircles(matches, color):
        srcXYArr = [match[1].get(srcCentroidKey) for match in matches]
        srcXArr, srcYArr = zip(*srcXYArr)
        plt.plot(srcXArr, srcYArr, "o", mec=color, mfc="none", ms=10,)

    if callable(matchColor):
        # different matches have different colors
        matchColorArr = [matchColor(match) for match in matches]

        # plot circles for each color separately
        matchColorSet = set(matchColorArr)
        for color in matchColorSet:
            subMatches = [match for match in matches if matchColor(match) == color]
            plotSourceCircles(subMatches, color=color)
    else:
        matchColorArr = [matchColor]*len(refXYArr)
        plotSourceCircles(matches, color=matchColor)

    lineSegData = makeLineSegmentData(refXYArr, srcXYArr, matchColorArr)
    plt.plot(*lineSegData)

    plt.show()
