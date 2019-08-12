# This file is part of meas_astrom.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["displayAstrometry", "plotAstrometry"]

import math

import numpy as np

import lsst.afw.display as afwDisplay
from lsst.afw.image import ExposureF
from lsst.afw.table import Point2DKey


def displayAstrometry(refCat=None, sourceCat=None, distortedCentroidKey=None, bbox=None, exposure=None,
                      matches=None, frame=1, title="", pause=True):
    """Show an astrometry debug image.

    Parameters
    ----------
    refCat : `lsst.afw.table.SimpleCatalog`
        reference object catalog; must have fields "centroid_x" an
        "centroid_y"
    sourceCat : `lsst.afw.table.SourceCatalg`
        source catalog; must have field "slot_Centroid_x" and "slot_Centroid_y"
    distortedCentroidKey : `lsst.afw.table.Key`
        key for sourceCat with field to use for distorted positions
    exposure : `lsst.afw.image.Exposure`
        exposure to display
    bbox : `lsst.geom.Box2I`
        bounding box of exposure; Used if the exposure is `None`
    matches : `list` of `lsst.afw.table.ReferenceMatch`
        List of matched objects
    frame : `int`
        frame number for display
    title : `str`
        title for display
    pause : `bool`
        pause for inspection of display? This is done by dropping into pdb.

    Notes
    -----

    - reference objects in refCat are shown as red X
    - sources in sourceCat are shown as green +
    - distorted sources in sourceCat (position given by distortedCentroidKey)
      are shown as green o
    - matches are shown as a yellow circle around the source and a yellow line
      connecting the reference object and source
    - if both exposure and bbox are `None`, no image is displayed

    """
    disp = afwDisplay.Display(frame=frame)
    disp.scale('asinh', 'zscale')

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
                disp.dot("+", sx, sy, size=10, ctype=afwDisplay.GREEN)
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
        print(f"""Meaning of symbols on display:
  Red x:    Reference catalogue [{len(refCat)}]
  Green +:  Source catalogue    [{len(sourceCat)}]""")
        if distortedCentroidKey is not None:
            print("""\
  Green o:  Source position allowing for distortion (with line to undistorted position)""")
        if matches is not None:
            print("""\
  Yellow o: Matched source""")

        print("")
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
    matchColor="y"
):
    """Plot reference objects, sources and matches

    Parameters
    ----------
    matches : `list` of `lsst.afw.table.ReferenceMatch`
        list of matches
    refCat : `lsst.afw.table.SimpleCatalog`
        reference object catalog, or None to not plot reference objects
    sourceCat : `lsst.afw.table.SourceCatalog`
        source catalog, or None to not plot sources
    refMarker : `str`
        pyplot marker for reference objects
    refColor : `str`
        pyplot color for reference objects
    sourceMarker : `str`
        pyplot marker for sources
    sourceColor : `str`
        pyplot color for sources
    matchColor : `str`
        color for matches; can be a constant
        or a function taking one match and returning a string

    Notes
    -----
    By default:

    - reference objects in refCat are shown as red X
    - sources in sourceCat are shown as green +
    - matches are shown as a yellow circle around the source and a line
      connecting the reference object to the source
    """
    # delay importing plt to give users a chance to set the backend before calling this function
    import matplotlib.pyplot as plt
    refSchema = matches[0][0].schema
    refCentroidKey = Point2DKey(refSchema["centroid"])
    srcSchema = matches[0][1].schema
    srcCentroidKey = Point2DKey(srcSchema["slot_Centroid"])

    if refCat is not None:
        refXArr, refYArr = list(zip(*[ref.get(refCentroidKey) for ref in refCat]))
        plt.plot(refXArr, refYArr, marker=refMarker, color=refColor, linestyle="")

    if sourceCat is not None:
        srcXArr, srcYArr = list(zip(*[src.get(srcCentroidKey) for src in sourceCat]))
        plt.plot(srcXArr, srcYArr, marker=sourceMarker, color=sourceColor, linestyle="")

    def makeLineSegmentData(refXYArr, srcXYArr, colorArr):
        """Make a list of line segement data

        This is used to draw line segements between ref and src positions in the specified color

        Notes
        -----
        The returned data has the format:
         [(refX0, srcX0), (refY0, srcY0), color0, (refX1, srcX1), (refY1, srcY1), color1,...]
        """
        if len(refXYArr) != len(srcXYArr):
            raise RuntimeError("len(refXYArr) = %d != %d = len(srcXYArr)" %
                               (len(refXYArr), len(srcXYArr)))
        if len(refXYArr) != len(colorArr):
            raise RuntimeError("len(refXYArr) = %d != %d = len(colorArr)" %
                               (len(refXYArr), len(colorArr)))

        refXArr, refYArr = list(zip(*refXYArr))
        srcXArr, srcYArr = list(zip(*srcXYArr))
        refSrcXArr = list(zip(refXArr, srcXArr))
        refSrcYArr = list(zip(refYArr, srcYArr))
        dataList = []
        for xycolor in zip(refSrcXArr, refSrcYArr, colorArr):
            for val in xycolor:
                dataList.append(val)
        return dataList

    refXYArr, srcXYArr = \
        list(zip(*[(match[0].get(refCentroidKey), match[1].get(srcCentroidKey)) for match in matches]))

    def plotSourceCircles(matches, color):
        srcXYArr = [match[1].get(srcCentroidKey) for match in matches]
        srcXArr, srcYArr = list(zip(*srcXYArr))
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
