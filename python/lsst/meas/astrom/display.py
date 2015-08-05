from __future__ import absolute_import, division, print_function

import math

import numpy as np

from lsst.afw.image import ExposureF
from lsst.afw.table import Point2DKey

__all__ = ["displayAstrometry", "plotAstrometry"]

def displayAstrometry(refCat=None, sourceCat=None, distortedCentroidKey=None, bbox=None, exposure=None,
                      matches=None, frame=1, title=""):
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
    """
    import lsst.afw.display.ds9 as ds9

    if exposure is not None:
        ds9.mtv(exposure, frame=frame, title=title)
    elif bbox is not None:
        ds9.mtv(exposure=ExposureF(bbox), frame=frame, title=title)

    with ds9.Buffering():
        if refCat is not None:
            refCentroidKey = Point2DKey(refCat.schema["centroid"])
            for refObj in refCat:
                rx, ry = refObj.get(refCentroidKey)
                ds9.dot("x", rx, ry, size=10, frame=frame, ctype=ds9.RED)

        if sourceCat is not None:
            sourceCentroidKey = Point2DKey(sourceCat.schema["slot_Centroid"])
            for source in sourceCat:
                sx, sy = source.get(sourceCentroidKey)
                ds9.dot("+", sx,  sy, size=10, frame=frame, ctype=ds9.GREEN)
                if distortedCentroidKey is not None:
                    dx, dy = source.get(distortedCentroidKey)
                    ds9.dot("o", dx, dy, size=10, frame=frame, ctype=ds9.GREEN)
                    ds9.line([(sx, sy), (dx, dy)], ctype=ds9.GREEN, frame=frame)

        if matches is not None:
            refCentroidKey = Point2DKey(matches[0].first.schema["centroid"])
            sourceCentroidKey = Point2DKey(matches[0].second.schema["slot_Centroid"])
            radArr = np.ndarray(len(matches))

            for i, m in enumerate(matches):
                refCentroid = m.first.get(refCentroidKey)
                sourceCentroid = m.second.get(sourceCentroidKey)
                radArr[i] = math.hypot(*(refCentroid - sourceCentroid))
                sx, sy = sourceCentroid
                ds9.dot("o", sx, sy, size=10, frame=frame, ctype=ds9.YELLOW)
                ds9.line([refCentroid, sourceCentroid], ctype=ds9.YELLOW)
                
            print("<match radius> = %.4g +- %.4g [%d matches]" %
                (radArr.mean(), radArr.std(), len(matches)))

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
