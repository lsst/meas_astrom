from __future__ import absolute_import, division, print_function

import lsst.afw.table as afwTable

__all__ = ["plotMatches"]

def plotMatches(
    matches,
    refCat=None,
    sourceCat=None,
    refMarker="+",
    refColor="b",
    sourceMarker="x",
    sourceColor="g",
    matchColor="b"):
    """Plot reference objects, sources and matches

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
    refCentroidKey = afwTable.Point2DKey(refSchema["centroid"])
    srcSchema = matches[0][1].schema
    srcCentroidKey = afwTable.Point2DKey(srcSchema["slot_Centroid"])

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

    def plotRefCircles(matches, color):
        refXYArr = [match[0].get(refCentroidKey) for match in matches]
        refXArr, refYArr = zip(*refXYArr)
        plt.plot(refXArr, refYArr, "o", mec=color, mfc="none", ms=10,)

    if callable(matchColor):
        # different matches have different colors
        matchColorArr = [matchColor(match) for match in matches]

        # plot circles for each color separately
        matchColorSet = set(matchColorArr)
        for color in matchColorSet:
            subMatches = [match for match in matches if matchColor(match) == color]
            plotRefCircles(subMatches, color=color)
    else:
        matchColorArr = [matchColor]*len(refXYArr)
        plotRefCircles(matches, color=matchColor)

    lineSegData = makeLineSegmentData(refXYArr, srcXYArr, matchColorArr)
    plt.plot(*lineSegData)

    plt.show()
