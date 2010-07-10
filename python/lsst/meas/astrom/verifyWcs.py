import math, sys
import numpy

import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg

def checkMatches(srcMatchSet, exposure, log=None):
    if not exposure:
        return {}
    
    im = exposure.getMaskedImage().getImage()
    width, height = im.getWidth(), im.getHeight()
    nx, ny = 3, 3
    w, h = width//nx, height//ny

    if w == 0:
        w = 1
    while nx*w < width:
        w += 1

    if h == 0:
        h = 1
    while ny*h < height:
        h += 1
        
    cellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0, 0), width, height), w, h)
    #
    # Populate cellSet
    #
    i = -1
    for m in srcMatchSet:
        i += 1
        
        src = m.second
        csrc = afwDetection.Source()
        csrc.setId(i)
        csrc.setXAstrom(src.getXAstrom())
        csrc.setYAstrom(src.getYAstrom())

        try:
            cellSet.insertCandidate(measAlg.PsfCandidateF(csrc, exposure.getMaskedImage()))
        except Exception, e:
            log.log(log.WARN, str(e))

    ncell = len(cellSet.getCellList())
    nobj = numpy.ndarray(ncell, dtype='i')
        
    for i in range(ncell):
        cell = cellSet.getCellList()[i]

        nobj[i] = cell.size()

        dx = numpy.ndarray(cell.size())
        dy = numpy.ndarray(cell.size())   

        j = 0
        for cand in cell:
            #
            # Swig doesn't know that we're a SpatialCellImageCandidate;  all it knows is that we have
            # a SpatialCellCandidate so we need an explicit (dynamic) cast
            #
            cand = measAlg.cast_PsfCandidateF(cand)
            
            mid = cand.getSource().getId()
            dx[j] = srcMatchSet[mid].first.getXAstrom() - srcMatchSet[mid].second.getXAstrom()
            dy[j] = srcMatchSet[mid].first.getYAstrom() - srcMatchSet[mid].second.getYAstrom()

            j += 1

        if log:
            log.log(log.DEBUG, "%s %-30s  %8s  dx,dy = %5.2f,%5.2f  rms_x,y = %5.2f,%5.2f" % \
                    (cell.getLabel(), cell.getBBox(), ("nobj=%d" % cell.size()),
                     dx.mean(), dy.mean(), dx.std(), dy.std()))


    nobj.sort()

    values = {}
    values["minObjectsPerCell"] = int(nobj[0]) # otherwise it's a numpy integral type
    values["maxObjectsPerCell"] = int(nobj[-1])
    values["meanObjectsPerCell"] = nobj.mean()
    values["stdObjectsPerCell"] = nobj.std()

    return values
