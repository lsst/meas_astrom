
import numpy as np

def sourceMatchStatistics(matchList, log=None):
    """ Compute statistics on the accuracy of a wcs solution, using a precomputed list 
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
        
        dist[i] = np.hypot( cx-sx, cy-cx)
        i = i+1
        
    
    dist.sort()

    quartile = []
    quartile.append(dist[0])
    
    i = int(size/4.)
    quartile.append(dist[i])
    
    i = int(size/2.)
    quartile.append(dist[i])
    
    i = int(3*size/4.)
    quartile.append(dist[i])
    
    quartile.append(dist[-1])
    
    returnObject = dict()
    returnObject['quartilesOfDiffInPixels'] = quartile
    returnObject['meanOfDiffInPixels'] = dist.mean()
    returnObject['rmsOfDiffInPixels'] = dist.std()
    
    
    return returnObject
    
