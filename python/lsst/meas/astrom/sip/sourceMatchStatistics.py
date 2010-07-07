
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
    values['diffInPixels_q25'] = quartiles[0]
    values['diffInPixels_Q50'] = quartiles[1]
    values['diffInPixels_Q75'] = quartiles[2]
    values['diffInPixels_mean'] = dist.mean()
    values['dDiffInPixels_std'] = dist.std()
    
    return values
    
