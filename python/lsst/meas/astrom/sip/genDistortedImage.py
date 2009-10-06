
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet


def linearXDistort(src, frac=.001):
    """Increase the xAstrom value in a Source object by frac. E.g
    src.xAstrom = 1000 --> 1001 if frac=.001
    
    Input:
    src     A Source object
    frac    How much to change X by
    
    Output:
    A deep copy of src, with the value of xAstrom changed
    """
    

    out = afwDet.Source(src)

       
    out.setXAstrom( out.getXAstrom()*(1+frac) )
    return out



def distortList(srcList, function):
    """Create a copy of srcList, and apply function to distort the 
    values of xAstrom and yAstrom. 
    
    Input:
    srcList     a SourceSet object
    function:   A function that does a deep copy of a single Source
    """
    
    out = afwDet.SourceSet()
    
    for src in srcList:
        out.push_back( function(src) )

    return out
