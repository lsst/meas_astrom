
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet

def noDistort(src):
    """Do no distortion. Used for sanity checking"""
    
    out = afwDet.Source(src)
    return out


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


def linearYDistort(src, frac=.001):
    """Increase the yAstrom value in a Source object by frac. E.g
    src.xAstrom = 1000 --> 1001 if frac=.001
    
    Input:
    src     A Source object
    frac    How much to change Y by
    
    Output:
    A deep copy of src, with the value of yAstrom changed
    """

    out = afwDet.Source(src)
    out.setYAstrom( out.getYAstrom()*(1+frac) )
    return out


def quadraticDistortX(src, frac=1e-6):
    """Distort image by terms with power <=2
    i.e y, y^2, x, xy, x^2
    """
    
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val = x**2
    
    out.setXAstrom(x + val*frac)
    out.setYAstrom(y)
    return out

def quadraticDistort(src, frac=1e-6):
    """Distort image by terms with power <=2
    i.e y, y^2, x, xy, x^2
    """
    
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val =  y + 2*y**2
    val += 3*x + 4*x*y
    val += x**2
    
    out.setXAstrom(x + val*frac)
    out.setYAstrom(y)
    return out
    
def T2DistortX(src, frac=1e-6):
    """Distort image by a 2nd order Cheby polynomial"""

    out = afwDet.Source(src)
    x = src.getXAstrom()
    val = 2*(x**2) -1
    out.setXAstrom( x + frac*val)
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
