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

import math

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


def cubicDistortX(src, frac=1e-9):
    """Distort image by terms with power <=2
    i.e y, y^2, x, xy, x^2
    """
    
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val = x**3
    
    out.setXAstrom(x + val*frac)
    out.setYAstrom(y)
    return out


def manyTermX(src, frac=1e-9):
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val = x**3 - 2*x**2 + 4*x - 9
    
    out.setXAstrom(x + val*frac)
    out.setYAstrom(y)
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


def quadraticDistortY(src, frac=1e-6):
    """Distort image by terms with power <=2
    i.e y, y^2, x, xy, x^2
    """
    
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val = y**2
    
    out.setXAstrom(x)
    out.setYAstrom(y + val*frac)
    return out


def cubicDistortY(src, frac=1e-9):
    """Distort image by terms with power <=2
    i.e y, y^2, x, xy, x^2
    """
    
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val = x**3
    
    out.setXAstrom(x)
    out.setYAstrom(y + val*frac)
    return out


def manyTermY(src, frac=1e-9):
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val = y**3 - 2*y**2 + 4*y - 9
    
    out.setXAstrom(x)
    out.setYAstrom(y + val*frac)
    return out


def crossTerms1(src, frac=1e-11):
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val = x**3 - 2*x**2 #+ 4*x - 9
    
    out.setXAstrom(x)
    out.setYAstrom(y + val*frac)
    return out


def crossTerms2(src, frac=1e-11):
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    val = y**3 - 2*y**2 + 4*y - 9
    
    out.setXAstrom(x+ val*frac)
    out.setYAstrom(y)
    return out


def crossTerms3(src, frac=1e-9):
    out = afwDet.Source(src)
    x = out.getXAstrom()
    y = out.getYAstrom()
    valx = x**3 - 2*x**2 + 4*x - 9
    valy = y**3 - 2*y**2 + 4*y - 9
    
    out.setXAstrom(x + valy*frac)
    out.setYAstrom(y + valx*frac)
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

    maxDiff=0
    for i in range(len(srcList)):
        s = srcList[i]
        o = out[i]
        
        x1, y1 = s.getXAstrom(), s.getYAstrom()
        x2, y2 = o.getXAstrom(), o.getYAstrom()
        
        diff = math.hypot(x1-x2, y1-y2)
        maxDiff = max(diff, maxDiff)
        
    print "Max deviation is %e pixels" %(maxDiff)
    
    return out
