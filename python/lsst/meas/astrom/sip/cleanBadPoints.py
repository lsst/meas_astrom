
import os
import numpy as np

import eups
import lsst.afw.detection as det

import sipLib as sip


def clean(srcMatch, wcs, order=3, nsigma=3):
    """Remove bad points from srcMatch
    
    Input:
    srcMatch : std::vector<det::SourceMatch>
    order:      Order of polynomial to use in robust fitting
    nsigma:    Sources more than this far away from the robust best fit
                polynomial are removed
                
    Return:
    std::vector<det::SourceMatch> of the good data points
    """
    
    #Cataloge ra and dec
    raFunc = lambda x: (x.first).getRa()
    catRa = map(raFunc, srcMatch)
    decFunc = lambda x: (x.first).getDec()
    catDec = map(decFunc, srcMatch)
    
    #Wcs ra and dec
    wcsX = np.zeros(len(catRa))
    wcsY = np.zeros(len(catRa))
    for i in range(len(catRa)):
        tmp1, tmp2 = wcs.raDecToXY(catRa[i], catDec[i])
        wcsX[i] = tmp1
        wcsY[i] = tmp2


    getx = lambda x: (x.second).getXAstrom()
    x = np.array(map(getx, srcMatch))
    dx = x - wcsX
    s = np.zeros( (len(x)) ) + .1
    
    idx = indicesOfGoodPoints(x, dx, s, order=order, nsigma=nsigma)

    clean = []
    for i in idx:
        clean.append(srcMatch[i])
    return clean
            
    

def indicesOfGoodPoints(x, y, s, order=1, nsigma=3, maxiter=100):
    """Return a list of indices in the range [0, len(x)]
    of points that lie less than nsigma away from the robust
    best fit polynomial
    """
    
    #Indices of elements of x sorted in order of increasing value
    idx = x.argsort()
    newidx=[]
    flag=True
    niter = 0
    while flag and niter < maxiter:
        rx = chooseRx(x, idx, order)
        ry = chooseRy(y, idx, order)
        rs = np.ones((len(rx)))
          
        lsf = sip.LeastSqFitter1dPoly(list(rx), list(ry), list(rs), order)
        fit = map(lambda x: lsf.valueAt(x), rx)

        f = map(lambda x: lsf.valueAt(x), x)
        
        sigma = (y-f).std()
        deviance = np.fabs( (y - f) /sigma)
        newidx = np.where(deviance < nsigma)

        if False:
            import matplotlib.pyplot as mpl
            mpl.plot(x, y, 'ks')
            mpl.plot(rx, ry, 'b-')        
            mpl.plot(rx, ry, 'bs')        
            mpl.plot(rx, fit, 'ms')        
            mpl.plot(rx, fit, 'm-')        

            #mpl.plot(x[newidx], y[newidx], 'rs')
            mpl.show()
        
        #If we haven't culled any points we're finished cleaning
        if len(newidx) == len(idx):
            flag = True
        
        #Otherwise, try another iteration    
        niter=niter+1
        flag=False
    
    #We get here because we either a) stopped finding bad points
    #or b) ran out of iterations. Eitherway, we just return our
    #list of indices of good points.
    #Somewhere along the line, newidx becomes a tuple, only the 0th
    #elt of which is the list we need
    if len(newidx[0]) == 0:
        raise RuntimeError("All points cleaned out. This is probably a bug")
    return newidx[0]
    

def chooseRx(x, idx, order):
    """Create order+1 values of the ordinate based on the median of groups of elements of x"""
    rSize = len(idx)/float(order+1)  #Note, a floating point number
    rx = np.zeros((order+1))
    
    for i in range(order+1):
        rng = range(int(rSize*i), int(rSize*(i+1)))
        rx[i] = np.mean(x[idx[rng]])
    return rx
    

def chooseRy(y, idx, order):
    """Create order+1 values of the ordinate based on the median of groups of elements of y"""
    rSize = len(idx)/float(order+1)  #Note, a floating point number
    ry = np.zeros((order+1))
    
    for i in range(order+1):
        rng = range(int(rSize*i), int(rSize*(i+1)))
        ry[i] = np.median(y[idx[rng]])
    return ry

