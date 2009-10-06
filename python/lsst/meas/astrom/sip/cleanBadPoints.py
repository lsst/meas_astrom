
import os
import matplotlib.pyplot as mpl
import numpy as np

import eups
import lsst.meas.astrom.net as net
import lsst.afw.detection as det

import sipLib as sip

def clean(srcMatch, order=1, nsigma=3):
    """Remove bad points from srcMatch
    
    Input:
    srcMatch : std::vector<det::SourceMatch>
    order:      Order of polynomial to use in robust fitting
    nsigma:    Sources more than this far away from the robust best fit
                polynomial are removed
                
    Return:
    std::vector<det::SourceMatch> of the good data points
    """
    
    getx = lambda x: x[0].getXAstrom()
    getdx = lambda x: getx(x) - x[1].getXAstrom()
    
    x = np.array(map(getx, srcMatch))
    dx = np.array(map(getdx, srcMatch))
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

    plot=False
    
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
        deviance = np.fabs( (y - f) /s)
        newidx = np.where(deviance < nsigma)

        if plot:
            mpl.plot(rx, ry, 'b-')        
            mpl.plot(rx, ry, 'bs')        
            mpl.plot(rx, fit, 'ms')        
            mpl.plot(rx, fit, 'm-')        
            for i in range(len(rx)):
                print rx[i], ry[i], lsf.valueAt(rx[i])

            mpl.plot(x[newidx], y[newidx], 'bs')
        
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
    return newidx[0]
    

def chooseRx(x, idx, order):
    rSize = len(idx)/float(order+1)  #Note, a floating point number
    rx = np.zeros((order+1))
    
    for i in range(order+1):
        #rng = idx[range(int(rSize*i),int(rSize*(i+1))]
        rng = range(int(rSize*i), int(rSize*(i+1)))
        rx[i] = np.mean(x[idx[rng]])
    return rx
    

def chooseRy(y, idx, order):
    rSize = len(idx)/float(order+1)  #Note, a floating point number
    ry = np.zeros((order+1))
    
    for i in range(order+1):
        #rng = idx[range(int(rSize*i),int(rSize*(i+1))]
        rng = range(int(rSize*i), int(rSize*(i+1)))
        ry[i] = np.median(y[idx[rng]])
    return ry

