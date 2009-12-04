#!/usr/bin/evn python

import os
import math
import matplotlib.pyplot as mpl

import eups
import lsst.meas.astrom.net as net
import lsst.afw.image as afwImg
import lsst.afw.display.ds9 as ds9

import lsst.meas.astrom.sip as sip
import lsst.meas.astrom.sip.sourceSetIO as sourceSetIO
import lsst.meas.astrom.sip.cleanBadPoints as cleanBadPoints


import getSourceSet 


import pdb

################################################


def doBase(path, basename, gas, rerun=False):
    """Run on all the chips and exposures (e=>exposure?) for a given basename
    if reRun=True, looks for output of previous runs
    """
    

    for e in range(1):
        for c in range(3):
        
            #Get wcs for a single amplifier
            filename = "%s-e%d-c%03d-a%02d.sci" % (basename, e, c, 4)
            filename = os.path.join(path, "IPSD", "output", "sci", "%s-e%d" % (basename, e), filename)
            exp = afwImg.ExposureF(filename)
            wcs=exp.getWcs()
            print wcs.getOriginRaDec()

            if rerun:
                #Try to load info from a text file
                filename="%s-e%d-c%03d.xy.txt" % (basename, e, c)
                try:
                    srcSet = sourceSetIO.read(filename)
                except:
                    rerun=False

            if rerun==False:
                #Get source set by extracting sources from images
                srcSet = extractSourceSetForCCD(path, basename, e, c)
            
            if astromVerifyWcs(gas, srcSet, wcs):
                print "Astrometry solution found"
            else:
                print "Warning, no solution found. Trying blind solve"
                if astromBlindSolve(gas, srcSet):
                    print "Blind solution found"
                    print gas.getWcs().getOriginRaDec(), gas.getSolvedImageScale()
                else:
                    print "Failed again. Giving up"
                
            

def extractSourceSetForCCD(path, basename, e, c, threshold=150):
    """Extracts sources from every amplifier on a chip, writes to a file and returns result
    """

    #Create an list of star positions. Function changes the value of path
    path2=path
    a=range(7)  #Every amplifier
    srcSet = getSourceSet.makeSourceList(path2, basename, e, c, a, threshold, verbose=True)

    filename="%s-e%i-c%03i.xy.txt" %(basename, e, c)
    path=os.path.join(path,filename)
    sourceSetIO.write(srcSet, path)

    return srcSet

def plotSourceSet(exp, srcSet, offSetToAmplifier=True):

    mi=exp.getMaskedImage()
    x0=mi.getX0()
    y0=mi.getY0()
    for i in range(srcSet.size()):
        if offSetToAmplifier:
            x=srcSet[i].getXAstrom()-x0
            y=srcSet[i].getYAstrom()-y0
        else:
            x=srcSet[i].getXAstrom()
            y=srcSet[i].getYAstrom()
            
        ds9.dot('o', x, y, ctype=ds9.RED)


################################################

def loadGas():

    andDir =eups.productDir("astrometry_net_data")
    print andDir
    if andDir == None:
        raise RuntimeError("astrometry_net_data not set up")
        
    paf = os.path.join(eups.productDir("astrometry_net_data"), "metadata.paf")
    gas = net.GlobalAstrometrySolution(paf)
    return gas


xyfile="v695856-e0-c000-a00.sci.xy.txt"
ra0=215.604025685476
dec0=53.1595451514076
def astromSolve(gas, srcSet, ra0, dec0):
    
    gas.reset()
    gas.setLogLevel(3)
    
    gas.setMinimumImageScale(.1)
    gas.setMaximumImageScale(.2)
    gas.setStarlist(srcSet)
    gas.setNumBrightObjects(50)
    
    return gas.solve(ra0, dec0)


def astromVerifyWcs(gas, srcSet, wcs):
    
    gas.reset()
    gas.setLogLevel(0)
    
    gas.setMinimumImageScale(.1)
    gas.setMaximumImageScale(.2)
    gas.setStarlist(srcSet)
    gas.setNumBrightObjects(50)
    
    return gas.solve(wcs)


def astromBlindSolve(gas, srcSet):
    
    gas.reset()
    gas.setLogLevel(0)
    
    gas.setMinimumImageScale(.1)
    gas.setMaximumImageScale(.2)
    gas.setStarlist(srcSet)
    gas.setNumBrightObjects(50)
    
    return gas.solve()
    



"""
import cfhtdata as cf
import sourceSetIO as ssi
srcSet = ssi.read('v695833-e0-c001-a06.sci.xy.txt')
gas = cf.loadGas()
cf.astromSolve(gas, srcSet, cf.ra0, cf.dec0)
"""

    




    
if __name__ == "__main__":
    do()

    
    
    
