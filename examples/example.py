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


#example.py
#An example of how to create a Wcs for an image
#Fergal Mullally

#Usage: python example.py


import re
import os
import glob
import math
from datetime import datetime

import eups
import lsst.afw.image as afwImage
import lsst.afw.detection.detectionLib as detect

import lsst.meas.astrom.net as net
import lsst.meas.astrom.sip as sip



def loadXYFromFile(filename):
    """Load a list of positions from a file and store in a Source"""
    f= open(filename)
    
    s1=detect.SourceSet()
    i=0
    for line in f:
        #Split the row into an array
        line = re.sub("^\s+", "", line)
        elts = re.split("\s+", line)
        
        #Swig requires floats, but elts is an array of strings
        x=float(elts[0])
        y=float(elts[1])
        flux=float(elts[2])

        source = detect.Source()

        source.setSourceId(i)
        source.setXAstrom(x); source.setXAstromErr(0.1)
        source.setYAstrom(y); source.setYAstromErr(0.1)

        s1.append(source)
        
        i=i + 1
    f.close()
    
    return s1


#Create a new Global Solution object. The object is constructed using a policy file
#stored in the $ASTROMETRY_NET_DATA directory
metadataFile = os.path.join(eups.productDir("astrometry_net_data"), "metadata.paf")
gas = net.GlobalAstrometrySolution(metadataFile)

#An image is abstracted as SourceSet, i.e a list of Source objects. To speed the example
#we're pre-extracted a list of sources, and stored in in a text file
path = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")
starlist = loadXYFromFile(path)
#Add the image to the solver
gas.setStarlist(starlist)



#Set the platescale of the image in arcsec per pixel. If you're not sure
#you can set a range of platescales. This step is not mandatory, but 
#it can speed up the code dramatically.
opt=1
if opt == 1:
    gas.setImageScaleArcsecPerPixel(.185)
elif opt == 2:
    gas.setMinimumImageScale(.1)
    gas.setMaximumImageScale(.2)
    

#Find an astrometric solution. 
opt=3
success=0
ra = 334.332500 
dec = -17.3263611 
if opt == 1:
    #If you have no idea where on the sky this image is, you can solve
    #with no positional arguments. However, this can take 10 minutes 
    #or more, so be patient. There is also no garauntee that a solution will
    #be found, so always check the return value.
    success = gas.solve()
elif opt == 2:
    #With the platescale set, and a good initial guess, a solution can
    #be found in less than a second. Note that if your initial guesses
    #are too far off, then the solver will fail
    success = gas.solve(ra, dec)
elif opt == 3:
    #Position can also be passed as a PointD object
    pos = afwImage.PointD(ra, dec)
    success = gas.solve(pos)
elif opt == 4:
    #You can also pass in a wcs, which is used as an initial guess
    #success = gas.solve(wcs)
    pass
    
if success:
    print "Solution found"
    
    #Return a wcs object.
    wcs = gas.getWcs()

    #The next step is to fit the distortion in the image. The distortion
    #is expressed as SIP polynomials (Shupe et al. 2005)
    
    #The first step is to extract a position catalogue from the solver
    imgSizeInArcsec = 8*60  #No harm in having the catalogue too big
    cat = gas.getCatalogue(imgSizeInArcsec)
    
    #Now produce a list of matching pairs of objects between the image
    #and the catalogue
    matcher = sip.MatchSrcToCatalogue(cat, starlist, wcs)    
    matchList = matcher.getMatches()

    #Now that we know where the objects are, and where the catalogue
    #says they should be, we can compute the distortion terms.
    opt=1
    if opt==1:   
        
        #The easiest way to do this is to specify how good a wcs you 
        #want (i.e what is the typical difference between where wcs says
        #the image is in radec space, and where the cataloge says it should be)
        #Bear in mind, that there is some jitter in the positions due to 
        #centroiding error which the SIP corrections can't account for,
        #so don't set your tolerance too low.
        maxScatter=1 #Arcseconds
        maxSipOrder=5
        
        sipObject = None
        sipObject = sip.CreateWcsWithSip(matchList, wcs, maxScatter, maxSipOrder)
    else:
        #Alternatively, you can specify the number of terms in the correction directly
        order=4
        sipObject = sup.CreateWcsWithSip(matchList, wcs, order)
        
    distortedWcs = sipObject.getNewWcs()            
    
else:
    print "Warning, no solution found"

#Finally, if you wish to solve another field, first reset the solver
gas.reset()

#Now you can set a new field and re-solve
