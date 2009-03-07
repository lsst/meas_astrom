
#example.py
#An example of how to use the GlobalAstrometrySolution class in Python
#Fergal Mullally

#Usage: python example.py


import re
import os
import glob
import math
from datetime import datetime

import eups
import lsst.afw.image as afwImage
import lsst.meas.astrom.net as net
import lsst.afw.detection.detectionLib as detect


def loadXYFromFile(filename):
    """Load a list of positions from a file and store in a Source"""
    f= open(filename)
    
    s1=detect.SourceContainer()
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
        source.setPsfMag(flux)

        s1.append(source)
        
        i=i + 1
    f.close()
    
    return s1


#Create a new Global Solution object
gas = net.GlobalAstrometrySolution()

#Read in a starlist from a file, and set it as the image to solve
#The sourcelist is a Source object. See loadXYFromFile() above.
path = os.path.join(eups.productDir("meas_astrom"), "tests", "cfht.xy.txt")
starlist = loadXYFromFile(path)
gas.setStarlist(starlist)

#An alternative method is to create the object with starlist as the argument
#gas2=net.GlobalAstrometrySolution(starlist)

#Read in the index files. These files list the various star patterns that
#will be matched against to determine where in the sky your image is.
#These files are big, and this step takes about 20-30 seconds
path= os.path.join(eups.productDir("astrometry_net_data"), "index-*.fits")

indices = glob.glob(path)
for f in indices:
    print "Loading index from %s" % (f)
    gas.addIndexFile(f)

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


if success:
    print "Solution found"
    
    #Find, for example, the pixel position of the intial guess ra/dec
    #In general, you will want to return a Wcs object and interrogate
    #it, rather than the gas object itself, but this function is handy
    #for debugging.
    xy = gas.raDecToXY(ra, dec)
    print xy.getX(), xy.getY()
    
    #Two types of Wcs object can be returned, with or without distortion.
    #Distortion corrections are implemented via SIP polynomials (see
    #Shupe et al. 2005, Astron. Data Anal. Software & Systems XIV,
    #ASP Conf. Series Vol XXX). You will usually want the distortion 
    #corrected Wcs.
    wcs = gas.getDistortedWcs()
    type(wcs)
    
    #Note that the centre of the wcs solution is not necesarily the centre
    #of the image
    centre = wcs.getOriginXY()
    type(centre)
    print centre.getX(), centre.getY()
else:
    print "Warning, no solution found"

#Finally, if you wish to solve another field, first reset the solver
gas.reset()

#Now you can set another sourcelist to solve and repeat the process.
#The index files need to be loaded only once for each gas object, so 
#the overhead cost can be spread over many images
