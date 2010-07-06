import os
import math

import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Log, Debug, LogRec, Prop
from lsst.pex.exceptions import LsstCppException
import lsst.afw.image as afwImg

import net as astromNet
import sip as astromSip
import sip.cleanBadPoints as cleanBadPoints

import lsst.afw.display.ds9 as ds9

try:
    import lsstDebug

    display = lsstDebug.Info(__name__).display
except ImportError, e:
    try:
        type(display)
    except NameError:
        display = False

def determineWcs(policy, exposure, sourceSet, log=None, solver=None, doTrim=False, forceImageSize=None):
    """Top level function for calculating a Wcs. 
    
    Given an initial guess at a Wcs (hidden inside an exposure) and a set of
    sources (sourceSet), use astrometry.net to confirm the Wcs, then calculate
    distortion terms.
    
    Input:
    policy:     An lsst.pex.policy.Policy object containing the parameters for the solver
    exposure    lsst.afw.image.Exposure representation of an image and a wcs 
                this provides the initial guess at position and plate scale
    sourceSet   A list of lsst.afw.detection.Source objects, indicating the pixel positions of 
                stars in the field
    log         A lsst.pex.logging.Log object (optional), used for printing progress
    doTrim        Check that all sources lie within the image, and remove those that don't.
    solver      Optionally provide a previously created astrometry.net solver. If not provided
                one will be created.
    forceImageSize  tuple of (W,H): force this image size, rather than getting it from the Exposure.
    """
	sipwcs = None
	tanMatches = None
	sipMatches = None

    if log is None:
        log = StdoutLog()   #Write log messages to stdout
    log.log(Log.INFO, "In determineWcs")


    #Short names
    exp = exposure
    srcSet = sourceSet

    if display:
        frame = 1
        ds9.mtv(exposure, frame=frame, title="wcsDet")

    if doTrim:
        nStart = len(srcSet)
        srcSet = trimBadPoints(exp, srcSet)
        if log:
            nEnd = len(srcSet)
            log.log(log.DEBUG, "Kept %i of %i sources after trimming" %(nEnd, nStart))

    if display:
        for s in srcSet:
            ds9.dot("o", s.getXAstrom(), s.getYAstrom(), size=3, ctype=ds9.RED, frame=frame)
        
    #Extract an initial guess wcs if available    
    wcsIn = exp.getWcs() #May be None
    if wcsIn is None:
        log.log(log.WARN, "No wcs found on exposure. Doing blind solve")
    
    #Setup solver
    if solver is None:
        path=os.path.join(os.environ['ASTROMETRY_NET_DATA_DIR'], "metadata.paf")
        solver = astromNet.GlobalAstrometrySolution(path)
        matchThreshold = policy.get('matchThreshold')
        solver.setMatchThreshold(matchThreshold)
    else:
        solver.reset()

    #Set solving params
    log.log(log.DEBUG, "Setting starlist")
    solver.setStarlist(srcSet)
    log.log(log.DEBUG, "Setting numBrightObj")
    solver.setNumBrightObjects( min(policy.get('numBrightStars'), len(srcSet)))

    if forceImageSize is not None:
        W,H = forceImageSize
    else:
        W,H = exp.getWidth(), exp.getHeight()

	solver.setImageSize(W, H)
    solver.setLogLevel(2)
    #solver.printSolverSettings(stdout)

    # FIXME -- add policy entry for this...
    #dscale = policy.get('pixelScaleUncertainty')
    dscale = None

    #Do a blind solve if we're told to, or if we don't have an input wcs
    doBlindSolve = policy.get('blindSolve') or (wcsIn is None)
    if doBlindSolve:
        log.log(log.DEBUG, "Solving with no initial guess at position")
        isSolved = solver.solve()
    elif dscale is not None:
        isSolved = solver.solve(wcsIn, dscale)
    else:
        isSolved = solver.solve(wcsIn)

    #Did we solve?
    log.log(log.DEBUG, "Finished Solve step.")
    if not isSolved:
        log.log(log.WARN, "No solution found, using input Wcs")
        return [], wcsIn
    wcs = solver.getWcs()
	tanwcs = wcs

    #
    # Generate a list of catalogue objects in the field.
    #
    
    #First obtain the catalogue-listed positions of stars
    log.log(log.DEBUG, "Determining match objects")
    imgSizeInArcsec = getImageSizeInArcsec(srcSet, wcs)
    
    #Do we want magnitude information
    filterName = chooseFilterName(exposure, policy, solver, log)
    try:
        cat = solver.getCatalogue(2*imgSizeInArcsec, filterName) 
    except LsstCppException, e:
        log.log(Log.WARN, str(e))
        log.log(Log.WARN, "Attempting to access catalogue positions and fluxes")
        version = os.environ['ASTROMETRY_NET_DATA_DIR']
        log.log(Log.WARN, "Catalogue version: %s" %(version))
        log.log(Log.WARN, "Requested filter: %s" %(filterName))
        log.log(Log.WARN, "Available filters: " + str(solver.getCatalogueMetadataFields()))
        raise
            

    matchList=[]    #Make sure this stays in scope
    if True:
        #Now generate a list of matching objects
        distInArcsec = policy.get('distanceForCatalogueMatchinArcsec')
        cleanParam = policy.get('cleaningParameter')

        matchList = matchSrcAndCatalogue(cat=cat, img=srcSet, wcs=wcs, 
            distInArcsec=distInArcsec, cleanParam=cleanParam)
		tanMatches = matchList

        if len(matchList) == 0:
            log.log(Log.WARN, "No matches found between input source and catalogue.")
            log.log(Log.WARN, "Something in wrong. Defaulting to input wcs")
            return [], wcsIn
            
        log.log(Log.DEBUG, "%i catalogue objects match input source list using linear Wcs" %(len(matchList)))
    else:
        #Use list of matches returned by astrometry.net
        log.log(Log.DEBUG, "Getting matched sources: Fluxes in band %s " %(filterName))
        matchList = solver.getMatchedSources(filterName)
    

    if policy.get('calculateSip'):
        sipOrder = policy.get('sipOrder')
        wcs, matchList = calculateSipTerms(wcs, cat, srcSet, distInArcsec, cleanParam, sipOrder, log)
		sipwcs = wcs
		sipMatches = matchList
    else:
        log.log(Log.DEBUG, "Updating wcs in input exposure with linear wcs")
        
    exposure.setWcs(wcs)

	# We want, for diagnostic plots:
	# -TAN WCS
	# -SIP WCS, if computed
	# -all image sources (sourceSet)
	# -image size (W,H)
	# -all reference sources ("cat")
	# -matchList / Astrometry.net matches
	# - (tan / sip)
	# -solver?
	# -MatchObj?

	#  (better grab this before solver.reset() !)

	if True:
		wcsPlots.wcsPlots(tanwcs, sourceSet, cat, tanMatches, W, H, 'tan-')
		if sipwcs and sipMatches:
			wcsPlots.wcsPlots(sipwcs, sourceSet, cat, sipMatches, W, H, 'sip-')


    solver.reset()

    if display:
        for s1, s2, d in matchList:
            # plot the catalogue positions
            ds9.dot("+", s1.getXAstrom(), s1.getYAstrom(), size=3, ctype=ds9.BLUE, frame=frame)

    return [matchList, wcs]

class StdoutLog():
    """If no log is passed, this class just writes the output to stdout, regardless of
    log verbosity"""
    
    def __init__(self):
        self.DEBUG="DEBUG"
        self.INFO="INFO"
        self.WARN="WARN"
        
    def log(self, arg1, arg2):
        print "%s" %(arg2)


def trimBadPoints(exp, srcSet):
    """Remove elements from srcSet whose xy positions aren't within the boundaries of exp
    
    Input:
    exp:    an Exposure object
    srcSet  A list of Source objects
    """
    
    x0, y0 = exp.getMaskedImage().getXY0()
    h, w = float(exp.getHeight()), float(exp.getWidth())

    goodSet = []
    for s in srcSet:
        if x0 < s.getXAstrom() < x0+w:
            if y0 < s.getYAstrom() < y0+h:
                goodSet.append(s)
    
    return goodSet
    
    
def chooseFilterName(exposure, policy, solver, log):
    """When extracting catalogue magnitudes, which colour filter should we request
    e.g U,B,V etc."""
    
    if log is None:
        log = StdoutLog()
    
    filterName = exposure.getFilter().getName()
    if filterName == "_unknown_":   #No symbol for this in afw
        log.log(log.DEBUG, "Exposure has no filter info. Using default")
        if not policy.exists("defaultFilterName"):
            log.log(log.DEBUG, "No default filter is set")
            return ""
        filterName = policy.get("defaultFilterName")

    log.log(Log.DEBUG, "Exposure was taken in %s band" %(filterName))
    
    availableFilters = solver.getCatalogueMetadataFields()
    availableFiltersStr = ", ".join(availableFilters) #Expressed as a string

    
    if filterName in availableFilters:
        log.log(Log.DEBUG, "Have catalogue magnitudes for %s" %(filterName))
        return filterName
    else:
        log.log(Log.DEBUG, "Catalogue doesn't contain %s, only [%s]" %(filterName, availableFiltersStr))
        log.log(Log.DEBUG, "Searching for default filter")
        
        if not policy.exists("defaultFilterName"):
            log.log(log.DEBUG, "No default filter is set")
            return ""
        defaultFilter = policy.get("defaultFilterName")

        if defaultFilter in availableFilters:
            log.log(log.DEBUG, "Using default filter name (%s)" %(defaultFilter))
            return defaultFilter
        else:
            raise ValueError("Default filter %s not included in catalogue [%s]" \
                    %(defaultFilter, availableFiltersStr))
            
    raise RuntimeError("This function should have returned before getting to this point")



    
def getImageSizeInArcsec(srcSet, wcs):
    """ Get the approximate size of the image in arcseconds
    
    Input: 
    srcSet List of detected objects in the image (with pixel positions)
    wcs    Wcs converts pixel positions to ra/dec
    
    """
    xfunc = lambda x: x.getXAstrom()
    yfunc = lambda x: x.getYAstrom()
    
    x = map(xfunc, srcSet)
    y = map(yfunc, srcSet)
    
    minx = min(x)
    maxx = max(x)
    miny = min(y)
    maxy = max(y)
    
    llc = wcs.pixelToSky(minx, miny)
    urc = wcs.pixelToSky(maxx, maxy)
    
    deltaRa = urc[0] - llc[0]
    deltaDec = urc[1] - llc[1]
    
    #Approximately right
    dist = math.sqrt(deltaRa**2 + deltaDec**2)
    return 3600*math.degrees(dist)  #arcsec


def calculateSipTerms(inputWcs, cat, srcSet, distInArcsec, cleanParam, sipOrder, log=None):
    """Iteratively calculate sip distortions and regenerate matchList based on improved wcs"""

    if log is None:
        log = StdoutLog()

    wcs = inputWcs

    #Create a first pass at a set of matching objects
    matchList = matchSrcAndCatalogue(cat=cat, img=srcSet, wcs=wcs, 
        distInArcsec=distInArcsec, cleanParam=cleanParam)

    i=0
    while True:
        try:
            sipObject = astromSip.CreateWcsWithSip(matchList, wcs, sipOrder)
            proposedWcs = sipObject.getNewWcs()
        except LsstCppException, e:
            log.log(Log.WARN, "Failed to calculate distortion terms. Error:")
            log.log(Log.WARN, str(e))
            log.log(Log.WARN, "Using best guess wcs")
            break

        matchSize = len(matchList)
        msg="Sip Iteration %i: %i objects match. rms scatter is %g arcsec or %g pixels" \
                %(i, matchSize, sipObject.getScatterInArcsec(), sipObject.getScatterInPixels())
        log.log(Log.DEBUG, msg)

        #Use the new wcs to update the match list        
        proposedMatchlist = matchSrcAndCatalogue(cat=cat, img=srcSet, wcs=proposedWcs, 
            distInArcsec=distInArcsec, cleanParam=cleanParam)


        if len(proposedMatchlist) <= matchSize:
            #We're regressing, so stop
            break

        wcs = proposedWcs
        matchList = proposedMatchlist 
        matchSize = len(matchList)
        i=i+1

    if not wcs.hasDistortion():
        log.log(Log.WARN, "Distortion fitter failed to improve on linear WCS")
        
    return wcs, matchList
            

def matchSrcAndCatalogue(cat=None, img=None, wcs=None, distInArcsec=1.0, cleanParam=3):
    """Given an input catalogue, match a list of objects in an image, given
    their x,y position and a wcs solution.
    
    Return: A list of x, y, dx and dy. Each element of the list is itself a list
    """
    
    if cat is None:
        raise RuntimeError("Catalogue list is not set")
    if img is None:
        raise RuntimeError("Image list is not set")
    if wcs is None:
        raise RuntimeError("wcs is not set")
    
        
    matcher = astromSip.MatchSrcToCatalogue(cat, img, wcs, distInArcsec)    
    matchList = matcher.getMatches()
    

    if matchList is None:
        raise RuntimeError("No matches found between image and catalogue")

    matchList = cleanBadPoints.clean(matchList, wcs, nsigma=cleanParam)
    return matchList


