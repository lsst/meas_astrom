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

import os
from math import hypot

import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Log, Debug, LogRec, Prop
from lsst.pex.exceptions import LsstCppException
import lsst.afw.image as afwImg
import lsst.daf.base as dafBase
import lsst.afw.coord as afwCoord

import net as astromNet
import sip as astromSip
import sip.cleanBadPoints as cleanBadPoints

import lsst.afw.display.ds9 as ds9
import numpy

try:
    import lsstDebug

    display = lsstDebug.Info(__name__).display
except ImportError, e:
    try:
        type(display)
    except NameError:
        display = False

def createSolver(policy, log):
    path=os.path.join(os.environ['ASTROMETRY_NET_DATA_DIR'], "metadata.paf")
    solver = astromNet.GlobalAstrometrySolution(path, log)
    matchThreshold = policy.get('matchThreshold')
    solver.setMatchThreshold(matchThreshold)
    # FIXME -- this could go in policy... or we could use new Astrometry.net logging
    # callbacks to put those messages to their own pexLogging channel.
    solver.setLogLevel(2)
    return solver

def getIdColumn(policy):
    '''Returns the column name of the ID field in the reference catalog'''
    idName = ''
    colname = 'defaultIdColumnName'
    if policy.exists(colname):
        idName = policy.get(colname)
    return idName

def joinMatchListWithCatalog(matchlist, matchmeta, policy, log=None, solver=None,
                             filterName=None, idName=None):
    if log is None:
        log = Log.getDefaultLog()

    if solver is None:
        solver = createSolver(policy, log)

    filterName = chooseFilterName(None, policy, solver, log, filterName)
    if idName is None:
        idName = getIdColumn(policy)

    version = matchmeta.getInt('SMATCHV')
    if version != 1:
        raise ValueError('SourceMatchVector version number is %i, not 1.' % version)

    # EUPS
    myandata = os.environ.get('ASTROMETRY_NET_DATA_DIR')
    andata = matchmeta.getString('ANEUPS')
    #if os.path.basename(myandata) != os.path.basename(andata):
    #    raise ValueError('Need ASTROMETRY_NET_DATA_DIR = "%s"' % 
    log.log(Log.DEBUG, 'Astrometry.net dir was "%s", now "%s"' %
            (os.path.basename(myandata), os.path.basename(andata)))

    anid = matchmeta.getInt('ANINDID')
    anhp = matchmeta.getInt('ANINDHP')
    anindexname = os.path.basename(matchmeta.getString('ANINDNM'))
    log.log(Log.DEBUG, 'Astrometry.net index was "%s" (id %i, healpix %i)' %
            (anindexname, anid, anhp))

    # all in deg.
    ra = matchmeta.getDouble('RA')
    dec = matchmeta.getDouble('DEC')
    rad = matchmeta.getDouble('RADIUS')
    log.log(Log.DEBUG, 'Searching RA,Dec %.3f,%.3f, radius %.1f arcsec, filter "%s", id column "%s", indexid %i' %
            (ra, dec, rad * 3600., filterName, idName, anid))
    #myinds = solver.getIndexList()

    # FIXME -- need anid?  Not necessarily... ref ids are supposed to be unique!
    cat = solver.getCatalogue(ra, dec, rad * 3600., filterName, idName, anid)
    log.log(Log.DEBUG, 'Found %i reference catalog sources in range' % len(cat))

    # build map of reference id to reference object.
    idtoref = dict([[s.getSourceId(), s] for s in cat])

    # Join.
    nmatched = 0
    for i in xrange(len(matchlist)):
        mid = matchlist[i].first.getSourceId()
        if mid in idtoref:
            ref = idtoref[mid]
            matchlist[i].first = ref
            nmatched += 1
    log.log(Log.DEBUG, 'Joined %i of %i matchlist reference IDs to reference objects' %
            (nmatched, len(matchlist)))

# Object returned by determineWcs.
class InitialAstrometry(object):
    self __init__(self):
        self.matches = None
        self.wcs = None
    def getMatches(self):
        return self.matches
    def getWcs(self):
        return self.wcs
    def getMatchMetadata(self):
        return getattr(self, 'matchMetadata', None)

def determineWcs(policy, exposure, sourceSet, log=None, solver=None, doTrim=False,
                 forceImageSize=None, filterName=None):
    '''Top level function for calculating an initial (per-chip) astrometric solution.

    Get an initial World Coordinate System (WCS) from Astrometry.net,
    then calculate SIP distortion terms.


    Input:
    policy:     An lsst.pex.policy.Policy object containing the parameters for the solver
    exposure    lsst.afw.image.Exposure representation of an image and a WCS 
                this provides the initial guess at position and plate scale
    sourceSet   A list of lsst.afw.detection.Source objects, indicating the pixel positions of
                stars in the field
    log         A lsst.pex.logging.Log object (optional), used for printing progress
    doTrim      Remove sources that are not inside the image.
    solver      Optionally provide a previously created astrometry.net solver. If not provided
                one will be created.
    forceImageSize  tuple of (W,H): force this image size, rather than getting it from the Exposure.
    filterName  Use this filter name, rather than getting it from the exposure.
    '''

    astrom = InitialAstrometry()

    if log is None:
        log = Log.getDefaultLog()

    if display:
        frame = 1
        ds9.mtv(exposure, frame=frame, title="wcsDet")

    if doTrim:
        nStart = len(sourceSet)
        sourceSet = trimBadPoints(exposure, sourceSet)
        if log:
            nEnd = len(sourceSet)
            log.log(log.DEBUG, "Kept %i of %i sources after trimming" %(nEnd, nStart))

    if display:
        for s in sourceSet:
            ds9.dot("o", s.getXAstrom(), s.getYAstrom(), size=3, ctype=ds9.RED, frame=frame)

    #Extract an initial guess WCS if available    
    wcsIn = exposure.getWcs() #May be None
    # Exposure uses the special object "NoWcs" instead of NULL.  Because they're special.
    haswcs = exposure.hasWcs()
    if not haswcs:
        log.log(log.WARN, "No WCS found in exposure. Doing blind solve")

    # Setup solver
    if solver is None:
        solver = createSolver(policy, log)
    else:
        solver.reset()

    # Set solving params
    log.log(log.DEBUG, "Setting starlist")
    solver.setStarlist(sourceSet)
    log.log(log.DEBUG, "Setting numBrightObj")
    solver.setNumBrightObjects( min(policy.get('numBrightStars'), len(sourceSet)))
    if forceImageSize is not None:
        (W,H) = forceImageSize
    else:
        (W,H) = (exposure.getWidth(), exposure.getHeight())
    solver.setImageSize(W, H)
    #solver.printSolverSettings(stdout)

    key = 'pixelScaleUncertainty'
    if policy.exists(key):
        dscale = float(policy.get(key))
    else:
        dscale = None

    # Do a blind solve if we're told to, or if we don't have an input WCS
    doBlindSolve = policy.get('blindSolve') or (not haswcs)
    if doBlindSolve:
        log.log(log.DEBUG, "Solving with no initial guess at position")
        isSolved = solver.solve()
    elif dscale is not None:
        isSolved = solver.solve(wcsIn, dscale)
    else:
        isSolved = solver.solve(wcsIn)

    moreMeta = dafBase.PropertyList()

    # Did we solve?
    log.log(log.DEBUG, "Finished Solve step.")
    if not isSolved:
        log.log(log.WARN, "No solution found, using input WCS")
        return astrom
    wcs = solver.getWcs()

    # Generate a list of catalogue objects in the field.
    imgSizeInArcsec = wcs.pixelScale() * hypot(W,H)
    filterName = chooseFilterName(exposure, policy, solver, log, filterName)
    idName = getIdColumn(policy)
    try:
        #print 'size, fiter, id', imgSizeInArcsec, filterName, idName
        #cat = solver.getCatalogue(2*imgSizeInArcsec, filterName, idName)
        margin = 50 # pixels
        cat = solver.getCatalogueForSolvedField(filterName, idName, margin)
    except LsstCppException, e:
        log.log(Log.WARN, str(e))
        log.log(Log.WARN, "Attempting to access catalogue positions and fluxes")
        version = os.environ['ASTROMETRY_NET_DATA_DIR']
        log.log(Log.WARN, "Catalogue version: %s" %(version))
        log.log(Log.WARN, "ID column: %s" %(idName))
        log.log(Log.WARN, "Requested filter: %s" %(filterName))
        log.log(Log.WARN, "Available filters: " + str(solver.getCatalogueMetadataFields()))
        raise

    if True:
        # Now generate a list of matching objects
        distInArcsec = policy.get('distanceForCatalogueMatchinArcsec')
        cleanParam = policy.get('cleaningParameter')

        matchList = matchSrcAndCatalogue(cat=cat, img=sourceSet, wcs=wcs,
            distInArcsec=distInArcsec, cleanParam=cleanParam)

        uniq = set([sm.second.getId() for sm in matchList])
        if len(matchList) != len(uniq):
            log.log(Log.WARN, "The list of matches stars contains duplicated reference sources (%i sources, %i unique ids)"
                    % (len(matchList), len(uniq)))

        if len(matchList) == 0:
            log.log(Log.WARN, "No matches found between input source and catalogue.")
            log.log(Log.WARN, "Something is wrong. Defaulting to input WCS")
            return astrom

        log.log(Log.DEBUG, "%i catalogue objects match input source list using linear WCS" %(len(matchList)))
    else:
        # Use list of matches returned by Astrometry.net
        log.log(Log.DEBUG, "Getting matched sources: Fluxes in column %s; Ids in column" % (filterName, idName))
        matchList = solver.getMatchedSources(filterName, idName)

    astrom.tanWcs = wcs
    astrom.tanMatches = matchList

    srcids = [s.getSourceId() for s in sourceSet]
    #print 'srcids:', srcids
    for m in matchList:
        #print 'Matchlist entry ids:', m.first.getSourceId(), m.second.getSourceId()
        assert(m.second.getSourceId() in srcids)
        assert(m.second in sourceSet)

    if policy.get('calculateSip'):
        sipOrder = policy.get('sipOrder')
        wcs, matchList = calculateSipTerms(wcs, cat, sourceSet, distInArcsec, cleanParam, sipOrder, log)

        astrom.sipWcs = wcs
        astrom.sipMatches = matchList
    else:
        log.log(Log.DEBUG, "Updating WCS in input exposure with linear WCS")

    exposure.setWcs(wcs)

    # add current EUPS astrometry_net_data setup.
    andata = os.environ.get('ASTROMETRY_NET_DATA_DIR')
    if andata is None:
        moreMeta.add('ANEUPS', 'none', 'ASTROMETRY_NET_DATA_DIR')
    else:
        andata = os.path.basename(andata)
        moreMeta.add('ANEUPS', andata, 'ASTROMETRY_NET_DATA_DIR')

    # cache: field center and size.  These may be off by 1/2 or 1 or 3/2 pixels.
    # dstn does not care.
    cx,cy = W/2.,H/2.
    radec = wcs.pixelToSky(cx, cy)
    ra,dec = radec.getLongitude(afwCoord.DEGREES), radec.getLatitude(afwCoord.DEGREES)
    moreMeta.add('RA', ra, 'field center in degrees')
    moreMeta.add('DEC', dec, 'field center in degrees')
    moreMeta.add('RADIUS', imgSizeInArcsec/2./3600.,
            'field radius in degrees, approximate')
    moreMeta.add('SMATCHV', 1, 'SourceMatchVector version number')

    if display:
        for s1, s2, d in matchList:
            # plot the catalogue positions
            ds9.dot("+", s1.getXAstrom(), s1.getYAstrom(), size=3, ctype=ds9.BLUE, frame=frame)

    matchListMeta = solver.getMatchedIndexMetadata()
    moreMeta.combine(matchListMeta)

    astrom.matchMetadata = moreMeta
    astrom.wcs = wcs
    astrom.matches = matchList

    return astrom

def trimBadPoints(exposure, sourceSet):
    """Remove elements from sourceSet whose xy positions aren't within the boundaries of exposure

    Input:
    exposure:    an Exposure object
    sourceSet  A list of Source objects
    """

    x0, y0 = exposure.getMaskedImage().getXY0()
    h, w = float(exposure.getHeight()), float(exposure.getWidth())

    goodSet = []
    for s in sourceSet:
        if x0 < s.getXAstrom() < x0+w:
            if y0 < s.getYAstrom() < y0+h:
                goodSet.append(s)

    return goodSet


def chooseFilterName(exposure, policy, solver, log, filterName=None):
    """When extracting catalogue magnitudes, which colour filter should we request
    e.g U,B,V etc."""

    if log is None:
        log = Log.getDefaultLog()

    if filterName is None:
        filterName = exposure.getFilter().getName()

    if filterName == "_unknown_":
        log.log(log.DEBUG, "Exposure has no filter name set. Using default.")
    else:
        log.log(Log.DEBUG, 'Exposure was taken with filter "%s"' % (filterName))

    availableFilters = solver.getCatalogueMetadataFields()
    if filterName in availableFilters:
        log.log(Log.DEBUG, 'Have catalogue magnitudes for filter: "%s"' %(filterName))
        return filterName

    log.log(Log.DEBUG, 'Catalogue doesn\'t contain filter "%s"; using default (available filters: "%s")' %
            (filterName, '", "'.join(availableFilters)))

    if not policy.exists("defaultFilterName"):
        log.log(log.DEBUG, "No default filter name is set")
        return ""
    defaultFilter = policy.get("defaultFilterName")
    if defaultFilter in availableFilters:
        log.log(log.DEBUG, 'Using default filter name "%s"' % (defaultFilter))
        return defaultFilter

    raise ValueError('Default filter "%s" not included in catalogue (available filters: "%s")' \
                     % (defaultFilter, '", "'.join(availableFilters)))


def calculateSipTerms(inputWcs, cat, sourceSet, distInArcsec, cleanParam, sipOrder, log=None):
    """Iteratively calculate sip distortions and regenerate matchList based on improved wcs"""

    if log is None:
        log = Log.getDefaultLog()

    wcs = inputWcs

    #Create a first pass at a set of matching objects
    matchList = matchSrcAndCatalogue(cat=cat, img=sourceSet, wcs=wcs,
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
        proposedMatchlist = matchSrcAndCatalogue(cat=cat, img=sourceSet, wcs=proposedWcs,
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


