#!/usr/bin/env python

import os, sys
from math import *
import eups
import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as logging
import lsst.pex.policy as policy
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils
import lsst
import lsst.meas.astrom as measAstrom
from lsst.pex.logging import Log
from lsst.afw.detection import Source, PersistableSourceVector

import lsst.meas.algorithms

if __name__ == '__main__':
    matchfilename = 'examples/imsim-v85501867-fi-R11-S11-matches.fits'
    sourcefilename = 'examples/imsim-v85501867-fi-R11-S11-sources.boost'

    pol = policy.Policy()
    pol.set('matchThreshold', 30)

    persistence = dafPersist.Persistence.getPersistence(pol)
    additionalData = dafBase.PropertySet()

    # match list
    loc = dafPersist.LogicalLocation(matchfilename)
    storageList2 = dafPersist.StorageList()
    storage2 = persistence.getRetrieveStorage('FitsStorage', loc)
    storageList2.append(storage2)
    matchlistptr = persistence.unsafeRetrieve("PersistableSourceMatchVector", storageList2, additionalData)
    psmv2 = afwDet.PersistableSourceMatchVector.swigConvert(matchlistptr)
    smv2 = psmv2.getSourceMatches()
    extra2 = psmv2.getSourceMatchMetadata()

    # sources
    loc3 = dafPersist.LogicalLocation(sourcefilename)
    storageList3 = dafPersist.StorageList()
    storage3 = persistence.getRetrieveStorage('BoostStorage', loc3)
    storageList3.append(storage3)
    psvptr = persistence.unsafeRetrieve('PersistableSourceVector', storageList3, additionalData)
    #psvptr = persistence.unsafeRetrieve('lsst.afw.detection.PersistableSourceVector', storageList3, additionalData)
    psv = afwDet.PersistableSourceVector.swigConvert(psvptr)
    sources = psv.getSources()

    print 'Got sources:', sources

    print 'Got SMV:', smv2
    print 'Got metadata:', extra2

    filtername = 'i'
    idname = 'id'

    measAstrom.joinMatchListWithCatalog(smv2, extra2, pol, filterName=filtername, idName=idname)

    # join sources.
    idtosrc = dict([[s.getSourceId(), s] for s in sources])
    nmatched = 0
    for i in xrange(len(smv2)):
        s = smv2[i].second
        sid = s.getSourceId()
        if sid in idtosrc:
            smv2[i].second = idtosrc[sid]
            nmatched += 1
    print 'Joined', nmatched, 'sources'

    print 'Got', len(sources), 'sources'
    print 'ids:', [s.getSourceId() for s in sources]

    print 'Matchlist source ids:', [m.second.getSourceId() for m in smv2]

    
    if False:
        for i in xrange(len(smv2)):
            r = smv2[i].first
            s = smv2[i].second
            print r.getSourceId(), r.getRa(), r.getDec(), r.getPsfFlux(), s.getXAstrom(), s.getYAstrom()
        

