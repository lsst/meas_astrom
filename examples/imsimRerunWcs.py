import os
import sys
from optparse import OptionParser

#import lsst.ip.pipeline as ipPipe
from lsst.datarel import runStage
import lsst.meas.pipeline as measPipe
from stageCtrl import *

from lsst.obs.lsstSim import LsstSimMapper
import lsst.daf.persistence as dafPersist

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-i', '--input', dest='inRoot', default='.', help='input root')
    parser.add_option('-o', '--output', dest='outRoot', default='.', help='output root')
    parser.add_option('-R', '--registry', help='registry', dest='registry')
    parser.add_option('-v', '--visit', type='int', dest='visit', action='append', type='int', default=[])
    parser.add_option('-r', '--raft', dest='raft', action='append', default=[])
    parser.add_option('-s', '--sensor', dest='sensor', action='append', default=[])
    (opt, args) = parser.parse_args()

    inmapper = LsstSimMapper(root=opt.inRoot, registry=opt.registry)
    bf = dafPersist.ButlerFactory(mapper=inmapper)
    inButler = bf.create()

    outmapper = LsstSimMapper(root=opt.outRoot, registry=opt.registry)
    bf = dafPersist.ButlerFactory(mapper=outmapper)
    outButler = bf.create()

    allkeys = []

    if not len(opt.visit):
        # Grab all available visits.
        print 'Grabbing all available visits...'
        #visits = inButler.queryMetadata('visitim', 'visit')
        visits = inButler.queryMetadata('raw', 'visit')
        print 'Got visits', visits
        opt.visit = visits

    for visit in opt.visit:
        if not len(opt.raft):
            print 'Grabbing all available rafts for visit', visit
            rafts = inButler.queryMetadata('raw', 'raft', visit=visit)
            print 'Got rafts', rafts
        else:
            rafts = opt.raft

        for raft in rafts:
            if not len(opt.sensor):
                print 'Grabbing all available sensors for visit', visit, 'raft', raft
                sensors = inButler.queryMetadata('raw', 'sensor', visit=visit, raft=raft)
                print 'Got sensors', sensors
            else:
                sensors = opt.sensor
            for sensor in sensors:
                allkeys.append(dict(visit=visit, raft=raft, sensor=sensor))

    for keys in allkeys:
        print 'Processing', keys

        #visitim = inButler.get('visitim', **keys)
        # HACK!!!
        visitim = inButler.get('calexp', **keys)
        sourceset_p = inButler.get('icSrc', **keys)
        sourceset = sourceset_p.getSources()

        clip = {
            'visitExposure': visitim,
            'sourceSet': sourceset,
            }

        clip = runStage(measPipe.WcsDeterminationStage,
                        """#<?cfg paf policy?>
                        inputExposureKey: visitExposure
                        inputSourceSetKey: sourceSet
                        outputWcsKey: measuredWcs
                        outputMatchListKey: matchList
                        numBrightStars: 150
                        defaultFilterName: mag
                        """, clip)

        clip = runStage(measPipe.PhotoCalStage,
                        """#<?cfg paf policy?>
                        sourceMatchSetKey: matchList
                        outputValueKey: photometricMagnitudeObject
                        """, clip)

        outButler.put(clip['matchList_persistable'], 'icMatch', **keys)
        outButler.put(clip['visitExposure'], 'calexp', **keys)

