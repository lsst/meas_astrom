from lsst.obs.lsstSim import LsstSimMapper
import lsst.daf.persistence as dafPersist


def addOptions(parser, input=True, output=False):
    if input:
        parser.add_option('-i', '--input', dest='inRoot', default='.', help='input root')
    if output:
        parser.add_option('-o', '--output', dest='outRoot', default='.', help='output root')
    parser.add_option('-R', '--registry', help='registry', dest='registry')
    parser.add_option('-v', '--visit', dest='visit', action='append', type='int', default=[])
    parser.add_option('-r', '--raft', dest='raft', action='append', default=[])
    parser.add_option('-s', '--sensor', dest='sensor', action='append', default=[])

def getInputButler(opt):
    inmapper = LsstSimMapper(root=opt.inRoot, registry=opt.registry)
    bf = dafPersist.ButlerFactory(mapper=inmapper)
    inButler = bf.create()
    return inButler

def getOutputButler(opt):
    outmapper = LsstSimMapper(root=opt.outRoot, registry=opt.registry)
    bf = dafPersist.ButlerFactory(mapper=outmapper)
    outButler = bf.create()
    return outButler

def getAllKeys(opt, inButler):
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
    return allkeys
