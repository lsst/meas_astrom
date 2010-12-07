import os
import sys
from optparse import OptionParser

#import lsst.ip.pipeline as ipPipe
from lsst.datarel import runStage
import lsst.meas.pipeline as measPipe
from stageCtrl import *

from lsst.obs.lsstSim import LsstSimMapper
import lsst.daf.persistence as dafPersist

import imsimUtils

def process(keys, inButler, outButler):
    # HACK!!!
    # visitim = inButler.get('visitim', **keys)
    print 'Processing', keys
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

def main():
    parser = OptionParser()
    imsimUtils.addOptions(parser, input=True, output=True)
    parser.add_option('-T', '--threads', dest='threads', default=None, help='run N processes at once', type='int')
    (opt, args) = parser.parse_args()

    inButler  = imsimUtils.getInputButler(opt)
    outButler = imsimUtils.getOutputButler(opt)

    allkeys = imsimUtils.getAllKeys(opt, inButler)

    if opt.threads is None:
        for keys in allkeys:
            process(keys, inButler, outButler)
    else:
        import multiprocessing
        p = multiprocessing.Pool(opt.threads)
        p.map(process, [(k, inButler, outButler) for k in allkeys])

if __name__ == '__main__':
    main()

