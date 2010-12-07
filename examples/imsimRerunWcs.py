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

if __name__ == '__main__':
    parser = OptionParser()
    imsimUtils.addOptions(parser, input=True, output=True)
    (opt, args) = parser.parse_args()

    inButler  = imsimUtils.getInputButler(opt)
    outButler = imsimUtils.getOutputButler(opt)

    allkeys = imsimUtils.getAllKeys(opt, inButler)

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

