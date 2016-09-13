import sys

import lsst.pex.policy as policy
import lsst.meas.astrom as measAstrom
from lsst.pex.logging import Log
import lsst.afw.geom as afwGeom


def main():
    from optparse import OptionParser

    parser = OptionParser(usage='%(program) [args] RA Dec radius')
    parser.add_option('-o', dest='outfn', help='FITS table output filename', default=None)
    (opt, args) = parser.parse_args()

    if len(args) != 3:
        parser.print_help()
        return -1

    ra = float(args[0])
    dec = float(args[1])
    radius = float(args[2])

    log = Log.getDefaultLog()
    log.setThreshold(Log.DEBUG)
    pol = policy.Policy()
    pol.set('matchThreshold', 30)
    solver = measAstrom.createSolver(pol, log)

    solver.setLogLevel(3)

    ids = solver.getIndexIdList()
    print 'Index IDs:', ids
    indexid = ids[0]

    idName = 'id'
    X = solver.getCatalogue(ra * afwGeom.degrees, dec * afwGeom.degrees,
                            radius * afwGeom.degrees, '', idName, indexid)
    ref = X.refsources
    inds = X.inds
    print 'Got', len(ref), 'reference catalog sources'
    print '  got indices:', len(inds)

    print 'Tag-along columns:'
    cols = solver.getTagAlongColumns(indexid)
    # print cols
    for c in cols:
        print '  column:', c.name, c.fitstype, c.ctype, c.units, c.arraysize
    colnames = [c.name for c in cols]

    tagdata = []
    for c in cols:
        fname = 'getTagAlong' + c.ctype
        func = getattr(solver, fname)
        data = func(indexid, c.name, inds)
        # print 'called', fname, 'to get', c.name, c.ctype, '(len %i)' % len(data)
        tagdata.append(data)

    if opt.outfn is None:
        # SSV
        print 'ra dec',
        for c in cols:
            if c.arraysize > 1:
                for a in len(c.arraysize):
                    print('%s_%i' % (c.name, a)),
            else:
                print c.name,
        print

        for i, r in enumerate(ref):
            print r.getRa().asDegrees(), r.getDec().asDegrees(),
            for c, d in zip(cols, tagdata):
                if c.arraysize > 1:
                    for a in len(c.arraysize):
                        print d[c.arraysize * i + a],
                else:
                    print d[i],
            print

    else:
        import pyfits
        import numpy as np

        fitscols = []
        fitscols.append(pyfits.Column(name='RA', array=np.array([r.getRa().asDegrees() for r in ref]),
                                      format='D', unit='deg'))
        fitscols.append(pyfits.Column(name='DEC', array=np.array([r.getDec().asDegrees() for r in ref]),
                                      format='D', unit='deg'))
        for c, d in zip(cols, tagdata):
            fmap = {'Int64': 'K',
                    'Int': 'J',
                    'Bool': 'L',
                    'Double': 'D',
                    }
            if c.arraysize > 1:
                # May have to reshape the array as well...
                fitscols.append(pyfits.Column(name=c.name, array=np.array(d),
                                              format='%i%s' % (c.arraysize, fmap.get(c.ctype, 'D'))))
            else:
                fitscols.append(pyfits.Column(name=c.name, array=np.array(d),
                                              format=fmap.get(c.ctype, 'D')))

        pyfits.new_table(fitscols).writeto(opt.outfn, clobber=True)
        print 'Wrote FITS table', opt.outfn

    return 0


if __name__ == '__main__':
    sys.exit(main())
