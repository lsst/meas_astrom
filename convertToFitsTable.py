#!/usr/bin/env python
"""Script to convert meas/astrom test data files to the new table paradigm.

The original files are text files with id, x, y, ra, dec, psfFlux, flags columns.

The new files are FITS binary tables that can be loaded into afw.table.SourceCatalog
objects, and contain roughly equivalent fields that mimic those produced by common
measurement algorithms in meas/alg.
"""
from __future__ import print_function

import lsst.afw.table

oldFlags = {
    "EDGE": 0x000001,
    "SHAPE_SHIFT": 0x000002,
    "SHAPE_MAXITER": 0x000004,
    "SHAPE_UNWEIGHTED": 0x000008,
    "SHAPE_UNWEIGHTED_PSF": 0x000010,
    "SHAPE_UNWEIGHTED_BAD": 0x000020,
    "PEAKCENTER": 0x000040,
    "BINNED1": 0x000080,
    "INTERP": 0x000100,
    "INTERP_CENTER": 0x000200,
    "SATUR": 0x000400,
    "SATUR_CENTER": 0x000800,
    "DETECT_NEGATIVE": 0x001000,
    "STAR": 0x002000,
    "PSFSTAR": 0x004000,
    "PHOTOM_NO_PSF": 0x008000,
    "PHOTOM_NO_PEAK": 0x010000,
    "PHOTOM_NO_SOURCE": 0x020000,
    "PHOTOM_NO_FOOTPRINT": 0x040000,
}

flagMapping = {
    "EDGE": "flags.pixel.edge",
    "SHAPE_SHIFT": "shape.sdss.flags.shift",
    "SHAPE_MAXITER": "shape.sdss.flags.maxiter",
    "SHAPE_UNWEIGHTED": "shape.sdss.flags.unweighted",
    "SHAPE_UNWEIGHTED_BAD": "shape.sdss.flags.unweightedbad",
    "PEAKCENTER": "flags.badcentroid",
    "INTERP": "flags.pixel.interpolated.any",
    "INTERP_CENTER": "flags.pixel.interpolated.center",
    "SATUR": "flags.pixel.saturated.any",
    "SATUR_CENTER": "flags.pixel.saturated.center",
}


def main(inputName, outputName):
    schema = lsst.afw.table.SourceTable.makeMinimalSchema()
    centroidKey = schema.addField("centroid", type=lsst.afw.geom.Point2D)
    schema.addField("centroid.flags", type="Flag")
    schema.addField("centroid.err", type="CovPointF")
    psfFluxKey = schema.addField("flux.psf", type=float)
    schema.addField("flux.psf.flags", type="Flag")
    schema.addField("flux.psf.err", type=float)
    for name in ("edge", "interpolated.any", "interpolated.center", "saturated.any", "saturated.center"):
        schema.addField("flags.pixel." + name, type="Flag")
    for name in ("unweightedbad", "unweighted", "shift", "maxiter"):
        schema.addField("shape.sdss.flags." + name, type="Flag")
    schema.addField("flags.badcentroid", type="Flag")
    outputCat = lsst.afw.table.SourceCatalog(schema)
    outputCat.table.defineCentroid("centroid")
    outputCat.table.definePsfFlux("flux.psf")
    with open(inputName, 'r') as inputFile:
        lineno = 0
        for line in inputFile:
            lineno += 1
            try:
                id, x, y, ra, dec, cts, flags = line.split()
            except Exception as e:
                print("Line %d: %s: %s" % (lineno, e, line), end=' ')
            record = outputCat.addNew()
            record.setId(int(id))
            record.setRa(float(ra) * lsst.afw.geom.degrees)
            record.setDec(float(dec) * lsst.afw.geom.degrees)
            record.set(centroidKey.getX(), float(x))
            record.set(centroidKey.getY(), float(y))
            record.set(psfFluxKey, float(cts))
            flags = int(flags)
            for oldName, mask in oldFlags.iteritems():
                if oldName == "BINNED1":
                    continue
                if mask & flags:
                    newName = flagMapping.get(oldName, None)
                    if newName is None:
                        print("Ignoring flag bit '%s'" % oldName)
                    else:
                        record.set(newName, True)
    outputCat.writeFits(outputName)

if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
