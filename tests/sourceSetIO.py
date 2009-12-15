import math, os, sys
import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage

def write(sourceSet, outfile="-"):
    if outfile == "-":
        fd = sys.stdout
    else:
        fd = open(outfile, "w")

    for s in sourceSet:
        print >> fd, s.getId(), s.getXAstrom(), s.getYAstrom(), s.getRa(), s.getDec(), s.getPsfFlux(), s.getFlagForDetection()



def read(fileName):
    fd = open(fileName, "r")

    sourceSet = afwDetection.SourceSet()
    lineno = 0
    for line in fd.readlines():
        lineno += 1
        try:
            id, x, y, ra, dec, cts, flags = line.split()
        except Exception, e:
            print "Line %d: %s: %s" % (lineno, e, line),

        s = afwDetection.Source()
        sourceSet.append(s)

        s.setId(int(id))
        s.setFlagForDetection(int(flags))
        s.setRa(float(ra))
        s.setXAstrom(float(x))
        s.setYAstrom(float(y))
        s.setDec(float(dec))
        s.setPsfFlux(float(cts))

    return sourceSet
