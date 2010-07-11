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

import math, os, sys
import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage

"""This is not a test file, but a convenient way of  reading/writing source sets to and from an ascii file"""


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
