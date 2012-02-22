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
import numpy as np

import lsst.meas.astrom as measAst
import lsst.meas.astrom.sip as sip
import lsst.meas.algorithms.utils as malgUtil
import lsst.pex.logging as pexLog

from lsst.meas.photocal.PhotometricMagnitude import PhotometricMagnitude

try:
    import matplotlib.pyplot as pyplot
except ImportError:
    pyplot = None

def calcPhotoCal(sourceMatch, log=None, magLimit=22, useCatalogClassification=True,
                 goodFlagValue=malgUtil.getDetectionFlags()['BINNED1'],
                 badFlagValue=malgUtil.getDetectionFlags()['BAD'],
                 usePsfFlux=True
                 ):
    """Calculate photometric calibration, i.e the zero point magnitude

    If useCatalogClassification is true, use the star/galaxy classification from the reference catalogue, otherwise
    use the value from the measured sources (specifically, the STAR bit in the detection flags)
    """

    if log is None:
        log = pexLog.Log.getDefaultLog()

    global display, fig
    import lsstDebug
    display = lsstDebug.Info(__name__).display
    if display and pyplot:
        try:
            fig.clf()
        except:
            fig = pyplot.figure()
            
    if len(sourceMatch) == 0:
        raise ValueError("sourceMatch contains no elements")

    origSourceMatch = sourceMatch

    # Only use stars for which the flags indicate the photometry is good.
    log.logdebug("Number of sources: %d" % (len(sourceMatch)))

    Iflags = [i for i,m in enumerate(sourceMatch) if
              (m.second.getFlagForDetection() & goodFlagValue) == goodFlagValue and
              (m.second.getFlagForDetection() & badFlagValue ) == 0]
    I = Iflags
    sourceMatch = [origSourceMatch[i] for i in I]
    log.logdebug("Number of sources with good flag settings: %d" % (len(sourceMatch)))
    if len(sourceMatch) == 0:
        raise ValueError("flags indicate all elements of sourceMatch have bad photometry")

    STAR = malgUtil.getDetectionFlags()['STAR']
    #
    # See if any catalogue objects are labelled as stars; if not use the measured object's classifier
    #
    if useCatalogClassification:
        Istar = [i for i,m in zip(I,sourceMatch)
                 if (m.first.getFlagForDetection() & STAR)]
        if len(Istar) == 0:
            log.warn('No catalog sources are classified as stars; using measured source STAR flag')
            useCatalogClassification = False
        else:
            I = Istar
            sourceMatch = [origSourceMatch[i] for i in I]

    if not useCatalogClassification:
        Istar = [i for i,m in zip(I,sourceMatch)
                 if (m.second.getFlagForDetection() & STAR)]
        if len(Istar) == 0:
            log.warn("No image objects are classified as stars; using all sources")
        else:
            I = Istar
            sourceMatch = [origSourceMatch[i] for i in I]

    log.logdebug("Number of sources after stellar cuts: %d" % (len(sourceMatch)))
    if len(sourceMatch) == 0:
        raise RuntimeError("No sources remaining in match list after cuts")
 
    # Convert fluxes to magnitudes
    fluxCat    = np.array([m.first.getPsfFlux()    for m in sourceMatch])
    fluxCatErr = np.array([m.first.getPsfFluxErr() for m in sourceMatch])
    if usePsfFlux:
        fluxSrc =    np.array([m.second.getPsfFlux()    for m in sourceMatch])
        fluxSrcErr = np.array([m.second.getPsfFluxErr() for m in sourceMatch])
    else:
        fluxSrc =    np.array([m.second.getApFlux()     for m in sourceMatch])
        fluxSrcErr = np.array([m.second.getApFluxErr()  for m in sourceMatch])

    keep = np.logical_and(fluxCat > 0, fluxSrc > 0)
    Iflux = [I[k] for k in np.flatnonzero(keep)]
    I = Iflux
    sourceMatch = [origSourceMatch[i] for i in I]

    # Catalogue may not have flux uncertainties; HACK
    if sum(fluxCatErr) == 0.0:
        fluxCatErr = np.sqrt(fluxCat)
    magSrc = -2.5 * np.log10(fluxSrc[keep])
    magCat = -2.5 * np.log10(fluxCat[keep])

    # Fitting with error bars in both axes is hard, so transfer all
    # the error to src, then convert to magnitude
    fluxErr = np.hypot(fluxSrcErr[keep], fluxCatErr[keep])
    magErr = fluxErr/fluxSrc[keep]/np.log(10)

    magSrcErr = fluxSrcErr[keep]/fluxSrc[keep]/np.log(10)
    magCatErr = fluxCatErr[keep]/fluxCat[keep]/np.log(10)

    Ibright = None
    if magLimit is not None:
        keep = (magCat < magLimit)
        Ibright = [I[k] for k in np.flatnonzero(keep)]
        I = Ibright
        magSrc = magSrc[keep]
        magCat = magCat[keep]
        magErr = magErr[keep]
        magSrcErr = magSrcErr[keep]
        magCatErr = magCatErr[keep]
        
    # Fit for zeropoint.  We can run the code more than once, so as to
    # give good stars that got clipped by a bad first guess a second
    # chance.
    sigma_max = [0.25]                  # maximum sigma to use when clipping
    nsigma = [3.0]                      # clip at nsigma
    useMedian = [True]
    niter = [20]                        # number of iterations

    zp = None                           # initial guess
    for i in range(len(nsigma)):
        zp, sigma, ngood = getZeroPoint(magSrc, magCat, magErr, zp0=zp,
                                        useMedian=useMedian[i], sigma_max=sigma_max[i],
                                        nsigma=nsigma[i], niter=niter[i], log=log)
        log.info("Magnitude zero point: %f +/- %f from %d stars" % (zp, sigma, ngood))

    photocal = PhotometricMagnitude(zeroFlux=1.0, zeroMag=zp)
    # add debugging data
    photocal.srcMag = magSrc
    photocal.srcMagErr = magSrcErr
    photocal.refMag = magCat
    photocal.refMagErr = magCatErr
    photocal.Iflags = Iflags
    photocal.Istar = Istar
    photocal.Iflux = Iflux
    photocal.Ibright = Ibright

    return photocal

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def getZeroPoint(src, ref, srcErr=None, zp0=None, useMedian=True, sigma_max=None, nsigma=2, niter=3, log=None):
    """Flux calibration code, returning (ZeroPoint, Distribution Width, Number of stars)

We perform niter iterations of a simple sigma-clipping algorithm with a a couple of twists:
  1.  We use the median/interquartile range to estimate the position to clip around, and the
  "sigma" to use
  2.  We never allow sigma to go _above_ a critical value sigma_max --- if we do, a sufficiently large estimate
  will prevent the clipping from ever taking effect
  3.  Rather than start with the median we start with a crude mode.  This means that a set of magnitude
  residuals with a tight core and asymmetrical outliers will start in the core.  We use the width of
  this core to set our maximum sigma (see 2.)  
    """

    dmag = ref - src
    i = np.argsort(dmag)
    dmag = dmag[i]
    if srcErr is not None:
        dmagErr = srcErr[i]
    else:
        dmagErr = np.ones(len(dmag))

    IQ_TO_STDEV = 0.741301109252802;    # 1 sigma in units of interquartile (assume Gaussian)

    npt = len(dmag)
    ngood = npt
    for i in range(niter):
        if i > 0:
            npt = sum(good)

        center = None
        if i == 0:
            #
            # Start by finding the mode
            #
            nhist = 20
            try:
                hist, edges = np.histogram(dmag, nhist, new=True)
            except TypeError:
                hist, edges = np.histogram(dmag, nhist) # they removed new=True around numpy 1.5
            imode = np.arange(nhist)[np.where(hist == hist.max())]

            if imode[-1] - imode[0] + 1 == len(imode): # Multiple modes, but all contiguous
                if zp0:
                    center = zp0
                else:
                    center = 0.5*(edges[imode[0]] + edges[imode[-1] + 1])

                peak = sum(hist[imode])/len(imode) # peak height

                # Estimate FWHM of mode
                j = imode[0]
                while j >= 0 and hist[j] > 0.5*peak:
                    j -= 1
                j = max(j, 0)
                q1 = dmag[sum(hist[range(j)])]

                j = imode[-1]
                while j < nhist and hist[j] > 0.5*peak:
                    j += 1
                j = min(j, nhist - 1)
                j = min(sum(hist[range(j)]), npt - 1)
                q3 = dmag[j]

                if q1 == q3:
                    q1 = dmag[int(0.25*npt)]
                    q3 = dmag[int(0.75*npt)]

                sig = (q3 - q1)/2.3 # estimate of standard deviation (based on FWHM; 2.358 for Gaussian)

                if sigma_max is None:
                    sigma_max = 2*sig   # upper bound on st. dev. for clipping. multiplier is a heuristic

                if log:
                    log.logdebug("Photo calibration histogram: center = %.2f, sig = %.2f" % (center, sig))

            else:
                if sigma_max is None:
                    sigma_max = dmag[-1] - dmag[0]

                center = np.median(dmag)
                q1 = dmag[int(0.25*npt)]
                q3 = dmag[int(0.75*npt)]
                sig = (q3 - q1)/2.3 # estimate of standard deviation (based on FWHM; 2.358 for Gaussian)

        if center is None:              # usually equivalent to (i > 0)
            gdmag = dmag[good]
            if useMedian:
                center = np.median(gdmag)
            else:
                gdmagErr = dmagErr[good]
                center = np.average(gdmag, weights=gdmagErr)
            
            q3 = gdmag[min(int(0.75*npt + 0.5), npt - 1)]
            q1 = gdmag[min(int(0.25*npt + 0.5), npt - 1)]
        
            sig = IQ_TO_STDEV*(q3 - q1)     # estimate of standard deviation

        good = abs(dmag - center) < nsigma*min(sig, sigma_max) # don't clip too softly

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        if display:
            if not pyplot:
                print >> sys.stderr, "I am unable to plot as I failed to import matplotlib"
            else:
                try:
                    fig.clf()

                    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

                    axes.plot(ref[good], dmag[good] - center, "b+")
                    axes.errorbar(ref[good], dmag[good] - center, yerr=dmagErr[good], linestyle='', color='b')

                    bad = np.logical_not(good)
                    if len(ref[bad]) > 0:
                        axes.plot(ref[bad], dmag[bad] - center, "r+")
                        axes.errorbar(ref[bad], dmag[bad] - center, yerr=dmagErr[bad], linestyle='', color='r')

                    axes.plot((-100, 100), (0, 0), "g-")
                    for x in (-1, 1):
                        axes.plot((-100, 100), x*0.05*np.ones(2), "g--")

                    axes.set_ylim(-1.1, 1.1)
                    axes.set_xlim(24, 13)
                    axes.set_xlabel("Reference")
                    axes.set_ylabel("Reference - Instrumental")

                    fig.show()
                    
                    if display > 1:
                        while i == 0 or reply != "c":
                            try:
                                reply = raw_input("Next iteration? [ynhpc] ")
                            except EOFError:
                                reply = "n"

                            if reply == "h":
                                print >> sys.stderr, "Options: c[ontinue] h[elp] n[o] p[db] y[es]"
                                continue

                            if reply in ("", "c", "n", "p", "y"):
                                break
                            else:
                                print >> sys.stderr, "Unrecognised response: %s" % reply

                        if reply == "n":
                            break
                        elif reply == "p":
                            import pdb; pdb.set_trace()
                except Exception, e:
                    print >> sys.stderr, "Error plotting in PhotoCal.getZeroPoint: %s" % e

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        old_ngood = ngood
        ngood = sum(good)
        if ngood == 0:
            msg = "PhotoCal.getZeroPoint: no good stars remain"

            if i == 0:                  # failed the first time round -- probably all fell in one bin
                center = np.average(dmag, weights=dmagErr)
                msg += " on first iteration; using average of all calibration stars"

            if log:
                log.log(log.WARN, msg)
            else:
                print >> sys.stderr, msg

            return center, sig, len(dmag)
        elif ngood == old_ngood:
            break

        if False:
            ref = ref[good]
            dmag = dmag[good]
            dmagErr = dmagErr[good]

    dmag = dmag[good]
    dmagErr = dmagErr[good]
    return np.average(dmag, weights=dmagErr), np.std(dmag, ddof=1), len(dmag)
