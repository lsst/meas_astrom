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
# \package lsst.meas.photocal
import math, os, sys
import numpy as np

from lsst.meas.photocal.colorterms import Colorterm
import lsst.meas.algorithms.utils as malgUtil
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConf
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
from lsst.meas.base.base import Version0FlagMapper

def checkSourceFlags(source, keys):
    """!Return True if the given source has all good flags set and none of the bad flags set.

    \param[in] source    SourceRecord object to process.
    \param[in] keys      Struct of source catalog keys, as returned by PhotCalTask.getKeys()
    """
    for k in keys.goodFlags:
        if not source.get(k): return False
    if source.getPsfFluxFlag(): return False
    for k in keys.badFlags:
        if source.get(k): return False
    return True

class PhotoCalConfig(pexConf.Config):

    magLimit = pexConf.Field(dtype=float, doc="Don't use objects fainter than this magnitude", default=22.0)
    applyColorTerms = pexConf.Field(
        dtype=bool, default=True,
        doc= "Apply photometric colour terms (if available) to reference stars",
        )
    doWriteOutput = pexConf.Field(
        dtype=bool, default=True,
        doc= "Write a field name astrom_usedByPhotoCal to the schema",
        )
    goodFlags = pexConf.ListField(
        dtype=str, optional=False,
        default=[],
        doc="List of source flag fields that must be set for a source to be used."
        )
    badFlags = pexConf.ListField(
        dtype=str, optional=False,
        default=["base_PixelFlags_flag_edge", "base_PixelFlags_flag_interpolated", "base_PixelFlags_flag_saturated"], 
        doc="List of source flag fields that will cause a source to be rejected when they are set."
        )
    sigmaMax = pexConf.Field(dtype=float, default=0.25, optional=True,
                              doc="maximum sigma to use when clipping")
    nSigma = pexConf.Field(dtype=float, default=3.0, optional=False, doc="clip at nSigma")
    useMedian = pexConf.Field(dtype=bool, default=True,
                              doc="use median instead of mean to compute zeropoint")
    nIter = pexConf.Field(dtype=int, default=20, optional=False, doc="number of iterations")

## \addtogroup LSST_task_documentation
## \{
## \page photoCalTask
## \ref PhotoCalTask_ "PhotoCalTask"
##      Detect positive and negative sources on an exposure and return a new SourceCatalog.
## \}

class PhotoCalTask(pipeBase.Task):
    """!
\anchor PhotoCalTask_

\brief Calculate the zero point of an exposure given a lsst.afw.table.ReferenceMatchVector.

\section meas_photocal_photocal_Contents Contents

 - \ref meas_photocal_photocal_Purpose
 - \ref meas_photocal_photocal_Initialize
 - \ref meas_photocal_photocal_IO
 - \ref meas_photocal_photocal_Config
 - \ref meas_photocal_photocal_Debug
 - \ref meas_photocal_photocal_Example

\section meas_photocal_photocal_Purpose	Description

\copybrief PhotoCalTask

Calculate an Exposure's zero-point given a set of flux measurements of stars matched to an input catalogue.
The type of flux to use is specified by PhotoCalConfig.fluxField.

The algorithm clips outliers iteratively, with parameters set in the configuration.

\note This task can adds fields to the schema, so any code calling this task must ensure that
these columns are indeed present in the input match list; see \ref meas_photocal_photocal_Example

\section meas_photocal_photocal_Initialize	Task initialisation

\copydoc init

\section meas_photocal_photocal_IO		Inputs/Outputs to the run method

\copydoc run

\section meas_photocal_photocal_Config       Configuration parameters

See \ref PhotoCalConfig

\section meas_photocal_photocal_Debug		Debug variables

The \link lsst.pipe.base.cmdLineTask.CmdLineTask command line task\endlink interface supports a
flag \c -d to import \b debug.py from your \c PYTHONPATH; see \ref baseDebug for more about \b debug.py files.

The available variables in PhotoCalTask are:
<DL>
  <DT> \c display
  <DD> If True enable other debug outputs
  <DT> \c displaySources
  <DD> If True, display the exposure on ds9's frame 1 and overlay the source catalogue:
    <DL>
      <DT> red x
      <DD> Bad objects
      <DT> blue +
      <DD> Matched objects deemed unsuitable for photometric calibration.
            Additional information is:
	    - a cyan o for galaxies
	    - a magenta o for variables
      <DT> magenta *
      <DD> Objects that failed the flux cut
      <DT> green o
      <DD> Objects used in the photometric calibration
    </DL>
  <DT> \c scatterPlot
  <DD> Make a scatter plot of flux v. reference magnitude as a function of reference magnitude.
    - good objects in blue
    - rejected objects in red
  (if \c scatterPlot is 2 or more, prompt to continue after each iteration)
</DL>

\section meas_photocal_photocal_Example	A complete example of using PhotoCalTask

This code is in \link photoCalTask.py\endlink in the examples directory, and can be run as \em e.g.
\code
examples/photoCalTask.py
\endcode
\dontinclude photoCalTask.py

Import the tasks (there are some other standard imports; read the file for details)
\skipline from lsst.pipe.tasks.astrometry
\skipline measPhotocal

We need to create both our tasks before processing any data as the task constructors
can add extra columns to the schema which we get from the input catalogue, \c scrCat:
\skipline getSchema

Astrometry first:
\skip AstrometryTask.ConfigClass
\until aTask
(that \c filterMap line is because our test code doesn't use a filter that the reference catalogue recognises,
so we tell it to use the \c r band)

Then photometry:
\skip measPhotocal
\until pTask

If the schema has indeed changed we need to add the new columns to the source table
(yes; this should be easier!)
\skip srcCat
\until srcCat = cat

We're now ready to process the data (we could loop over multiple exposures/catalogues using the same
task objects):
\skip matches
\until result

We can then unpack and use the results:
\skip calib
\until np.log

<HR>
To investigate the \ref meas_photocal_photocal_Debug, put something like
\code{.py}
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
        if name == "lsst.meas.photocal.PhotoCal":
            di.display = 1

        return di

    lsstDebug.Info = DebugInfo
\endcode
into your debug.py file and run photoCalTask.py with the \c --debug flag.
    """
    ConfigClass = PhotoCalConfig
    _DefaultName = "photoCal"

    # Need init as well as __init__ because "\copydoc __init__" fails (doxygen bug 732264)
    def init(self, schema, tableVersion=0, **kwds):
        """!Create the photometric calibration task.  Most arguments are simply passed onto pipe.base.Task.

        \param schema An lsst::afw::table::Schema used to create the output lsst.afw.table.SourceCatalog
        \param tableVersion argument to indicate which afw::table table version 
        \param **kwds keyword arguments to be passed to the lsst.pipe.base.task.Task constructor

        """
        self.__init__(schema, tableVersion, **kwds)

    def __init__(self, schema, tableVersion=0, **kwds):
        """!Create the photometric calibration task.  See PhotoCalTask.init for documentation
        """
        pipeBase.Task.__init__(self, **kwds)
        if self.config.doWriteOutput:
            if tableVersion == 0:
                self.outputField = schema.addField("classification.photometric", type="Flag",
                                          doc="set if source was used in photometric calibration")
            else:
                self.outputField = schema.addField("photocal_photometricStandard", type="Flag",
                                          doc="set if source was used in photometric calibration")
        else:
            self.outputField = None

    def getKeys(self, schema, tableVersion=0):
        """!Return a struct containing the source catalog keys for fields used by PhotoCalTask."""

        if tableVersion == 0:
            goodFlags = [schema.find(name).key for name in Version0FlagMapper(self.config.goodFlags)]
            badFlags = [schema.find(name).key for name in Version0FlagMapper(self.config.badFlags)]
        else:
            goodFlags = [schema.find(name).key for name in self.config.goodFlags]
            badFlags = [schema.find(name).key for name in self.config.badFlags]
        return pipeBase.Struct(goodFlags=goodFlags, badFlags=badFlags)

    @pipeBase.timeMethod
    def selectMatches(self, matches, keys, frame=None):
        """!Select reference/source matches according the criteria specified in the config.

        \param[in] matches ReferenceMatchVector (not modified)
        \param[in] keys    Struct of source catalog keys, as returned by getKeys()
        \param[in] frame   ds9 frame number to use for debugging display
        if frame is non-None, display information about trimmed objects on that ds9 frame:
         - Bad:               red x
         - Unsuitable objects: blue +  (and a cyan o if a galaxy)
         - Failed flux cut:   magenta *

        \return a \link lsst.afw.table.ReferenceMatchVector\endlink that contains only the selected matches.
        If a schema was passed during task construction, a flag field will be set on sources
        in the selected matches.

        \throws ValueError There are no valid matches.
        """

        self.log.logdebug("Number of input matches: %d" % (len(matches)))
        if len(matches) == 0:
            raise ValueError("No input matches")

        # Only use stars for which the flags indicate the photometry is good.
        afterFlagCutInd = [i for i, m in enumerate(matches) if checkSourceFlags(m.second, keys)]
        afterFlagCut = [matches[i] for i in afterFlagCutInd]
        self.log.logdebug("Number of matches after source flag cuts: %d" % (len(afterFlagCut)))

        if len(afterFlagCut) != len(matches):
            if frame is not None:
                with ds9.Buffering():
                    for i, m in enumerate(matches):
                        if i not in afterFlagCutInd:
                            x, y = m.second.getCentroid()
                            ds9.dot("x", x,  y, size=4, frame=frame, ctype=ds9.RED)

            matches = afterFlagCut

        if len(matches) == 0:
            raise ValueError("All matches eliminated by source flags")

        refSchema = matches[0].first.schema
        try:
            refKey = refSchema.find("photometric").key
            try:
                stargalKey = refSchema.find("stargal").key
            except:
                stargalKey = None

            try:
                varKey = refSchema.find("var").key
            except:
                varKey = None
        except:
            self.log.warn("No 'photometric' flag key found in reference schema.")
            refKey = None

        if refKey is not None:
            afterRefCutInd = [i for i, m in enumerate(matches) if m.first.get(refKey)]
            afterRefCut = [matches[i] for i in afterRefCutInd]

            if len(afterRefCut) != len(matches):
                if frame is not None:
                    with ds9.Buffering():
                        for i, m in enumerate(matches):
                            if i not in afterRefCutInd:
                                x, y = m.second.getCentroid()
                                ds9.dot("+", x,  y, size=4, frame=frame, ctype=ds9.BLUE)

                                if stargalKey and not m.first.get(stargalKey):
                                    ds9.dot("o", x,  y, size=6, frame=frame, ctype=ds9.CYAN)
                                if varKey and m.first.get(varKey):
                                    ds9.dot("o", x,  y, size=6, frame=frame, ctype=ds9.MAGENTA)

                matches = afterRefCut

        self.log.logdebug("Number of matches after reference catalog cuts: %d" % (len(matches)))
        if len(matches) == 0:
            raise RuntimeError("No sources remain in match list after reference catalog cuts.")
        fluxKey = refSchema.find("flux").key
        if self.config.magLimit is not None:
            fluxLimit = 10.0**(-self.config.magLimit/2.5)

            afterMagCutInd = [i for i, m in enumerate(matches) if (m.first.get(fluxKey) > fluxLimit
                                                                   and m.second.getPsfFlux() > 0.0)]
        else:
            afterMagCutInd = [i for i, m in enumerate(matches) if m.second.getPsfFlux() > 0.0]

        afterMagCut = [matches[i] for i in afterMagCutInd]

        if len(afterMagCut) != len(matches):
            if frame is not None:
                with ds9.Buffering():
                    for i, m in enumerate(matches):
                        if i not in afterMagCutInd:
                            x, y = m.second.getCentroid()
                            ds9.dot("*", x,  y, size=4, frame=frame, ctype=ds9.MAGENTA)

            matches = afterMagCut

        self.log.logdebug("Number of matches after magnitude limit cuts: %d" % (len(matches)))

        if len(matches) == 0:
            raise RuntimeError("No sources remaining in match list after magnitude limit cuts.")

        if frame is not None:
            with ds9.Buffering():
                for m in matches:
                    x, y = m.second.getCentroid()
                    ds9.dot("o", x,  y, size=4, frame=frame, ctype=ds9.GREEN)

        result = afwTable.ReferenceMatchVector()
        for m in matches:
            if self.outputField is not None:
                m.second.set(self.outputField, True)
            result.append(m)
        return result

    @pipeBase.timeMethod
    def extractMagArrays(self, matches, filterName, keys):
        """!Extract magnitude and magnitude error arrays from the given matches.

        \param[in] matches     \link lsst::afw::table::ReferenceMatchVector\endlink object containing reference/source matches
        \param[in] filterName Name of filter being calibrated
        \param[in] keys       Struct of source catalog keys, as returned by getKeys()

        \return Struct containing srcMag, refMag, srcMagErr, refMagErr, and magErr numpy arrays
        where magErr is an error in the magnitude; the error in srcMag - refMag.
        \note These are the \em inputs to the photometric calibration, some may have been
        discarded by clipping while estimating the calibration (https://jira.lsstcorp.org/browse/DM-813)
        """
        srcFlux = np.array([m.second.getPsfFlux() for m in matches])
        srcFluxErr = np.array([m.second.getPsfFluxErr() for m in matches])
        if not np.all(np.isfinite(srcFluxErr)):
            self.log.warn("Source catalog does not have flux uncertainties; using sqrt(flux).")
            srcFluxErr = np.sqrt(srcFlux)

        if not matches:
            raise RuntimeError("No reference stars are available")

        refSchema = matches[0].first.schema
        if self.config.applyColorTerms:
            ct = Colorterm.getColorterm(filterName)
        else:
            ct = None

        if ct:                          # we have a colour term to worry about
            fluxNames = [ct.primary, ct.secondary]
            missingFluxes = []
            for flux in fluxNames:
                try:
                    refSchema.find(flux).key
                except KeyError:
                    missingFluxes.append(flux)

            if missingFluxes:
                self.log.warn("Source catalog does not have fluxes for %s; ignoring color terms" %
                              " ".join(missingFluxes))
                ct = None

        if not ct:
            fluxNames = ["flux"]

        refFluxes = []
        refFluxErrors = []
        for flux in fluxNames:
            fluxKey = refSchema.find(flux).key
            refFlux = np.array([m.first.get(fluxKey) for m in matches])
            try:
                fluxErrKey = refSchema.find(flux + ".err").key
                refFluxErr = np.array([m.first.get(fluxErrKey) for m in matches])
            except KeyError:
                # Catalogue may not have flux uncertainties; HACK
                self.log.warn("Reference catalog does not have flux uncertainties for %s; using sqrt(flux)."
                              % flux)
                refFluxErr = np.sqrt(refFlux)

            refFluxes.append(refFlux)
            refFluxErrors.append(refFluxErr)

        srcMag = -2.5*np.log10(srcFlux)
        if ct:                          # we have a colour term to worry about
            refMag =  -2.5*np.log10(refFluxes[0]) # primary
            refMag2 = -2.5*np.log10(refFluxes[1]) # secondary

            refMag = ct.transformMags(filterName, refMag, refMag2)
            refFluxErr = ct.propagateFluxErrors(filterName, refFluxErrors[0], refFluxErrors[1])
        else:
            refMag = -2.5*np.log10(refFluxes[0])

        # Fitting with error bars in both axes is hard, so transfer all
        # the error to src, then convert to magnitude
        fluxErr = np.hypot(srcFluxErr, refFluxErr)
        magErr = fluxErr / (srcFlux * np.log(10))

        srcMagErr = srcFluxErr / (srcFlux * np.log(10))
        refMagErr = refFluxErr / (refFlux * np.log(10))

        return pipeBase.Struct(
            srcMag = srcMag,
            refMag = refMag,
            magErr = magErr,
            srcMagErr = srcMagErr,
            refMagErr = refMagErr
            )

    @pipeBase.timeMethod
    def run(self, exposure, matches):
        """!Do photometric calibration - select matches to use and (possibly iteratively) compute
        the zero point.

        \param[in]  exposure   Exposure upon which the sources in the matches were detected.
        \param[in]  matches    Input lsst.afw.table.ReferenceMatchVector
        (\em i.e. a list of lsst.afw.table.Match with
        \c first being of type lsst.afw.table.SimpleRecord and \c second type lsst.afw.table.SourceRecord ---
        the reference object and matched object respectively).
        (will not be modified  except to set the outputField if requested.).

        \return Struct of:
         - calib -------  \link lsst::afw::image::Calib\endlink object containing the zero point
         - arrays ------ Magnitude arrays returned be PhotoCalTask.extractMagArrays
         - matches ----- Final ReferenceMatchVector, as returned by PhotoCalTask.selectMatches.

The exposure is only used to provide the name of the filter being calibrated (it may also be
used to generate debugging plots).

The reference objects:
 - Must include a field \c photometric; True for objects which should be considered as photometric standards
 - Must include a field \c flux; the flux used to impose a magnitude limit and also to calibrate the data to (unless a colour term is specified, in which case ColorTerm.primary is used;  See https://jira.lsstcorp.org/browse/DM-933)
 - May include a field \c stargal; if present, True means that the object is a star
 - May include a field \c var; if present, True means that the object is variable

The measured sources:
- Must include PhotoCalConfig.fluxField; the flux measurement to be used for calibration

\throws RuntimeError with the following strings:

<DL>
<DT> `sources' schema does not contain the calibration object flag "XXX"`
<DD> The constructor added fields to the schema that aren't in the Sources
<DT> No input matches
<DD> The input match vector is empty
<DT> All matches eliminated by source flags
<DD> The flags specified by \c badFlags in the config eliminated all candidate objects
<DT> No sources remain in match list after reference catalog cuts
<DD> The reference catalogue has a column "photometric", but no matched objects have it set
<DT> No sources remaining in match list after magnitude limit cuts
<DD> All surviving matches are either too faint in the catalogue or have negative or \c NaN flux
<DT> No reference stars are available
<DD> No matches survive all the checks
</DL>

        """
        global scatterPlot, fig
        import lsstDebug

        display = lsstDebug.Info(__name__).display
        displaySources = display and lsstDebug.Info(__name__).displaySources
        scatterPlot = display and lsstDebug.Info(__name__).scatterPlot

        if scatterPlot:
            from matplotlib import pyplot
            try:
                fig.clf()
            except:
                fig = pyplot.figure()

        if displaySources:
            frame = 1
            ds9.mtv(exposure, frame=frame, title="photocal")
        else:
            frame = None

        keys = self.getKeys(matches[0].second.schema, matches[0].second.getTable().getVersion())
        matches = self.selectMatches(matches, keys, frame=frame)
        arrays = self.extractMagArrays(matches, exposure.getFilter().getName(), keys)

        if matches and self.outputField:
            try:
                # matches[].second is a measured source, wherein we wish to set outputField.
                # Check that the field is present in the Sources schema.
                matches[0].second.getSchema().find(self.outputField)
            except:
                raise RuntimeError("sources' schema does not contain the used-in-calibration flag \"%s\"" %
                                   self.config.outputField)

        # Fit for zeropoint.  We can run the code more than once, so as to
        # give good stars that got clipped by a bad first guess a second
        # chance.
        # FIXME: these should be config values

        calib = afwImage.Calib()
        zp = None                           # initial guess
        r = self.getZeroPoint(arrays.srcMag, arrays.refMag, arrays.magErr, zp0=zp)
        zp = r.zp
        self.log.info("Magnitude zero point: %f +/- %f from %d stars" % (r.zp, r.sigma, r.ngood))

        flux0 = 10**(0.4*r.zp) # Flux of mag=0 star
        flux0err = 0.4*math.log(10)*flux0*r.sigma # Error in flux0

        calib.setFluxMag0(flux0, flux0err)

        return pipeBase.Struct(
            calib = calib,
            arrays = arrays,
            matches = matches,
            zp = r.zp,
            sigma = r.sigma,
            ngood = r.ngood,
            )

    def getZeroPoint(self, src, ref, srcErr=None, zp0=None):
        """!Flux calibration code, returning (ZeroPoint, Distribution Width, Number of stars)

        We perform nIter iterations of a simple sigma-clipping algorithm with a a couple of twists:
        1.  We use the median/interquartile range to estimate the position to clip around, and the
        "sigma" to use.
        2.  We never allow sigma to go _above_ a critical value sigmaMax --- if we do, a sufficiently
        large estimate will prevent the clipping from ever taking effect.
        3.  Rather than start with the median we start with a crude mode.  This means that a set of magnitude
        residuals with a tight core and asymmetrical outliers will start in the core.  We use the width of
        this core to set our maximum sigma (see 2.)
        """

        sigmaMax = self.config.sigmaMax

        dmag = ref - src

        i = np.argsort(dmag)
        dmag = dmag[i]

        if srcErr is not None:
            dmagErr = srcErr[i]
        else:
            dmagErr = np.ones(len(dmag))

        # need to remove nan elements to avoid errors in stats calculation with numpy
        ind_noNan = np.array([ i for i in range(len(dmag)) if (not np.isnan(dmag[i]) and not np.isnan(dmagErr[i])) ])
        dmag = dmag[ind_noNan]
        dmagErr = dmagErr[ind_noNan]

        IQ_TO_STDEV = 0.741301109252802;    # 1 sigma in units of interquartile (assume Gaussian)

        npt = len(dmag)
        ngood = npt
        for i in range(self.config.nIter):
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

                    if sigmaMax is None:
                        sigmaMax = 2*sig   # upper bound on st. dev. for clipping. multiplier is a heuristic

                    self.log.logdebug("Photo calibration histogram: center = %.2f, sig = %.2f"
                                      % (center, sig))

                else:
                    if sigmaMax is None:
                        sigmaMax = dmag[-1] - dmag[0]

                    center = np.median(dmag)
                    q1 = dmag[int(0.25*npt)]
                    q3 = dmag[int(0.75*npt)]
                    sig = (q3 - q1)/2.3 # estimate of standard deviation (based on FWHM; 2.358 for Gaussian)

            if center is None:              # usually equivalent to (i > 0)
                gdmag = dmag[good]
                if self.config.useMedian:
                    center = np.median(gdmag)
                else:
                    gdmagErr = dmagErr[good]
                    center = np.average(gdmag, weights=gdmagErr)

                q3 = gdmag[min(int(0.75*npt + 0.5), npt - 1)]
                q1 = gdmag[min(int(0.25*npt + 0.5), npt - 1)]

                sig = IQ_TO_STDEV*(q3 - q1)     # estimate of standard deviation

            good = abs(dmag - center) < self.config.nSigma*min(sig, sigmaMax) # don't clip too softly

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            if scatterPlot:
                from matplotlib import pyplot
                try:
                    fig.clf()

                    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

                    axes.plot(ref[good], dmag[good] - center, "b+")
                    axes.errorbar(ref[good], dmag[good] - center, yerr=dmagErr[good],
                                  linestyle='', color='b')

                    bad = np.logical_not(good)
                    if len(ref[bad]) > 0:
                        axes.plot(ref[bad], dmag[bad] - center, "r+")
                        axes.errorbar(ref[bad], dmag[bad] - center, yerr=dmagErr[bad],
                                      linestyle='', color='r')

                    axes.plot((-100, 100), (0, 0), "g-")
                    for x in (-1, 1):
                        axes.plot((-100, 100), x*0.05*np.ones(2), "g--")

                    axes.set_ylim(-1.1, 1.1)
                    axes.set_xlim(24, 13)
                    axes.set_xlabel("Reference")
                    axes.set_ylabel("Reference - Instrumental")

                    fig.show()

                    if scatterPlot > 1:
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

                self.log.log(self.log.WARN, msg)

                return pipeBase.Struct(
                    zp = center,
                    sigma = sig,
                    ngood = len(dmag)
                    )
            elif ngood == old_ngood:
                break

            if False:
                ref = ref[good]
                dmag = dmag[good]
                dmagErr = dmagErr[good]

        dmag = dmag[good]
        dmagErr = dmagErr[good]
        return pipeBase.Struct(
            zp = np.average(dmag, weights=dmagErr),
            sigma = np.std(dmag, ddof=1),
            ngood = len(dmag)
            )
