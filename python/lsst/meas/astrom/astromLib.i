// -*- lsst-c++ -*-
%define astromLib_DOCSTRING
"
Python interface to lsst::meas::astrom
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom", docstring=astromLib_DOCSTRING) astromLib

%{
#include "lsst/afw/image.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/table.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/astrom/matchOptimisticB.h"
#include "lsst/meas/astrom/makeMatchStatistics.h"
#include "lsst/meas/astrom/PolynomialTransform.h"
#include "lsst/meas/astrom/SipTransform.h"
#include "lsst/meas/astrom/ScaledPolynomialTransformFitter.h"
%}

%include "lsst/p_lsstSwig.i"
%initializeNumPy(meas_astrom)
%{
#include "ndarray/converter.h"
#include "ndarray/converter/eigen.h"
%}
%include "ndarray.i"
%include "lsst/pex/config.h"

%shared_ptr(lsst::meas::astrom::MatchOptimisticBControl);

%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/wcs.i"
%import "lsst/afw/math/statistics.i"
%import "lsst/afw/geom/geomLib.i"

%include "lsst/meas/astrom/matchOptimisticB.h"

%include "std_pair.i"

namespace std { // have to do this in namespace block so Swig can find size_t
%template(DoubleSizePair) pair<double,size_t>;
}

%declareNumPyConverters(ndarray::Array<double const,2,2>)

// Avoid dangling references by returning by value in Python.
%returnCopy(lsst::meas::astrom::ScaledPolynomialTransform::getPoly)
%returnCopy(lsst::meas::astrom::ScaledPolynomialTransform::getInputScaling)
%returnCopy(lsst::meas::astrom::ScaledPolynomialTransform::getOutputScalingInverse)
%include "lsst/meas/astrom/PolynomialTransform.h"

%returnCopy(lsst::meas::astrom::SipTransformBase::getPoly)
%returnCopy(lsst::meas::astrom::SipTransformBase::getPixelOrigin)
%returnCopy(lsst::meas::astrom::SipTransformBase::getCDMatrix)
%include "lsst/meas/astrom/SipTransform.h"

%returnCopy(lsst::meas::astrom::ScaledPolynomialTransformFitter::getPoly)
%returnCopy(lsst::meas::astrom::ScaledPolynomialTransformFitter::getInputScaling)
%returnCopy(lsst::meas::astrom::ScaledPolynomialTransformFitter::getOutputScaling)
%returnCopy(lsst::meas::astrom::ScaledPolynomialTransformFitter::getTransform)
%returnCopy(lsst::meas::astrom::ScaledPolynomialTransformFitter::getData)
%include "lsst/meas/astrom/ScaledPolynomialTransformFitter.h"

%include "lsst/meas/astrom/makeMatchStatistics.h"

%define %declareMakeMatchStatisticsFunctions(MATCH)
%template(makeMatchStatistics) lsst::meas::astrom::makeMatchStatistics<MATCH>;
%template(makeMatchStatisticsInPixels) lsst::meas::astrom::makeMatchStatisticsInPixels<MATCH>;
%template(makeMatchStatisticsInRadians) lsst::meas::astrom::makeMatchStatisticsInRadians<MATCH>;
%enddef

%declareMakeMatchStatisticsFunctions(lsst::afw::table::ReferenceMatch);
%declareMakeMatchStatisticsFunctions(lsst::afw::table::SourceMatch);
