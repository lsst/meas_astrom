// -*- lsst-c++ -*-

%{
#include "lsst/afw/table.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/astrom/PolynomialTransform.h"
#include "lsst/meas/astrom/SipTransform.h"
#include "lsst/meas/astrom/ScaledPolynomialTransformFitter.h"
%}

%include "std_pair.i"

%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/wcs.i"
%import "lsst/afw/geom/geomLib.i"


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
