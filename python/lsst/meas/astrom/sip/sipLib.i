// -*- lsst-c++ -*-
%define sipLib_DOCSTRING
"
Python interface to lsst::afw::meas::astrom::sip classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.sip",docstring=sipLib_DOCSTRING) sipLib

%{
#   include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"
#   include "lsst/afw/math/FunctionLibrary.h"
#   include "lsst/meas/astrom/sip/LeastSqFitter1d.h"
#   include "lsst/meas/astrom/sip/LeastSqFitter2d.h"
#   include "lsst/meas/astrom/sip/createWcsSip.h"
%}

%include "lsst/p_lsstSwig.i"
%include "std_string.i"
%include "std_vector.i"


%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/afw/trunk/python/lsst/afw/sip/sipLib.i $"):
    """Return a version given a HeadURL string. If a different version is setup, return that too"""

    return HeadURL
%}

%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"

%lsst_exceptions();

%include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"
%include "lsst/afw/math/FunctionLibrary.h"
%include "lsst/meas/astrom/sip/LeastSqFitter1d.h"
%include "lsst/meas/astrom/sip/LeastSqFitter2d.h"
%include "lsst/meas/astrom/sip/createWcsSip.h"


%template(LeastSqFitter1dPoly) lsst::meas::astrom::sip::LeastSqFitter1d< lsst::afw::math::PolynomialFunction1<double> >;

%template(LeastSqFitter2dPoly) lsst::meas::astrom::sip::LeastSqFitter2d< lsst::afw::math::PolynomialFunction1<double> >;

