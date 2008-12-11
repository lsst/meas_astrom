// -*- lsst-c++ -*-
%define netLib_DOCSTRING
"
Python interface to lsst::afw::meas::astrom::net classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.net",docstring=netLib_DOCSTRING) netLib

%{
#   include "lsst/daf/base.h"
#   include "lsst/pex/policy/Policy.h"
#   include "lsst/pex/policy/PolicyFile.h"
#   include "lsst/afw/image.h"
#   include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
%}

%include "lsst/p_lsstSwig.i"

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/afw/trunk/python/lsst/afw/net/netLib.i $"):
    """Return a version given a HeadURL string. If a different version is setup, return that too"""

    return HeadURL
%}

//%template(vectorF) std::vector<float>;

%import "lsst/afw/image/imageLib.i"

%lsst_exceptions();

%include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"

