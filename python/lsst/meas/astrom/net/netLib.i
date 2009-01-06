// -*- lsst-c++ -*-
%define netLib_DOCSTRING
"
Python interface to lsst::afw::meas::astrom::net classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.net",docstring=netLib_DOCSTRING) netLib

%{
#   include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
%}

%include "lsst/p_lsstSwig.i"
%include "std_string.i"
%include "std_vector.i"
//Did you declare the std::vector<std::string> using %template?


%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/afw/trunk/python/lsst/afw/net/netLib.i $"):
    """Return a version given a HeadURL string. If a different version is setup, return that too"""

    return HeadURL
%}

%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"

%lsst_exceptions();

%include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
