// -*- lsst-c++ -*-
%define netLib_DOCSTRING
"
Python interface to lsst::afw::meas::astrom::net classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.net",docstring=netLib_DOCSTRING) netLib

%{
#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
#include "lsst/afw/geom.h" // should not be needed; ticket #1121
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
//template(PointD) lsst::afw::image::PointD<double,double>;
%include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"

%template(vectorSourceMatch) std::vector<boost::shared_ptr<lsst::afw::detection::SourceMatch> >;
