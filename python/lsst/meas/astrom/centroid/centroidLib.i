// -*- lsst-c++ -*-
%define centroidLib_DOCSTRING
"
Python interface to lsst::afw::meas::astrom::centroid classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.centroid",docstring=centroidLib_DOCSTRING) centroidLib

%{
#   include "lsst/meas/astrom/centroid/Centroid.h"
#   include "lsst/meas/astrom/centroid/NaiveCentroid.h"
%}

%include "lsst/p_lsstSwig.i"
%include "std_vector.i"

%ignore boost::noncopyable;
namespace boost {
    class noncopyable {};
}

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/afw/trunk/python/lsst/afw/centroid/centroidLib.i $"):
    """Return a version given a HeadURL string. If a different version is setup, return that too"""

    return HeadURL
%}

%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"

%lsst_exceptions();

%template(xyAndError) std::pair<double, double>;

%include "lsst/meas/astrom/centroid/Centroid.h"
%include "lsst/meas/astrom/centroid/NaiveCentroid.h"

%template(CentroiderF) lsst::meas::astrom::centroid::Centroider<lsst::afw::image::Image<float> >;
%template(NaiveCentroiderF) lsst::meas::astrom::centroid::NaiveCentroider<lsst::afw::image::Image<float> >;
%template(make_Centroider) lsst::meas::astrom::centroid::make_Centroider<lsst::afw::image::Image<float> >;
