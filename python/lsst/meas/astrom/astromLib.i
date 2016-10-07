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

%include "lsst/meas/astrom/transforms.i"

%include "lsst/meas/astrom/matchOptimisticB.h"

%include "lsst/meas/astrom/makeMatchStatistics.h"

%define %declareMakeMatchStatisticsFunctions(MATCH)
%template(makeMatchStatistics) lsst::meas::astrom::makeMatchStatistics<MATCH>;
%template(makeMatchStatisticsInPixels) lsst::meas::astrom::makeMatchStatisticsInPixels<MATCH>;
%template(makeMatchStatisticsInRadians) lsst::meas::astrom::makeMatchStatisticsInRadians<MATCH>;
%enddef

%declareMakeMatchStatisticsFunctions(lsst::afw::table::ReferenceMatch);
%declareMakeMatchStatisticsFunctions(lsst::afw::table::SourceMatch);
