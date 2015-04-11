// -*- lsst-c++ -*-
%define astromLib_DOCSTRING
"
Python interface to lsst::meas::astrom
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom", docstring=astromLib_DOCSTRING) astromLib

%{
#include "lsst/pex/logging.h"
#include "lsst/afw/image.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/table.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/astrom/matchOptimisticB.h"
%}

%include "lsst/p_lsstSwig.i"
%initializeNumPy(meas_astrom)

%include "lsst/pex/config.h"

%shared_ptr(lsst::meas::astrom::MatchOptimisticBControl);

%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/wcs.i"

%include "lsst/meas/astrom/matchOptimisticB.h"
