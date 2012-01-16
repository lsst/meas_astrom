%define astrometry_net_DOCSTRING
"
Python interface to Astrometry.net
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.astrometry_net",
		docstring=astrometry_net_DOCSTRING) astrometry_net


%{
	extern "C" {
#include "solver.h"
	}
	%}

%include "solver.h"

