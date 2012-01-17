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
#undef FALSE
#undef TRUE
	}
#include <vector>
#include "boost/shared_ptr.hpp"

#include "lsst/base.h"
#include "lsst/pex/logging.h"
#include "lsst/daf/persistence.h"
#include "lsst/daf/base.h"
#include "lsst/daf/base/Persistable.h"
//#include "lsst/afw/detection/Source.h"
//#include "lsst/afw/image.h"
#include "lsst/afw.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/AperturePhotometry.h"


	//namespace afwDet = lsst::afw::detection;

	%}

%include "std_string.i"
%include "std_vector.i"
%include "lsst/p_lsstSwig.i"
%import  "lsst/afw/utils.i" 
%include "lsst/daf/base/persistenceMacros.i"
%include "lsst/base.h"
 //%import "lsst/daf/base/baseLib.i"
%import "lsst/daf/base/Persistable.h"

%lsst_exceptions();

//%shared_ptr(lsst::afw::detection::Source);
//%template(SourceSet) std::vector<boost::shared_ptr<lsst::afw::detection::Source> >;

//%import "lsst/afw/detection/detectionLib.i"
//%import "lsst/afw/image/imageLib.i"

//%import "detlib.i"

//%import "lsst/afw/detection/Source.h"

//%ignore PersistableSourceVector;
%import "lsst/afw/detection/source.i"
 //%import "lsst/afw/detection/detectionLib.i"
 //%import "lsst/afw/detection/Source.h"


%include "solver.h"

%extend solver_t {

	//void setStars(afwDet::SourceSet srcs) {
	void setStars(lsst::afw::detection::SourceSet srcs) {
		// convert to Astrometry.net "starxy_t"
		const size_t N = srcs.size();
		starxy_t *starxy = starxy_new(N, true, false);
		for (size_t i=0; i<N; ++i) {
			double const x    = srcs[i]->getXAstrom();
			double const y    = srcs[i]->getYAstrom();
			double const flux = srcs[i]->getPsfFlux();
			starxy_set(starxy, i, x, y);
			starxy_set_flux(starxy, i, flux);
		}
		// Sort the array
		starxy_sort_by_flux(starxy);

		starxy_free(solver_get_field($self));
		solver_free_field($self);
		solver_set_field($self, starxy);
		solver_reset_field_size($self);
		// Find field boundaries and precompute kdtree
		solver_preprocess_field($self);
	}


 }
