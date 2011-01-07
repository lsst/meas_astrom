// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
%define netLib_DOCSTRING
"
Python interface to lsst::afw::meas::astrom::net classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.net",docstring=netLib_DOCSTRING) netLib


%{
#include "lsst/pex/logging/BlockTimingLog.h"    
#include "lsst/pex/logging/ScreenLog.h"    
#include "lsst/pex/logging/DualLog.h"    
#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
#include "lsst/afw/geom.h" // should not be needed; ticket #1121

#include "starkd_search_stars.c"

%}

// As in ap/.../apLib.i -- handle swig problems with boost::int64_t
namespace boost {
#if defined(SWIGWORDSIZE64)
    typedef long int64_t;
#else
    typedef long long int64_t;
#endif
}

%include "lsst/p_lsstSwig.i"
%include "std_string.i"
%include "std_vector.i"
 //%include "std_pair.i"

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
%include "qfits_table.h"

 // Ignore astrometry.net declarations that aren't implemented in the library we use
 // (which is a subset of the complete codebase)
%ignore matchobj_compute_overlap;
%ignore matchobj_compute_derived;
%ignore matchobj_get_index_name;
%include "matchobj.h"
%include "index.h"

PyObject* starkd_search_stars(startree_t* s, double ra, double dec, double radius);

%extend lsst::meas::astrom::net::GlobalAstrometrySolution {
    // Arbitrarily loads the first index.
    PyObject* getIndexStars(double ra, double dec, double radius) {
        $self->loadIndices();
        std::vector<const index_t*> indexList = $self->getIndexList();
        //printf("%i indices\n", indexList.size());
        if (indexList.size() == 0) {
            PyErr_SetString(PyExc_ValueError, "No index files are loaded.");
            return NULL;
        }
        startree_t* skdt = indexList[0]->starkd;
        //MatchObj* match = solver_get_best_match($self->_solver);
        //startree_t* skdt = match->index->starkd;
        return starkd_search_stars(skdt, ra, dec, radius);
    }

    PyObject* getIndexStarsInSolvedField(double pixelmargin) {
        MatchObj* mo = $self->getMatchObject();
        index_t* index = mo->index;
        index_reload(index);
        startree_t* skdt = index->starkd;
        tan_t* tanwcs = &(mo->wcstan);
        return starkd_search_stars_in_field(skdt, tanwcs, pixelmargin);
    }

 };

%template(vectorTagAlongColumn) std::vector<lsst::meas::astrom::net::TagAlongColumn>;

%template(vectorSourceMatch) std::vector<boost::shared_ptr<lsst::afw::detection::SourceMatch> >;

//%template(pair_SourceSet_vector_ints) std::pair<lsst::afw::detection::SourceSet, std::vector<int> >;
