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
 
%define sipLib_DOCSTRING
"
Python interface to lsst::afw::meas::astrom::sip classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.astrom.sip",docstring=sipLib_DOCSTRING) sipLib

%{
#include "lsst/pex/logging/FileDestination.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/meas/astrom/sip/LeastSqFitter1d.h"
#include "lsst/meas/astrom/sip/LeastSqFitter2d.h"
#include "lsst/meas/astrom/sip/CreateWcsWithSip.h"
#include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/Match.h"

#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_ASTROM_NUMPY_API
#include "numpy/arrayobject.h"
%}

%init %{ import_array(); %}

%include "lsst/p_lsstSwig.i"
%include "std_string.i"
%include "ndarray.i"

%declareNumPyConverters(Eigen::Matrix2d);
%declareNumPyConverters(Eigen::MatrixXd);

%import "lsst/afw/math/mathLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/table/Source.i"
%import "lsst/afw/table/Match.i"

%lsst_exceptions();

%include "lsst/base.h"
%include "lsst/meas/astrom/sip/LeastSqFitter1d.h"
%include "lsst/meas/astrom/sip/LeastSqFitter2d.h"
%include "lsst/meas/astrom/sip/CreateWcsWithSip.h"
%include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"

%template(LeastSqFitter1dPoly) lsst::meas::astrom::sip::LeastSqFitter1d< lsst::afw::math::PolynomialFunction1<double> >;

%template(LeastSqFitter2dPoly) lsst::meas::astrom::sip::LeastSqFitter2d< lsst::afw::math::PolynomialFunction1<double> >;


%define %declareCreateWcsWithSip(SUFFIX, MATCH)
%template(CreateWcsWithSip##SUFFIX) lsst::meas::astrom::sip::CreateWcsWithSip<MATCH>;
%template(makeCreateWcsWithSip) lsst::meas::astrom::sip::makeCreateWcsWithSip<MATCH>;
%enddef

%declareCreateWcsWithSip(Reference, lsst::afw::table::ReferenceMatch);
%declareCreateWcsWithSip(Source, lsst::afw::table::SourceMatch);
