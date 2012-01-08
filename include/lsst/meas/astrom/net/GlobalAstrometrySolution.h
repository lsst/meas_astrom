// -*- LSST-C++ -*-

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
 

#ifndef GLOBAL_ASTROMETRY_SOLUTION_H
#define GLOBAL_ASTROMETRY_SOLUTION_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cfloat>

#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/format.hpp"
#include "Eigen/Core"

#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/SourceMatch.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/utils/Utils.h"

#include "wcslib/wcs.h"

extern "C" {
#include "solver.h"
#include "index.h"
#include "tweak.h"
#include "healpix.h"
#include "bl.h"
#include "log.h"
#include "tic.h"
#include "verify.h"
}

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace net {


///Define variables to indicate tha parity of the image (i.e whether, when the image is rotated so that 
///North is up, East is to the left (Normal), or to the right, (FLIPPED)). The constants PARITY_* are defined
///in solver.h, but not visible to Python. The default setting is UNKNOWN_PARITY. Setting the parity correctly
///reduces the time taken to solve an image by about a factor of two.
enum {
    NORMAL_PARITY = PARITY_NORMAL,
    FLIPPED_PARITY = PARITY_FLIP,
    UNKNOWN_PARITY = PARITY_BOTH
};
    
struct TagAlongColumn_s {
    std::string name;
    tfits_type fitstype;
    std::string ctype;
    std::string units;
    int arraysize;
};
typedef struct TagAlongColumn_s TagAlongColumn;

// Ugly: an internal representation of reference sources, required for fetching
// "tag-along" data
class InternalRefSources;

typedef CONST_PTR(InternalRefSources) InternalRefSourcesCPtr;

struct ReferenceSources_s {
    lsst::afw::detection::SourceSet refsources;
    InternalRefSourcesCPtr intrefsources;
};
typedef struct ReferenceSources_s ReferenceSources;

//!\brief Solve for WCS based only on the positions of stars in an image 
///
///See the examples/ directory for an example of how to use
class GlobalAstrometrySolution {
public:
    typedef boost::shared_ptr<GlobalAstrometrySolution> Ptr;
    typedef boost::shared_ptr<GlobalAstrometrySolution const> ConstPtr;

    //Constructors
    explicit GlobalAstrometrySolution(const std::string policyPath,
                                      lsst::pex::logging::Log mylog=lsst::pex::logging::Log(
                                          lsst::pex::logging::Log::getDefaultLog(),
                                          "meas.astrom.net",
                                          lsst::pex::logging::Log::INFO
                                                                                           )
                                     );

    ~GlobalAstrometrySolution();

    //Tune the solver
    void setDefaultValues();
    void setStarlist(lsst::afw::detection::SourceSet vec);
    void setNumBrightObjects(int N);
    void setImageSize(int W, int H);
    inline void setMinimumImageScale(lsst::afw::geom::Angle scale){ _solver->funits_lower = scale.asArcseconds(); }
    inline void setMaximumImageScale(lsst::afw::geom::Angle scale){ _solver->funits_upper = scale.asArcseconds(); }
    void setImageScaleArcsecPerPixel(lsst::afw::geom::Angle scale);
    void setLogLevel(int level);
    void setMatchThreshold(double threshold);
    void setParity(int parity);

    //Solve for a wcs solution
    bool solve();

    bool solve(lsst::afw::coord::Coord::ConstPtr raDec);
    bool solve(const lsst::afw::geom::Angle ra, const lsst::afw::geom::Angle dec);
    bool solve(const lsst::afw::image::Wcs::Ptr wcsPtr, double imageScaleUncertaintyPercent = 5);

    //Return the solution
    lsst::afw::image::Wcs::Ptr getWcs();
    lsst::afw::image::Wcs::Ptr getDistortedWcs(int order = 3);
    std::vector<lsst::afw::detection::SourceMatch> getMatchedSources(std::string filterName="",
								     std::string idName="");
    lsst::afw::geom::Angle getSolvedImageScale();

    std::vector<std::string> getCatalogueMetadataFields();

    ReferenceSources
    getCatalogueForSolvedField(std::string filter, std::string idname, double margin);

    ReferenceSources
    getCatalogue(lsst::afw::geom::Angle ra, lsst::afw::geom::Angle dec,
                 lsst::afw::geom::Angle radius,
                 std::string filterName, std::string idName,
                 int indexId = -1,
                 bool useIndexHealpix = true,
                 bool resolveDuplicates = true,
                 bool resolveUsingId = true);

    // These are spelled out largely for the benefit of use from python.
    std::vector<double> getTagAlongDouble(InternalRefSourcesCPtr irefs,
                                          std::string columnName);

    std::vector<int> getTagAlongInt(InternalRefSourcesCPtr irefs,
                                    std::string columnName);

    std::vector<boost::int64_t> getTagAlongInt64(InternalRefSourcesCPtr irefs,
                                                 std::string columnName);

    std::vector<bool> getTagAlongBool(InternalRefSourcesCPtr irefs,
                                      std::string columnName);

    std::vector<TagAlongColumn> getTagAlongColumns(int indexId = -1);

    lsst::afw::detection::SourceSet getCatalogue(lsst::afw::geom::Angle radius, std::string filterName,
						 std::string idName);

    // Returns a vector of vectors, where the first two are RA and Dec [in degrees],
    // followed by each of the requested columns.
    std::vector<std::vector<double> > getCatalogueExtra(lsst::afw::geom::Angle ra, lsst::afw::geom::Angle dec,
                                                        lsst::afw::geom::Angle radius,
                                                        std::vector<std::string> columns, int indexId = -1);

    void loadIndices();

    std::vector<const index_t*> getIndexList();

    std::vector<int> getIndexIdList();

    MatchObj* getMatchObject();

    lsst::daf::base::PropertySet::Ptr getMatchedIndexMetadata();

    //Call this before performing a new match
    void reset();



    class RefSourceFilter;
private:
    lsst::pex::logging::Log _mylog;

    std::vector<index_t*> _indexList;

    solver_t *_solver;
    // List of sources
    starxy_t *_starxy;
    // Number of bright sources to use in the solve
    int _numBrightObjects;
        
    // Variables indicating the coordinate system of the solution
    double _equinox;
    std::string _raDecSys;

    // Did we get an astrometric solution?
    bool _isSolved;

    lsst::afw::coord::CoordSystem _getCoordSys();

    index_t *_loadIndexMeta(std::string filename);

    void _solverSetField();
    bool _callSolver(lsst::afw::geom::Angle ra=lsst::afw::geom::NullAngle,
                     lsst::afw::geom::Angle dec=lsst::afw::geom::NullAngle);
    int _addSuitableIndicesToSolver(lsst::afw::geom::Angle minImageSize,
                                    lsst::afw::geom::Angle maxImageSize,
                                    lsst::afw::geom::Angle ra=lsst::afw::geom::NullAngle,
                                    lsst::afw::geom::Angle dec=lsst::afw::geom::NullAngle);

    template <typename T>
    std::vector<T> _getTagAlongData(InternalRefSourcesCPtr irefs,
                                    std::string columnName,
                                    tfits_type ctype);

    index_t* _getIndex(int indexId);

    ReferenceSources
    searchCatalogue(const double* xyzcenter, double r2,
                    std::string filterName, std::string idName,
                    int indexId = -1,
                    bool useIndexHealpix = true,
                    bool resolveDuplicates = true,
                    bool resolveUsingId = true,
                    RefSourceFilter* filt = NULL);


};

}}}}
#endif
