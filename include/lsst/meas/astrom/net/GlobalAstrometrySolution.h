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
#include "Eigen/Core.h"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/SourceMatch.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
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
    

//When deciding whether an index file needs to be loaded from disk, we we want to use
//position on the sky as a factor (if we're doing a blind search we don't). 
//Used by _addSuitableIndicesToSolver. If ra or dec are set to this value, then
//the function does not use position as a factor
enum { NO_POSITION_SET=-360};


struct TagAlongColumn_s {
    std::string name;
    tfits_type fitstype;
    std::string ctype;
    std::string units;
    int arraysize;
};
typedef struct TagAlongColumn_s TagAlongColumn;

struct ReferenceSources_s {
    lsst::afw::detection::SourceSet refsources;
    int indexid;
    std::vector<int> inds;
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
    inline void setMinimumImageScale(double scale){   _solver->funits_lower = scale; }
    inline void setMaximumImageScale(double scale){   _solver->funits_upper = scale; }
    void setImageScaleArcsecPerPixel(double scale);
    void setLogLevel(int level);
    void setMatchThreshold(double threshold);
    void setParity(int parity);

    ///Finding a match requires a minimum number of objects in the field. astrometry.net
    ///recommends at least 20, which is the default value for the class. Setting
    ///the value any lower is probably not a good idea, and may lead to false matches.
    void setMinimumNumberOfObjectsToAccept(double num){
        _minimumNumberOfObjectsToAccept = num;
    }

    //Solve for a wcs solution
    bool solve();
    bool solve(const afw::image::PointD raDec);
    bool solve(double ra, double dec);
    bool solve(const lsst::afw::image::Wcs::Ptr wcsPtr, double imageScaleUncertaintyPercent = 5);

    //Return the solution
    lsst::afw::image::Wcs::Ptr getWcs();
    lsst::afw::image::Wcs::Ptr getDistortedWcs(int order = 3);
    std::vector<lsst::afw::detection::SourceMatch> getMatchedSources(std::string filterName="",
								     std::string idName="");
    double getSolvedImageScale();

    std::vector<std::string> getCatalogueMetadataFields();

    ReferenceSources
    getCatalogueForSolvedField(std::string filter, std::string idname, double margin);

    ReferenceSources
    getCatalogue(double ra, double dec, double radiusInArcsec, 
                 std::string filterName, std::string idName,
                 int indexId = -1);

    std::vector<double> getTagAlongDouble(int indexId, std::string columnName,
                                          std::vector<int> inds);
    std::vector<int> getTagAlongInt(int indexId, std::string columnName,
                                    std::vector<int> inds);
    std::vector<boost::int64_t> getTagAlongInt64(int indexId, std::string columnName,
                                                 std::vector<int> inds);
    std::vector<bool> getTagAlongBool(int indexId, std::string columnName,
                                      std::vector<int> inds);

    std::vector<TagAlongColumn> getTagAlongColumns(int indexId = -1);


    lsst::afw::detection::SourceSet getCatalogue(double radiusInArcsec, std::string filterName,
						 std::string idName);

    std::vector<std::vector<double> > getCatalogueExtra(double ra, double dec, double radiusInArcsec,
                                                       std::vector<std::string> columns, int indexId = -1);
    //const startree_t* getStarTree(int indexId);

    void loadIndices();

    std::vector<const index_t*> getIndexList();

    std::vector<int> getIndexIdList();

    MatchObj* getMatchObject();

    lsst::daf::base::PropertySet::Ptr getMatchedIndexMetadata();

    //Call this before performing a new match
    void reset();



private:
    lsst::pex::logging::Log _mylog;

    std::vector<index_t*> _indexList;

    solver_t *_solver;
    starxy_t *_starxy;   ///List of sources to solve for
    int _numBrightObjects;   //Only use the brightest objects in solve.
    //Refuse to add starlists with fewer objects than this value
    int _minimumNumberOfObjectsToAccept; 
        
    //Variables indicating the coordinate system of the solution
    double _equinox;
    std::string _raDecSys;
    
    bool _isSolved;

    index_t *_loadIndexMeta(std::string filename);

    void _solverSetField();
    bool _callSolver(double ra=NO_POSITION_SET, double dec=NO_POSITION_SET);
    int _addSuitableIndicesToSolver(double minImageSizeArcsec, double maxImageSizeArcsec, \
        double ra=NO_POSITION_SET, double dec=NO_POSITION_SET);    

    template <typename T>
    std::vector<T> _getTagAlongData(int indexId, std::string columnName,
                                    tfits_type ctype, std::vector<int> inds);

    index_t* _getIndex(int indexId);

};

}}}}
#endif
