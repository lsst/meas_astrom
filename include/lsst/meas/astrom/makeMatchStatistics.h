// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2015 LSST Corporation.
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

#ifndef MEAS_ASTROM_MAKE_MATCH_STATISTICS_H
#define MEAS_ASTROM_MAKE_MATCH_STATISTICS_H

#include <vector>

#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/table/Match.h"


namespace lsst { 
namespace meas { 
namespace astrom {

/**
 * Compute statistics of the distance field of a match list
 *
 * @param[in] matchList  list of matchList between reference objects and sources; fields read:
 *                  - distance: distance between source and reference object, in arbitrary units;
 *                      the resulting statistics have the same units as distance
 * @param[in] flags  what to calculate; OR constants such as lsst::afw::math::MEAN, MEANCLIP, STDDEV, MEDIAN,
 *                  defined in lsst/afw/math/Statitics.h's Property enum
 * @param[in] sctrl  statistics configuration
 */
template<typename MatchT>
afw::math::Statistics makeMatchStatistics(
    std::vector<MatchT> const & matchList,
    int const flags,  
    afw::math::StatisticsControl const & sctrl = afw::math::StatisticsControl()
);

/**
 * Compute statistics of on-detector radial separation for a match list, in pixels
 *
 * @param[in] wcs  WCS describing pixel to sky transformation
 * @param[in] matchList  list of matchList between reference objects and sources; fields read:
 *                  - first: reference object; only the coord is read
 *                  - second: source; only the centroid is read
 * @param[in] flags  what to calculate; OR constants such as lsst::afw::math::MEAN, MEANCLIP, STDDEV, MEDIAN,
 *                  defined in lsst/afw/math/Statitics.h's Property enum
 * @param[in] sctrl  statistics configuration
 */
template<typename MatchT>
afw::math::Statistics makeMatchStatisticsInPixels(
    afw::geom::SkyWcs const & wcs,
    std::vector<MatchT> const & matchList,
    int const flags,  
    afw::math::StatisticsControl const & sctrl = afw::math::StatisticsControl()
);

/**
 * Compute statistics of on-sky radial separation for a match list, in radians
 *
 * @param[in] wcs  WCS describing pixel to sky transformation
 * @param[in] matchList  list of matchList between reference objects and sources; fields read:
 *                  - first: reference object; only the coord is read
 *                  - second: source; only the centroid is read
 * @param[in] flags  what to calculate; OR constants such as lsst::afw::math::MEAN, MEANCLIP, STDDEV, MEDIAN,
 *                  defined in lsst/afw/math/Statitics.h's Property enum
 * @param[in] sctrl  statistics configuration
 */
template<typename MatchT>
afw::math::Statistics makeMatchStatisticsInRadians(
    afw::geom::SkyWcs const & wcs,
    std::vector<MatchT> const & matchList,
    int const flags,
    afw::math::StatisticsControl const & sctrl = afw::math::StatisticsControl()
);

}}}

#endif
