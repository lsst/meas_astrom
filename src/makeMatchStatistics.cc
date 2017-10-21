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
 
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/meas/astrom/makeMatchStatistics.h"

namespace lsst { 
namespace meas { 
namespace astrom {

template<typename MatchT>
afw::math::Statistics makeMatchStatistics(
    std::vector<MatchT> const & matchList,
    int const flags,  
    afw::math::StatisticsControl const & sctrl
) {
    if (matchList.empty()) {
        throw LSST_EXCEPT(pexExcept::RuntimeError, "matchList is empty");
    }
    std::vector<double> val;
    val.reserve(matchList.size());

    for (auto const & match: matchList) {
        val.push_back(match.distance);
    }
    return afw::math::makeStatistics(val, flags, sctrl);
}

template<typename MatchT>
afw::math::Statistics makeMatchStatisticsInPixels(
    afw::geom::SkyWcs const & wcs,
    std::vector<MatchT> const & matchList,
    int const flags,  
    afw::math::StatisticsControl const & sctrl
) {
    if (matchList.empty()) {
        throw LSST_EXCEPT(pexExcept::RuntimeError, "matchList is empty");
    }
    std::vector<double> val;
    val.reserve(matchList.size());

    for (auto const & match: matchList) {
        auto refPtr = match.first;
        auto srcPtr = match.second;
        auto srcX = srcPtr->getX();
        auto srcY = srcPtr->getY();
        auto refPos = wcs.skyToPixel(refPtr->getCoord());
        auto refX = refPos[0];
        auto refY = refPos[1];
        val.push_back(::hypot(srcX - refX, srcY - refY));
    }
    return afw::math::makeStatistics(val, flags, sctrl);
}

template<typename MatchT>
afw::math::Statistics makeMatchStatisticsInRadians(
    afw::geom::SkyWcs const & wcs,
    std::vector<MatchT> const & matchList,
    int const flags,  
    afw::math::StatisticsControl const & sctrl
) {
    if (matchList.empty()) {
        throw LSST_EXCEPT(pexExcept::RuntimeError, "matchList is empty");
    }
    std::vector<double> val;
    val.reserve(matchList.size());

    for (auto const & match: matchList) {
        auto refPtr = match.first;
        auto srcPtr = match.second;
        auto refCoord = refPtr->getCoord();
        auto srcCoord = wcs.pixelToSky(srcPtr->getCentroid());
        auto angSep = refCoord.angularSeparation(srcCoord);
        val.push_back(angSep.asRadians());
    }
    return afw::math::makeStatistics(val, flags, sctrl);
}

#define INSTANTIATE(MATCH) \
    template afw::math::Statistics makeMatchStatistics<MATCH>( \
        std::vector<MATCH> const & matchList, \
        int const flags,   \
        afw::math::StatisticsControl const & sctrl \
    ); \
    template afw::math::Statistics makeMatchStatisticsInPixels<MATCH>( \
        afw::geom::SkyWcs const & wcs, \
        std::vector<MATCH> const & matchList, \
        int const flags,   \
        afw::math::StatisticsControl const & sctrl \
    ); \
    template afw::math::Statistics makeMatchStatisticsInRadians<MATCH>( \
        afw::geom::SkyWcs const & wcs, \
        std::vector<MATCH> const & matchList, \
        int const flags,   \
        afw::math::StatisticsControl const & sctrl \
    );

INSTANTIATE(afw::table::ReferenceMatch);
INSTANTIATE(afw::table::SourceMatch);

}}}


