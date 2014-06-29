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
 
#include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"
#include "lsst/afw/image/Wcs.h"

namespace lsst { namespace meas { namespace astrom { namespace sip {

/// \brief Create a list of common objects from a catalogue and an image.
/// 
/// Take a catSet, i.e a SimpleCatalog giving the true positions of objects in
/// ra/dec space, and an imageSet, i.e. a SourceCatalog of the pixel positions
/// of sources detected on a CCD image.  Use the Wcs to calculate the ra/dec of
/// the image source, and then find the subset of objects withcommon positions
/// in the two sets. The returned list of matches is one-to-one, i.e no object
/// in the catalogue is matched twice, and no object in the imageSet is matched
/// twice.
/// 
/// \param catSet  List of catalogue positions with known ra and dec
/// \param imgSet  List of image positions with known x and y
/// \param wcs     A Wcs object to convert from xy to radec
/// \param dist    How close to objects need to be in order to be considered the same
/// 
MatchSrcToCatalogue::MatchSrcToCatalogue(afw::table::SimpleCatalog const& catSet,  
                                         afw::table::SourceCatalog const& imgSet, 
                                         CONST_PTR(lsst::afw::image::Wcs) wcs, 
                                         afw::geom::Angle dist)
{
    setImgSrcSet(imgSet);
    setCatSrcSet(catSet);
    setDist(dist);
    setWcs(wcs);
}    

/// Set a new value for the maximum allowed distance between two matching objects (in ra/dec space) 
void MatchSrcToCatalogue::setDist(afw::geom::Angle dist) {
    if (dist <= 0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "Distance must be > 0");
    }
    _dist = dist;
}

/// Set a different Wcs solution
void MatchSrcToCatalogue::setWcs(CONST_PTR(lsst::afw::image::Wcs) wcs) {
    _wcs = wcs;
}

/// sourceSet is a vector of pointers to Sources.
void MatchSrcToCatalogue::setImgSrcSet(afw::table::SourceCatalog const &srcSet) {
    _imgSet = srcSet;
}

void MatchSrcToCatalogue::setCatSrcSet(afw::table::SimpleCatalog const &srcSet) {
    _catSet = srcSet;
}

void MatchSrcToCatalogue::findMatches() {
    if (!_imgSet.getTable()->hasCentroid()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "SourceTable passed to MatchSrcToCatalogue does not have its centroid slot set."
        );
    }

    for (afw::table::SourceCatalog::const_iterator i = _imgSet.begin(); i != _imgSet.end(); ++i) {
        i->updateCoord(*_wcs);
    }
    
    _match = afw::table::matchRaDec(_catSet, _imgSet, _dist);

    _removeOneToMany();  
    _removeManyToOne();  
}


///We require that out matches be one to one, i.e any element matches no more than once for either 
///the catalogue or the image. However, our implementation of findMatches uses afw::table::matchRaDec()
///which does not garauntee that. This function does the (slow) search and removal.
void MatchSrcToCatalogue::_removeOneToMany() {
    std::size_t size = _match.size();
    for (std::size_t i=0; i<size; ++i) {
        for (std::size_t j=i+1; j<size; ++j) {
            //If the same Source appears twice keep the one with the smaller separation from its match
            if ( _match[i].first == _match[j].first ) {
                //Keep the one with the shorter match distance, and disgard the other
                if (_match[i].distance < _match[j].distance) {
                    _match.erase(_match.begin() + j);
                    size--;
                    j--;
                } else {
                    _match.erase(_match.begin() + i);
                    size--;
                    i--;
                    break;
                }
            }
        }
    }
}


/// This function is identical to
///_removeOneToMany() except first is replaced with second for the match structures
void MatchSrcToCatalogue::_removeManyToOne()  {
    std::size_t size = _match.size();
    for (std::size_t i=0; i<size; ++i) {
        for (std::size_t j=i+1; j<size; ++j) {
            //If the same Source appears twice
            if (_match[i].second == _match[j].second) {
                //Keep the one with the shorter match distance, and disgard the other
                if (_match[i].distance < _match[j].distance) {
                    _match.erase(_match.begin() + j);
                    size--;
                    j--;
                } else {
                    _match.erase(_match.begin() + i);
                    size--;
                    i--;
                    break;
                }
            }
        }
    }
}


afw::table::ReferenceMatchVector MatchSrcToCatalogue::getMatches() {
    findMatches();
    return _match;
}
    
}}}} // namespace lsst::meas::astrom::sip
