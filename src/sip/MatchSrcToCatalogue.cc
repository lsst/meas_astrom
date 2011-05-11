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

namespace except = lsst::pex::exceptions;
namespace afwCoord = lsst::afw::coord;
namespace afwGeom = lsst::afw::geom;
namespace det = lsst::afw::detection;
namespace sip = lsst::meas::astrom::sip;


/// \brief Create a list of common objects from a catalogue and an image.
/// 
/// Take a catalogueSet, i.e a SourceSet of Sources giving the true positions of objects in ra/dec space,
/// and an imageSet, i.e. the pixel positions of sources detected on a CCD image. Use the Wcs to 
/// calculate the ra/dec of the image source, and then find the subset of objects with common positions
/// in the two sets. The returned list of matches is one-to-one, i.e no object in the catalogue is matched
/// twice, and no object in the imageSet is matched twice.
/// 
/// \param catSet  List of catalogue positions with known ra and dec
/// \param imgSet  List of image positions with known x and y
/// \param wcs     A Wcs object to convert from xy to radec
/// \param distInArcsec How close to objects need to be in order to be considered the same
/// 
sip::MatchSrcToCatalogue::MatchSrcToCatalogue(det::SourceSet const& catSet,  
                                              det::SourceSet const& imgSet, 
                                              CONST_PTR(lsst::afw::image::Wcs) wcs, 
                                              afwGeom::Angle dist) {
    setImgSrcSet(imgSet);
    setCatSrcSet(catSet);
    setDist(dist);
    setWcs(wcs);
}


sip::MatchSrcToCatalogue::~MatchSrcToCatalogue() {}
    

/// Set a new value for the maximum allowed distance between two matching objects (in ra/dec space) 
void sip::MatchSrcToCatalogue::setDist(afwGeom::Angle dist) {
    if (dist <= 0) {
        throw LSST_EXCEPT(except::InvalidParameterException, "Distance must be > 0");
    }
    _dist = dist;
}

/// Set a different Wcs solution
void sip::MatchSrcToCatalogue::setWcs(CONST_PTR(lsst::afw::image::Wcs) wcs) {
    _wcs = wcs;
}

/// sourceSet is a vector of pointers to Sources.
void sip::MatchSrcToCatalogue::setImgSrcSet(const det::SourceSet &srcSet) {
    _imgSet = srcSet;
}

void sip::MatchSrcToCatalogue::setCatSrcSet(const det::SourceSet &srcSet) {
    _catSet = srcSet;
}

void sip::MatchSrcToCatalogue::findMatches() {
    //The design of this class ensures all private variables must be set at this point,
    //so assertions should be thrown if this is not the case

    //printf("findMatches: %i img and %i srces\n", (int)(_imgSet.size()), (int)(_catSet.size()));
    //Calculate ra and dec for every imgSrc
    for (unsigned int i = 0; i < _imgSet.size(); ++i) {
        /*
         printf("img source %i: X,Y Astrom: %.1f, %.1f\n", i,
         _imgSet[i]->getXAstrom(), _imgSet[i]->getYAstrom());
         */
        // set the RA,Dec{,Astrom} fields from {X,Y}Astrom
        _imgSet[i]->setRaDecFromXy(_wcs);
        //printf("  -> RA,Dec %.3f, %.3f\n", _imgSet[i]->getRa().asDegrees(), _imgSet[i]->getDec().asDegrees());
    }

    //For completeness, set x and y for the catSrc
    for (unsigned int i = 0; i < _catSet.size(); ++i) {
        _catSet[i]->setRaDecAstrom(_catSet[i]->getRaDec());
        _catSet[i]->setXyAstromFromRaDec(_wcs);
        /*
         printf("cat source %i: RA,Dec %.3f, %.3f\n", i, _catSet[i]->getRa().asDegrees(), _catSet[i]->getDec().asDegrees());
         printf("  RA,DecAstrom: (%.3f, %.3f)\n", _catSet[i]->getRaAstrom().asDegrees(), _catSet[i]->getDecAstrom().asDegrees());
         printf("  -> X,YAstrom %.1f, %.1f\n", _catSet[i]->getXAstrom(), _catSet[i]->getYAstrom());
         */
    }
    
    _match = det::matchRaDec(_catSet, _imgSet, _dist);

    _removeOneToMany();  
    _removeManyToOne();  
   
    if (_match.size() == 0) {
        std::cout << _imgSet.size() << " " << _catSet.size() << std::endl;
        throw LSST_EXCEPT(except::RuntimeErrorException, "No matching objects found");
    }
}


///We require that out matches be one to one, i.e any element matches no more than once for either 
///the catalogue or the image. However, our implementation of findMatches uses det::matchRaDec()
///which does not garauntee that. This function does the (slow) search and removal.
void sip::MatchSrcToCatalogue::_removeOneToMany() {
    unsigned int size = _match.size();
    for (unsigned int i=0; i<size; ++i) {
        for (unsigned int j=i+1; j<size; ++j) {
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
void sip::MatchSrcToCatalogue::_removeManyToOne()  {
    unsigned int size = _match.size();
    for (unsigned int i=0; i<size; ++i) {
        for (unsigned int j=i+1; j<size; ++j) {
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


std::vector<det::SourceMatch> sip::MatchSrcToCatalogue::getMatches() {
    findMatches();
    return _match;
}
    

