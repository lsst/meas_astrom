// -*- LSST-C++ -*-

#include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"


namespace except = lsst::pex::exceptions;
namespace afwImg = lsst::afw::image;
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
sip::MatchSrcToCatalogue::MatchSrcToCatalogue(const det::SourceSet &catSet,  
                                         const det::SourceSet &imgSet, 
                                         const lsst::afw::image::Wcs wcs, 
                                         double distInArcsec  
                                        ){
    setImgSrcSet(imgSet);
    setCatSrcSet(catSet);
    setDist(distInArcsec);
    setWcs(wcs);
}


sip::MatchSrcToCatalogue::~MatchSrcToCatalogue() {}
    

/// Set a new value for the maximum allowed distance between two matching objects (in ra/dec space) 
void sip::MatchSrcToCatalogue::setDist(double distInArcsec){
    if (distInArcsec <= 0){
        throw LSST_EXCEPT(except::InvalidParameterException, "Distance must be > 0");
    }

    _distInArcsec = distInArcsec;
}



/// Set a different Wcs solution
void sip::MatchSrcToCatalogue::setWcs(const lsst::afw::image::Wcs &wcs){
    _wcs = wcs;  //Shallow copy
}

//void MatchSrcToCatalogue::setCatSrcSet(const det::SourceSet &srcSet);

/// Perform a deep copy of a set of sources from the image into this object
///
/// sourceSet is a vector of pointers to Sources. We create deep copies of the objects
/// pointed to so that we can freely mutate them inside the object without affecting
/// the input argument.
void sip::MatchSrcToCatalogue::setImgSrcSet(const det::SourceSet &srcSet) {
    //Destroy the old imgSet
    //imgSet.~imgSet();
    _imgSet = _deepCopySourceSet(srcSet);
}

void sip::MatchSrcToCatalogue::setCatSrcSet(const det::SourceSet &srcSet) {
    _catSet = _deepCopySourceSet(srcSet);
}

void sip::MatchSrcToCatalogue::findMatches() {

    ///@todo Make sure everything is setup as it should be
    //The design of this class ensures all private variables must be set at this point,
    //so assertions should be thrown if this is not the case
    
    //Calculate ra and dec for every imgSrc
    for (unsigned int i = 0; i< _imgSet.size(); ++i) {
        double x = _imgSet[i]->getXAstrom();
        double y = _imgSet[i]->getYAstrom();
        
        afwImg::PointD raDec = _wcs.xyToRaDec(x, y);

        _imgSet[i]->setRa(raDec[0]);
        _imgSet[i]->setDec(raDec[1]);
    }
    
    _match = det::matchRaDec(_catSet, _imgSet, _distInArcsec);

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
    for (unsigned int i = 0; i< size; ++i) {
        for (unsigned int j = i + 1; j< size; ++j) {
            //If the same Source appears twice keep the one with the smaller separation from its match
            if ( _match[i].first == _match[j].first ) {
                //Keep the one with the shorter match distance, and disgard the other
                if ( _match[i].distance < _match[j].distance){
                    _match.erase(_match.begin() + j);
                    size--;
                }
                else {  
                    _match.erase(_match.begin() + i);
                    size--;
                    i--;    //Otherwise the for loop will skip an element
                    j = size + 1; //Nothing else to do for the deleted element
                }
            }
        }
    }
}


/// This function is identical to
///_removeOneToMany() except first is replaced with second for the match structures
void sip::MatchSrcToCatalogue::_removeManyToOne()  {

    
    unsigned int size = _match.size();
    for (unsigned int i = 0; i< size; ++i) {
        for (unsigned int j = i + 1; j< size; ++j) {
            //If the same Source appears twice
            if ( _match[i].second == _match[j].second ) {
                //Keep the one with the shorter match distance, and disgard the other
                if ( _match[i].distance < _match[j].distance ){
                    _match.erase(_match.begin() + j);
                    size--;
                }
                else {  
                    _match.erase(_match.begin() + i);
                    size--;
                    i--;    //Otherwise the for loop will skip an element
                    j = size + 1; //Nothing else to do for the deleted element
                }
            }
        }
    }
}


std::vector<det::SourceMatch> sip::MatchSrcToCatalogue::getMatches() {
    findMatches();
    return _match;
}
    

det::SourceSet sip::MatchSrcToCatalogue::_deepCopySourceSet(const det::SourceSet &in) {

    unsigned int size = in.size();
    det::SourceSet out;

    for (unsigned int i = 0; i<size; ++i){
        //Allocate heap memory for a new object, pointed to by tmp
        det::Source::Ptr tmp(new det::Source);
        //Deep copy the ith Source
        (*tmp) = *(in[i]);
        out.push_back(tmp);
    }

    return out;
}
        
    
