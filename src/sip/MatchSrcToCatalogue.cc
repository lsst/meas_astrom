

#include "lsst/meas/astrom/sip/MatchSrcToCatalogue.h"


namespace except = lsst::pex::exceptions;
namespace afwImg = lsst::afw::image;
using namespace lsst::meas::astrom::sip;

MatchSrcToCatalogue::MatchSrcToCatalogue(const det::SourceSet &catSet,  ///<Input list of objects from catalogue
                                         const det::SourceSet &imgSet, ///<Input list of objects from image
                                         const lsst::afw::image::Wcs wcs,   ///< World Coordinate System object
                                         double distInArcsec  
                                        )
{
    setImgSrcSet(imgSet);
    setCatSrcSet(catSet);
    setDist(distInArcsec);
    setWcs(wcs);

    findMatches();
}


MatchSrcToCatalogue::~MatchSrcToCatalogue() {
}
    

/// Set a new value for the maximum allowed distance between two matching objects (in ra/dec space) 
void MatchSrcToCatalogue::setDist(double distInArcsec)
{
    if(distInArcsec <= 0){
        throw LSST_EXCEPT(except::InvalidParameterException, "Distance must be > 0");
    }

    _distInArcsec = distInArcsec;
}



/// Set a different Wcs solution
void MatchSrcToCatalogue::setWcs(const lsst::afw::image::Wcs &wcs){
    _wcs = wcs;  //Shallow copy
}

//void MatchSrcToCatalogue::setCatSrcSet(const det::SourceSet &srcSet);

/// Perform a deep copy of a set of sources from the image into this object
///
/// sourceSet is a vector of pointers to Sources. We create deep copies of the objects
/// pointed to so that we can freely mutate them inside the object without affecting
/// the input argument.
void MatchSrcToCatalogue::setImgSrcSet(const det::SourceSet &srcSet) {
    //Destroy the old imgSet
    //imgSet.~imgSet();
    _imgSet = _deepCopySourceSet(srcSet);
}

void MatchSrcToCatalogue::setCatSrcSet(const det::SourceSet &srcSet) {
    _catSet = _deepCopySourceSet(srcSet);
}

void MatchSrcToCatalogue::findMatches() {

    ///@todo Make sure everything is setup as it should be
    //The design of this class ensures all private variables must be set at this point,
    //so assertions should be thrown if this is not the case
    
    //Calculate ra and dec for every imgSrc
    for(unsigned int i=0; i< _imgSet.size(); ++i) {
        double x = _imgSet[i]->getXAstrom();
        double y = _imgSet[i]->getYAstrom();
        
        afwImg::PointD raDec = _wcs.xyToRaDec(x, y);

        _imgSet[i]->setRa(raDec[0]);
        _imgSet[i]->setDec(raDec[1]);
    }
    
    _match = det::matchRaDec(_imgSet, _catSet, _distInArcsec);

    if(_match.size() == 0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "No matching objects found");
    }
}


std::vector<det::SourceMatch> MatchSrcToCatalogue::getMatches() {
    return _match;
}
    

det::SourceSet MatchSrcToCatalogue::_deepCopySourceSet(const det::SourceSet &in) {

    unsigned int size = in.size();
    det::SourceSet out;

    for(unsigned int i=0; i<size; ++i){
        //Allocate heap memory for a new object, pointed to by tmp
        det::Source::Ptr tmp(new det::Source);
        //Deep copy the ith Source
        (*tmp) = *(in[i]);
        out.push_back(tmp);
    }

    return out;
}
        
    
