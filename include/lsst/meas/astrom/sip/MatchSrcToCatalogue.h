// -*- LSST-C++ -*-

#ifndef MATCH_SRC_TO_CATALOGUE
#define MATCH_SRC_TO_CATALOGUE

#include <iostream>
#include <cmath>

#include "lsst/base.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/SourceMatch.h"

namespace lsst {
    namespace afw {
        namespace image {
            class Wcs;
        }
    }
namespace meas { 
namespace astrom { 
namespace sip {

namespace det = lsst::afw::detection;


/// Match a SourceSet of objects with known ra/dec with a SourceSet of objects with known xy positions
/// Take a catalogue of objects with known positions, a catalogue of objects with known xy, and a wcs
/// to convert one xy <--> ra/dec. This class then finds the set of objects with common ra/dec.
///
/// The simplest usage is to create an object of this class, then access the corresponence sets with
/// getMatchedImgSet() and getMatchedCatSet(). Creating the object automatically calculates the sets
/// of corresponences for you. If you are unhappy with these matches, you can change one or more of your
/// input arguments and redo the match with findMatches()
class MatchSrcToCatalogue{
public:

    typedef boost::shared_ptr<MatchSrcToCatalogue> Ptr;
    typedef boost::shared_ptr<MatchSrcToCatalogue const> ConstPtr;
    
    MatchSrcToCatalogue(det::SourceSet const& catSet,
                        det::SourceSet const& imgSet,       
                        CONST_PTR(lsst::afw::image::Wcs) wcs,   
                        double distInArcsec = 1.0     
                       );

    //Don't need a Copy constructor
    
    //Destructor
    ~MatchSrcToCatalogue();

    //Mutators
    void setDist(double distInArcsec);
    void setWcs(CONST_PTR(lsst::afw::image::Wcs) wcs);
    void setCatSrcSet(const det::SourceSet &srcSet);
    void setImgSrcSet(const det::SourceSet &srcSet);

    void findMatches();

    //Accessors
    std::vector<det::SourceMatch> getMatches();
    

private:
    det::SourceSet _imgSet, _catSet;      ///< Copies of input sets
    std::vector<det::SourceMatch> _match; ///List of tuples of matching indices
    CONST_PTR(lsst::afw::image::Wcs) _wcs;
    double _distInArcsec;               ///< How close must two objects be to match 

    det::SourceSet _deepCopySourceSet(const det::SourceSet &in);
    void _removeOneToMany();
    void _removeManyToOne();
};



}}}}

#endif

                

