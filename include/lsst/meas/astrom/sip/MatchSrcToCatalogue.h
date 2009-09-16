

#ifndef MATCH_SRC_TO_CATALOGUE
#define MATCH_SRC_TO_CATALOGUE

#include <iostream>
#include <cmath>

#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/SourceMatch.h"
#include "lsst/afw/image/Wcs.h"

namespace det = lsst::afw::detection;

namespace lsst { namespace meas { namespace astrom { namespace sip {



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
    
    MatchSrcToCatalogue(const det::SourceSet &catSet,       
                        const det::SourceSet &imgSet,       
                        const lsst::afw::image::Wcs wcs,   
                        double distInArcsec=1       
                       );

    //Copy constructor
    MatchSrcToCatalogue(MatchSrcToCatalogue const &);
    //Destructor
    ~MatchSrcToCatalogue();

    //Mutators
    void setDist(double distInArcsec);
    void setWcs(const lsst::afw::image::Wcs &wcs);
    void setCatSrcSet(const det::SourceSet &srcSet);
    void setImgSrcSet(const det::SourceSet &srcSet);

    void findMatches();

    //Accessors
    std::vector<det::SourceMatch> getMatches();
    

private:
    det::SourceSet _imgSet, _catSet;    ///< Copies of input sets
    std::vector<det::SourceMatch> _match;  ///List of tuples of matching indices
    lsst::afw::image::Wcs _wcs;
    double _distInArcsec;                        ///< How close must two objects be to match 

    det::SourceSet _deepCopySourceSet(const det::SourceSet &in);
};



}}}}

#endif

                

