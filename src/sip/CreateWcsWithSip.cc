// -*- LSST-C++ -*-

#include "lsst/meas/astrom/sip/CreateWcsWithSip.h"


namespace lsst { 
namespace meas { 
namespace astrom { 
namespace sip {


using namespace std;

namespace except = lsst::pex::exceptions;
namespace pexLog = lsst::pex::logging;
namespace afwCoord = lsst::afw::coord;
namespace afwGeom = lsst::afw::geom;
namespace afwImg = lsst::afw::image;
namespace det = lsst::afw::detection;
namespace math = lsst::afw::math;



/// Create a wcs including SIP polynomials of the requested order
/// 
/// \param match A vector of SourceMatches. Each source match consists of two
/// sources (one from the catalogue, one from the image), and the distance
/// between them
/// \param linearWcs A linear WCS that maps pixel position to ra/dec
/// \param order How many terms to compute for the SIP polynomial
CreateWcsWithSip::CreateWcsWithSip(const std::vector<det::SourceMatch> match, 
                                   const afwImg::Wcs &linearWcs, 
                                   int order) :
                                   _matchList(match), 
                                   _linearWcs(linearWcs) {
    
    _createWcs(order);
}


/// Create a wcs including SIP polynomials that reproduces positions to better than a given tolerance.
/// Sip polynomials of progressively higher order are tried until the tolerance is exceeded, or the
/// order of the correction becomes to large
/// 
/// Beware that the SIP can't reproduce positions better than the uncertainty in your centroiding.
/// 
/// \param match A vector of SourceMatches. Each source match consists of two
/// sources (one from the catalogue, one from the image), and the distance
/// between them
/// \param linearWcs A linear WCS that maps pixel position to ra/dec
/// \param maxScatterInArcsec Fit must be better than this value.
/// \param maxOrder If the order of the corrections exceeds this value give up.
CreateWcsWithSip::CreateWcsWithSip(const std::vector<det::SourceMatch> match,
                                   const afwImg::Wcs &linearWcs,
                                   double maxScatterInArcsec,
                                   int maxOrder):
                                   _matchList(match), 
                                   _linearWcs(linearWcs) {

    for (int order = 3; order<maxOrder; ++order) {
        _createWcs(order);
        double scatter = getScatterInArcsec();

        if (scatter < maxScatterInArcsec){
           return;
        }
    }
    
    throw LSST_EXCEPT(except::RuntimeErrorException, "Failed to reach required tolerance"); 
}


///Return a Wcs object including the SIP matrices
afwImg::TanWcs CreateWcsWithSip::getNewWcs() {
    return _newWcs;
}


///Get the scatter in position in pixel space 
double CreateWcsWithSip::getScatterInPixels() {
    unsigned int size = _matchList.size();
    
    vector<double> val;
    val.reserve(size);
    
    for (unsigned int i = 0; i< size; ++i) {
        det::Source::Ptr catSrc = _matchList[i].first;
        det::Source::Ptr imgSrc = _matchList[i].second;

        double imgX = imgSrc->getXAstrom();
        double imgY = imgSrc->getYAstrom();
        
        afwGeom::PointD xy = _newWcs.skyToPixel(catSrc->getRa(), catSrc->getDec());    
        double catX = xy[0];
        double catY = xy[1];
        
        val.push_back(hypot(imgX - catX, imgY - catY));
   }
    
    //Because we're looking at a hypotenuse which is always greater than zero
    //the median is a good metric of the typical size
    math::Statistics stat = math::makeStatistics(val, math::MEDIAN);
    double scatter = stat.getValue(math::MEDIAN);
    
    return scatter;
}
    

///Get the scatter in position in arcseconds
double CreateWcsWithSip::getScatterInArcsec() {
    unsigned int size = _matchList.size();
    
    vector<double> val;
    val.reserve(size);
    
    for (unsigned int i = 0; i< size; ++i) {
        det::Source::Ptr catSrc = _matchList[i].first;
        det::Source::Ptr imgSrc = _matchList[i].second;

        double catRa = catSrc->getRa();
        double catDec = catSrc->getDec();
        
        
        afwCoord::Coord::ConstPtr ad = _newWcs.pixelToSky(imgSrc->getXAstrom(), imgSrc->getYAstrom());    
        double imgRa = ad->getLongitude(afwCoord::DEGREES);
        double imgDec = ad->getLatitude(afwCoord::DEGREES);
        
        //This is not strictly the correct calculation for distance in raDec space,
        //but because we are dealing with distances hopefully << 1" it's a reasonable 
        //approximation
        val.push_back(hypot(imgRa - catRa, imgDec - catDec));
    }
    
    assert(val.size() > 0);
    
    //Because we're looking at a hypotenuse which is always greater than zero
    //the median is a good metric of the typical size
    math::Statistics stat = math::makeStatistics(val, math::MEDIAN);
    double scatter = stat.getValue(math::MEDIAN);
    
    return scatter*3600;    //Convert to arcsec
}



//
// Private functions
//

///Does the donkey work.
void CreateWcsWithSip::_createWcs(int order){

    if (order < 1) {
        string msg = "Order must be greater than or equal to 1";
        throw LSST_EXCEPT(except::RuntimeErrorException, msg); 
    }

    if (_matchList.size() ==0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Match vector is empty"); 
    }

    pexLog::Log mylog(pexLog::Log::getDefaultLog(), "meas.astrom.sip", pexLog::Log::DEBUG);
    mylog.log(pexLog::Log::DEBUG, "Determining Sip parameters");    
                              
    afwGeom::PointD wcsOrigin = _linearWcs.getPixelOrigin();
    
    //Note on naming convention.
    //u,v are as defined in the Shupe et al. (2000) paper as relative pixel coordinates
    //What Shupe calls U,V (undistorted, Linear, pixel coordinates), I call lu and lv. 
    //Having variables u and U seems confusing to me.
    //In a similar style, Shupe uses f,g and F,G to denote forward and reverse distortions,
    //I use f,g and lf, lg
    vector<double> u;   //Relative distorted x position (relative to wcsOrigin)
    vector<double> v;   //Relative y distorted position
    
    vector<double> lu;  //Relative linear x position (relative to wcsOrigin)
    vector<double> lv;  //Relative linear y position

    vector<double> f;  //Difference between image and catalogue position. 
    vector<double> g;  //Difference between image and catalogue position

    vector<double> sf;  //Uncertainty in x position
    vector<double> sg;  //Uncertainty in y position
   
    //Unpack the SourceMatch vector
    mylog.log(pexLog::Log::DEBUG, "Unpacking SourceMatchSet");    
    for (unsigned int i = 0; i< _matchList.size(); ++i) {
    
        det::Source::Ptr catSrc = _matchList[i].first;
        det::Source::Ptr imgSrc = _matchList[i].second;
    
        //Distorted pixel position
        u.push_back(imgSrc->getXAstrom() - wcsOrigin[0]);
        v.push_back(imgSrc->getYAstrom() - wcsOrigin[1]);
        
        //Linear pixel position
        afwGeom::PointD xy = _linearWcs.skyToPixel(catSrc->getRa(), catSrc->getDec());    
        lu.push_back(xy[0] - wcsOrigin[0]);
        lv.push_back(xy[1] - wcsOrigin[1]);
        
        //Forward distortion terms
        f.push_back(xy[0] - imgSrc->getXAstrom());
        g.push_back(xy[1] - imgSrc->getYAstrom());
        
        //Reverse distortion terms are calculated later
        
        //Uncertainties in distortion        
        //The +1s are temporary workarounds guarding against AstromErr's not being set
        sf.push_back(imgSrc->getXAstromErr() + 1);
        sg.push_back(imgSrc->getYAstromErr() + 1);
    }

    //Now fit the forward distortions
    mylog.log(pexLog::Log::DEBUG, "Calculating forward distortion coeffecients");    
    _sipA = _calculateSip(u, v, f, sf, order);
    _sipB = _calculateSip(u, v, g, sg, order);

    //Construct a new wcs from the old one
    mylog.log(pexLog::Log::DEBUG, "Creating new wcs structure");        
    afwGeom::PointD crval = _linearWcs.getSkyOrigin();
    afwGeom::PointD crpix = _linearWcs.getPixelOrigin();
    Eigen::Matrix2d CD = _linearWcs.getCDMatrix();
    
    //The zeroth element of the SIP matrices is just an offset, so the standard calls 
    //for this offset to be folded into a refined value for crpix
    crpix = afwGeom::PointD(crpix - afwGeom::makePointD(_sipA(0, 0), _sipB(0, 0)));
    _sipA(0, 0) = _sipB(0, 0) = 0;
    
    //The A01, A10, B01 and B10 terms in the SIP matrices are just linear corrections, so
    //the standard calls for these to be folded into the CD matrix. The algebra is
    //a little involved here.
    double c00 = CD(0, 0);
    double c01 = CD(0, 1);
    double c10 = CD(1, 0);
    double c11 = CD(1, 1);
    
    CD(0, 0) = c00*(1 + _sipA(1, 0)) + c01*_sipB(1, 0);
    CD(0, 1) = c01*(1 + _sipB(0, 1)) + c00*_sipA(0, 1);
    CD(1, 0) = c10*(1 + _sipA(1, 0)) + c11*_sipB(1, 0);
    CD(1, 1) = c11*(1 + _sipB(0, 1)) + c10*_sipA(0, 1);
    
    _sipA(0, 1) = _sipA(1, 0) = 0;
    _sipB(0, 1) = _sipB(1, 0) = 0;

    afwImg::Wcs tmpWcs = afwImg::Wcs(crval, crpix, CD);

    //Now that we've revised our Wcs, we can now calculate our reverse terms
    mylog.log(pexLog::Log::DEBUG, "Calculating reverse distortion coeffecients");        

    lu.clear();  //Relative linear x position (relative to wcsOrigin)
    lv.clear();  //Relative linear y position

    vector<double> lf;  //Reverse Distortion
    vector<double> lg;  //Reverse Distortion

    wcsOrigin = tmpWcs.getPixelOrigin();
    for (unsigned int i = 0; i< _matchList.size(); ++i) {
        det::Source::Ptr catSrc = _matchList[i].first;
        det::Source::Ptr imgSrc = _matchList[i].second;

        //Linear pixel position
        afwGeom::PointD xy = tmpWcs.skyToPixel(catSrc->getRa(), catSrc->getDec());    
        lu.push_back(xy[0] - wcsOrigin[0]);
        lv.push_back(xy[1] - wcsOrigin[1]);

        //Reverse distortion tersm
        lf.push_back((imgSrc->getXAstrom()-xy[0]));
        lg.push_back(imgSrc->getYAstrom()-xy[1]);
    }
    
    Eigen::MatrixXd sipAp = _calculateSip(lu, lv, lf, sf, order);
    Eigen::MatrixXd sipBp = _calculateSip(lu, lv, lg, sg, order);

    _newWcs = afwImg::TanWcs(crval, crpix, CD, _sipA, _sipB, sipAp, sipBp);
}


/// Calculates the Sip matrices A and B.
/// \param u Relative x position on image ( x - CRPIX1)
/// \param v Relative y position on image ( y - CRPIX2)
/// \param z Amount of distortion at (u,v)
/// \param s variances on z
/// \param order Order of polynomial to fit
/// \param range
/// 
/// For calculating forward distortion polynomials, A,B:
/// u and v are observed pixel positions on the chip, normalised to lie between -1 and 1. 
/// z is the amount of distortion parallel either to u or v. If parallel to u, the Sip A matrix
/// is returned. If parallel to v, the sip B matrix is returned. This parameter is not checked
/// so you have to make sure you are passing in the right thing.
///  
/// To get the reverse distortion polynomials, Ap, Bp, the arguments should be lu, lv, lf (or lg)
/// etc.
Eigen::MatrixXd CreateWcsWithSip::_calculateSip(const std::vector<double> u, const std::vector<double> v,
    const std::vector<double> z, const std::vector<double> s, int order) {

    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(u, v, z, s, order);
    return lsf.getParams();
    
}



/// Needed if we're fitting Chebychev coefficents to calculate the distortion
void CreateWcsWithSip::_scaleVector(std::vector<double>& v) {
    if (v.size() == 0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Trying to rescale an empty vector"); 
    }

    double min = *min_element(v.begin(), v.end());
    double max = *max_element(v.begin(), v.end());

    if (min == max) {
        string msg = "Trying to rescale a vector where all elements have the same value";
        throw LSST_EXCEPT(except::RuntimeErrorException, msg); 
    }
    
    for (unsigned int i = 0; i< v.size(); ++i) {
        v[i] -= min;
        v[i] /= max - min;
        v[i] = 2*v[i] - 1;
    }

}

 
///Convert a 2d matrix of Chebyshev coefficients (as produced by LeastSqFitter2d) into
///a SIP matrix.
Eigen::MatrixXd CreateWcsWithSip::_convertChebyToSip(Eigen::MatrixXd const & cheby) {
    int maxOrder = 5; // Because I've only defined the coeffiencents up to this order
    
    //Coeffs is a lookup table of the polynomial coeffecients of Chebychev functions
    //with order less than maxorder. For example, T2 = 2x^2 - 1, so coeff[2,0] == -1,
    //coeff[2,1] = 0 and coeff[2,2] = 2
    Eigen::MatrixXd coeff(Eigen::MatrixXd::Zero(maxOrder, maxOrder));
    
    
    coeff(0, 0) = 1; //T0
    coeff(1, 1) = 1; //T1
    coeff(2, 0) = -1; coeff(2, 2) = 2;    //T2
    coeff(3, 1) = -3; coeff(3, 3) = 4;
    coeff(4, 0) = 1; coeff(4, 2) = -8, coeff(4, 4) = 8;
    
    int rows = cheby.rows();
    
    if (rows > maxOrder) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Input matrix is too large. Maxorder==5"); 
    }
        
    if (rows != cheby.cols()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Input matrix is not square"); 
    }

    Eigen::MatrixXd out(rows, rows);
    for (int i = 0; i<rows; ++i) {
        for (int j = 0; j<rows; ++j) { //For the (i,j)th elt of the output matrix (x^i y^j
            out(i, j) = 0;
            for (int k = 0; k<rows; ++k) {
                for (int el = 0; el<rows; ++el) {    
                    out(i, j) += cheby(k, el) * coeff(k, i) * coeff(el, j);
                }
            }
        }
    }
    
    return out;
}


    
        
            
        
}}}}
