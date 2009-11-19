

using namespace std;

#include "lsst/meas/astrom/sip/createWcsSip.h"


namespace except = lsst::pex::exceptions;
namespace pexLog = lsst::pex::logging;
namespace afwImg = lsst::afw::image;
namespace det = lsst::afw::detection;
namespace math = lsst::afw::math;


namespace lsst { namespace meas { namespace astrom { namespace sip {


/// Create a wcs including SIP polynomials
/// 
/// Given a list of matching sources between a catalogue and an image,
/// and a linear Wcs that describes the mapping from pixel space in the image
/// and ra/dec space in the catalogue, calculate discrepancies between the two
/// and compute SIP distortion polynomials to describe the discrepancy
///
/// Note that the SIP standard insists* that the lowest three terms in the distortion
/// polynomials be zero (A00, A10, A01, B00, etc.). To achieve this, we need to 
/// adjust the values of CD and CRPIX from the input wcs. This may not be the 
/// behaviour you expect.
///
/// *The standard is detailed in Shupe et al. (2005, ASP Conf.), and this 
/// fact is mentioned very obliquly between Eqns 3 and 4. Many other implementations
/// do not have this requirement.
/// 
/// \param match. A vector of SourceMatches. Each source match consists of two
/// sources (one from the catalogue, one from the image), and the distance
/// between them
/// \param linearWcs A linear WCS that maps pixel position to ra/dec
/// \param order How many terms to compute for the SIP polynomial
afwImg::Wcs createWcsWithSip(const std::vector<det::SourceMatch> match,
                             const afwImg::Wcs &linearWcs,
                             int order) {
    
    if(order < 1) {
        string msg = "Order must be greater than or equal to 1";
        throw LSST_EXCEPT(except::RuntimeErrorException, msg); 
    }

    if(match.size() ==0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Match vector is empty"); 
    }

    pexLog::Log mylog(pexLog::Log::getDefaultLog(), "meas.astrom.sip", pexLog::Log::DEBUG);
    mylog.log(pexLog::Log::INFO, "Determining Sip parameters");    
                              
    afwImg::PointD wcsOrigin = linearWcs.getOriginXY();
    
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
    mylog.log(pexLog::Log::INFO, "Unpacking SourceMatchSet");    
    for(unsigned int i=0; i< match.size(); ++i) {
    
        det::Source::Ptr catSrc = match[i].first;
        det::Source::Ptr imgSrc = match[i].second;
    
        //Distorted pixel position
        u.push_back(imgSrc->getXAstrom() - wcsOrigin[0]);
        v.push_back(imgSrc->getYAstrom() - wcsOrigin[1]);
        
        //Linear pixel position
        afwImg::PointD xy = linearWcs.raDecToXY(catSrc->getRa(), catSrc->getDec());    
        lu.push_back(xy[0] - wcsOrigin[0]);
        lv.push_back(xy[1] - wcsOrigin[1]);
        
        //Forward distortion terms
        f.push_back(xy[0] - imgSrc->getXAstrom());
        g.push_back(xy[1] - imgSrc->getYAstrom());
        
        //Reverse distortion terms are calculated later
        
        //Uncertainties in distortion        
        //The +1s are temporary workarounds guarding against AstromErr's not being set
        sf.push_back(imgSrc->getXAstromErr()+1);
        sg.push_back(imgSrc->getYAstromErr()+1);
    }
    
    /*
    //Scale vectors to range from [-1,1] (so we can fit a Cheby to them)
    scaleVector(u);
    scaleVector(v);
    scaleVector(lu);
    scaleVector(lv);
    //scaleVector(f);
    // scaleVector(g);
    */

    //Now fit the forward distortions
    mylog.log(pexLog::Log::INFO, "Calculating forward distortion coeffecients");    
    Eigen::MatrixXd sipA = calculateSip(u, v, f, sf, order);
    Eigen::MatrixXd sipB = calculateSip(u, v, g, sg, order);

    //Construct a new wcs from the old one
    mylog.log(pexLog::Log::INFO, "Creating new wcs structure");        
    afwImg::PointD crval = linearWcs.getOriginRaDec();
    afwImg::PointD crpix = linearWcs.getOriginXY();
    Eigen::Matrix2d CD = linearWcs.getLinearTransformMatrix();
    
    //The zeroth element of the SIP matrices is just an offset, so the standard calls 
    //for this offset to be folded into a refined value for crpix
    crpix = crpix - afwImg::PointD(sipA(0,0), sipB(0,0));
    sipA(0,0) = sipB(0,0) = 0;
    
    //The A01, A10, B01 and B10 terms in the SIP matrices are just linear corrections, so
    //the standard calls for these to be folded into the CD matrix. The algebra is
    //a little involved here.
    double c00 = CD(0,0);
    double c01 = CD(0,1);
    double c10 = CD(1,0);
    double c11 = CD(1,1);
    
    CD(0,0) = c00*(1+sipA(1,0)) + c01*sipB(1,0);
    CD(0,1) = c01*(1+sipB(0,1)) + c00*sipA(0,1);
    CD(1,0) = c10*(1+sipA(1,0)) + c11*sipB(1,0);
    CD(1,1) = c11*(1+sipB(0,1)) + c10*sipA(0,1);
    
    sipA(0,1) = sipA(1,0) = 0;
    sipB(0,1) = sipB(1,0) = 0;

    afwImg::Wcs tmpWcs = afwImg::Wcs(crval, crpix, CD);

    //Now that we've revised our Wcs, we can now calculate our reverse terms
    mylog.log(pexLog::Log::INFO, "Calculating reverse distortion coeffecients");        

    lu.clear();  //Relative linear x position (relative to wcsOrigin)
    lv.clear();  //Relative linear y position

    vector<double> lf;  //Reverse Distortion
    vector<double> lg;  //Reverse Distortion

    wcsOrigin = tmpWcs.getOriginXY();
    for(unsigned int i=0; i< match.size(); ++i) {
        det::Source::Ptr catSrc = match[i].first;
        det::Source::Ptr imgSrc = match[i].second;

        //Linear pixel position
        afwImg::PointD xy = tmpWcs.raDecToXY(catSrc->getRa(), catSrc->getDec());    
        lu.push_back(xy[0] - wcsOrigin[0]);
        lv.push_back(xy[1] - wcsOrigin[1]);

        //Reverse distortion tersm
        lf.push_back((imgSrc->getXAstrom()-xy[0]));
        lg.push_back(imgSrc->getYAstrom()-xy[1]);
    }
    
    Eigen::MatrixXd sipAp = calculateSip(lu, lv, lf, sf, order);
    Eigen::MatrixXd sipBp = calculateSip(lu, lv, lg, sg, order);

    return afwImg::Wcs(crval, crpix, CD, sipA, sipB, sipAp, sipBp);
}



void scaleVector(vector<double>& v) {
    if(v.size() == 0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Trying to rescale an empty vector"); 
    }

    double min = *min_element(v.begin(), v.end());
    double max = *max_element(v.begin(), v.end());

    if(min == max) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Trying to rescale a vector where all elements have the same value"); 
    }
    
    for(unsigned int i=0; i< v.size(); ++i) {
        v[i] -= min;
        v[i] /= max-min;
        v[i] = 2*v[i] - 1;
    }

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
Eigen::MatrixXd calculateSip(const vector<double> u, const vector<double> v, const vector<double> z,
                             const vector<double> s, int order) {

    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(u, v, z, s, order);
    return lsf.getParams();
    
}



 
///Convert a 2d matrix of Chebyshev coefficients (as produced by LeastSqFitter2d) into
///a SIP matrix.
Eigen::MatrixXd convertChebyToSip(Eigen::MatrixXd cheby) {
    int maxOrder=5;
    
    //Coeffs is a lookup table of the polynomial coeffecients of Chebychev functions
    //with order less than maxorder. For example, T2 = 2x^2 - 1, so coeff[2,0] == -1,
    //coeff[2,1] = 0 and coeff[2,2] = 2
    Eigen::MatrixXd coeff(Eigen::MatrixXd::Zero(maxOrder, maxOrder));
    
    
    coeff(0,0) = 1; //T0
    coeff(1,1) = 1; //T1
    coeff(2,0) = -1; coeff(2,2) = 2;    //T2
    coeff(3,1) = -3; coeff(3,3) = 4;
    coeff(4,0) = 1; coeff(4,2) = -8, coeff(4,4)=8;
    
    int rows = cheby.rows();
    
    if(rows > maxOrder) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Input matrix is too large. Maxorder==5"); 
    }
        
    if(rows != cheby.cols()) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Input matrix is not square"); 
    }

    Eigen::MatrixXd out(rows, rows);
    for(int i=0; i<rows; ++i) {
        for(int j=0; j<rows; ++j) { //For the (i,j)th elt of the output matrix (x^i y^j
            out(i,j) = 0;
            for(int k=0; k<rows; ++k) {
                for(int el=0; el<rows; ++el) {    
                    out(i,j) += cheby(k,el) * coeff(k, i) * coeff(el, j);
                }
            }
        }
    }
    
    return out;
}


double getScatterInPixels(const std::vector<det::SourceMatch> match, const afwImg::Wcs &wcs) {
    unsigned int size = match.size();
    
    vector<double> val;
    val.reserve(size);
    
    for(unsigned int i=0; i< size; ++i) {
        det::Source::Ptr catSrc = match[i].first;
        det::Source::Ptr imgSrc = match[i].second;

        double imgX = imgSrc->getXAstrom();
        double imgY = imgSrc->getYAstrom();
        
        afwImg::PointD xy = wcs.raDecToXY(catSrc->getRa(), catSrc->getDec());    
        double catX = xy[0];
        double catY = xy[1];
        
        val.push_back(hypot(imgX-catX, imgY-catY));
   }
    
    //Because we're looking at a hypotenuse which is always greater than zero
    //the median is a good metric of the typical size
    math::Statistics stat = math::makeStatistics(val, math::MEDIAN);
    double scatter = stat.getValue(math::MEDIAN);
    
    return scatter;
}
    

double getScatterInArcsec(const std::vector<det::SourceMatch> match, const afwImg::Wcs &wcs) {
    unsigned int size = match.size();
    
    vector<double> val;
    val.reserve(size);
    
    for(unsigned int i=0; i< size; ++i) {
        det::Source::Ptr catSrc = match[i].first;
        det::Source::Ptr imgSrc = match[i].second;

        double catRa = catSrc->getRa();
        double catDec = catSrc->getDec();
        
        
        afwImg::PointD ad = wcs.xyToRaDec(imgSrc->getXAstrom(), imgSrc->getYAstrom());    
        double imgRa = ad[0];
        double imgDec = ad[1];
        
        //This is not strictly the correct calculation for distance in raDec space,
        //but because we are dealing with distances hopefully << 1" it's a reasonable 
        //approximation
        val.push_back(hypot(imgRa-catRa, imgDec-catDec));
    }
    
    assert(val.size() > 0);
    
    //Because we're looking at a hypotenuse which is always greater than zero
    //the median is a good metric of the typical size
    math::Statistics stat = math::makeStatistics(val, math::MEDIAN);
    double rms = stat.getValue(math::MEDIAN);
    
    return rms*3600;    //Convert to arcsec
}
    
        
            
        
}}}}
