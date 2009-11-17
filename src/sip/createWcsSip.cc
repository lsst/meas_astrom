

using namespace std;

#include "lsst/meas/astrom/sip/createWcsSip.h"

#include "cpgplot.h"

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
/// \param match. A vector of SourceMatches. Each source match consists of two
/// sources (one from the catalogue, one from the image), and the distance
/// between them
/// \param linearWcs A linear WCS that maps pixel position to ra/dec
/// \param order How many terms to compute for the SIP polynomial
afwImg::Wcs createWcsWithSip(const std::vector<det::SourceMatch> match,
                             const afwImg::Wcs &linearWcs,
                             int order) {
                             
    cpgopen("/xserve");
    
    if(order < 1) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Order must be greater than or equal to 1"); 
    }

    if(match.size() ==0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "match vector is empty"); 
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

    vector<double> lf;  //Reverse Distortion
    vector<double> lg;  //Reverse Distortion

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
        f.push_back(catSrc->getXAstrom() - imgSrc->getXAstrom());
        g.push_back(catSrc->getYAstrom() - imgSrc->getYAstrom());
        
        //Reverse distortion tersm
        lf.push_back(-f[i]);
        lg.push_back(-g[i]);
        
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
    
    //Calculate inverse matrices    
    mylog.log(pexLog::Log::INFO, "Calculating reverse distortion coeffecients");        
    Eigen::MatrixXd sipAp = calculateSip(lu, lv, lf, sf, order);
    Eigen::MatrixXd sipBp = calculateSip(lu, lv, lg, sg, order);
    
    //Construct a new wcs from the old one
    mylog.log(pexLog::Log::INFO, "Creating new wcs structure");        
    afwImg::PointD crval = linearWcs.getOriginRaDec();
    afwImg::PointD crpix = linearWcs.getOriginXY();
    Eigen::Matrix2d CD = linearWcs.getLinearTransformMatrix();
    
    cpgclos();
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
/// To get the reverse distortion polynomials, Ap, Bp, the arguments should be lu, lv, lf (or lg) etc.
Eigen::MatrixXd calculateSip(const vector<double> u, const vector<double> v, const vector<double> z,
                             const vector<double> s, int order) {

    sip::LeastSqFitter2d<math::PolynomialFunction1<double> > lsf(u, v, z, s, order);

    cpgpage();
    cpgswin(-600,600, -1, .5);
    cpgbox("bcnts",0,0,"bcnts",0,0);
    
    for(unsigned int i=0; i<u.size(); ++i) {
        cpgsci(1);
        cpgpt1(u[i], z[i], 2);
        cpgsci(2);
        cpgpt1(u[i], lsf.valueAt(u[i], v[i]), 1);
    }
    cpgsci(1);

    return lsf.getParams();
    
    /*
    sip::LeastSqFitter2d<math::Chebyshev1Function1<double> > lsf(u, v, z, s, order);

    Eigen::MatrixXd cheby = lsf.getParams();
    
    

    cout << "Cheby: " << endl << cheby << endl;
    return convertChebyToSip(cheby);
    */
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
                for(int el=0; el<rows; ++el) {    //Calc effect of (k,el)th element of the input matrix
                    out(i,j) += cheby(k,el) * coeff(k, i) * coeff(el, j);
                }
            }
        }
    }
    
    cout << "SIP: " << endl << out << endl;    
    //These values are explicitly set to zero in the SIP standard, because the linear wcs
    //is supposed to take care of them.
    //For the moment I leave them in, because I think I might need them to update the Wcs
    out(0,0) = out(0,1) = out(1,0) = 0;
    
    return out;
}


double getRmsInPixels(const std::vector<det::SourceMatch> match, const afwImg::Wcs &linearWcs) {
    unsigned int size = match.size();
    
    vector<double> val;
    val.reserve(size);
    
    for(unsigned int i=0; i< size; ++i) {
        det::Source::Ptr imgSrc = match[i].first;
        det::Source::Ptr catSrc = match[i].second;

        double imgX = imgSrc->getXAstrom();
        double imgY = imgSrc->getYAstrom();
        
        afwImg::PointD xy = linearWcs.raDecToXY(catSrc->getRa(), catSrc->getDec());    
        double catX = xy[0];
        double catY = xy[1];
        
        val[i] = hypot(imgX-catX, imgY-catY);
    }
    
    math::Statistics stat = math::makeStatistics(val, math::MEAN | math::STDEV);
    double rms = stat.getValue(math::STDEV);
    
    return rms;
}
    

double getRmsInArcsec(const std::vector<det::SourceMatch> match, const afwImg::Wcs &linearWcs) {
    unsigned int size = match.size();
    
    vector<double> val;
    val.reserve(size);
    
    for(unsigned int i=0; i< size; ++i) {
        det::Source::Ptr imgSrc = match[i].first;
        det::Source::Ptr catSrc = match[i].second;

        double catRa = catSrc->getRa();
        double catDec = catSrc->getDec();
        
        
        afwImg::PointD ad = linearWcs.xyToRaDec(imgSrc->getXAstrom(), imgSrc->getYAstrom());    
        double imgRa = ad[0];
        double imgDec = ad[1];
        
        //This is not strictly the correct calculation for distance in raDec space,
        //but because we are dealing with distances hopefully << 1" it's a reasonable 
        //approximation
        val[i] = hypot(imgRa-catRa, imgDec-catDec);
    }
    
    math::Statistics stat = math::makeStatistics(val, math::MEAN | math::STDEV);
    double rms = stat.getValue(math::STDEV);
    
    return rms;
}
    
        
            
        
}}}}
