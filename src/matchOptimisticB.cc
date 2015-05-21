#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include "boost/scoped_array.hpp"
#include "boost/shared_array.hpp"
#include "boost/multi_index_container.hpp"
#include "boost/multi_index/sequenced_index.hpp"
#include "boost/multi_index/ordered_index.hpp"
#include "boost/multi_index/global_fun.hpp"

#include "gsl/gsl_linalg.h"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/meas/astrom/matchOptimisticB.h"

namespace pexExcept = lsst::pex::exceptions;
namespace afwTable = lsst::afw::table;

namespace {
    using namespace lsst::meas::astrom;

    afwTable::RecordId getRefId(ProxyPair const &proxyPair) {
        return proxyPair.first.record->getId();
    }

    // a collection of ReferenceMatch indexed by insertion order and by reference object ID
    struct refIdTag{}; // tag for retrieval by reference object ID
    struct insertionOrderTag{}; // tag for retrieval by insertion order
    typedef boost::multi_index_container<
        ProxyPair,
        boost::multi_index::indexed_by<         // index by insertion order, so brightest sources first
            boost::multi_index::sequenced<
                boost::multi_index::tag<insertionOrderTag>
            >,
            boost::multi_index::ordered_unique< // index by reference object ID
                boost::multi_index::tag<refIdTag>,
                boost::multi_index::global_fun<
                    ProxyPair const &, afwTable::RecordId , &getRefId
            > >
        >
    > MultiIndexedProxyPairList;

    // Algorithm is based on V.Tabur 2007, PASA, 24, 189-198
    // "Fast Algorithms for Matching CCD Images to a Stellar Catalogue"

    /**
    Return |ang1-ang2| wrapped into the range [0, 2 pi]

    @param[in] ang1  angle 1 (rad)
    @param[in] ang2  angle 2 (rad)
    */
    inline double absDeltaAngle(double ang1, double ang2) {
        return std::fmod(std::fabs(ang1 - ang2), M_PI*2);
    }

    bool cmpPair(ProxyPair const &a, ProxyPair const &b) {
        return a.distance > b.distance;
    }

    /**
    Return true if a has brighter flux than b

    When used as a compare function for std::sort, the result is sorted by decreasing brightness
    */
    struct CompareProxyFlux {

        bool operator()(RecordProxy const & a, RecordProxy const & b) const {
            double aFlux = a.record->get(key);
            double bFlux = b.record->get(key);
            if (lsst::utils::isnan(aFlux)) {
                aFlux = 0.0;
            }
            if (lsst::utils::isnan(bFlux)) { 
                bFlux = 0.0;
            }
            return aFlux > bFlux;
        }

        afwTable::Key<double> key;
    };

    /**
    Sort a copy of a vector of ProxyVector by decreasing brightness, and return a subset

    @param[in] a  list of ProxyVector
    @param[in] key  key of value on which to sort, in descending order
    @param[in] num  maximum number of elements to return
    @param[in] startInd  starting index

    @throw pexExcept::InvalidParameterError if startInd is out of range
    */
    ProxyVector selectPoint(
        ProxyVector const &a,
        afwTable::Key<double> const & key,
        std::size_t num,
        std::size_t startInd=0
    ) {
        if (startInd >= a.size()) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "startInd too big");
        }
        CompareProxyFlux cmp = {key};
        ProxyVector b(a);
        std::sort(b.begin(), b.end(), cmp);
        std::size_t const endInd = std::min(startInd + num, b.size());
        return ProxyVector(b.begin() + startInd, b.begin() + endInd);
    }

    std::vector<ProxyPair> searchPair(
        std::vector<ProxyPair> const &a,
        ProxyPair const &p,
        double e,
        double e_dpa
    ) {
        std::vector<ProxyPair> v;

        for (size_t i = 0; i < a.size(); i++) {
            double dd = std::fabs(a[i].distance - p.distance);
            double dpa = absDeltaAngle(a[i].pa, p.pa);
            if (dd < e && dpa < e_dpa) {
                v.push_back(a[i]);
            }
        }

        return v;
    }

    std::vector<ProxyPair>::iterator searchPair3(
        std::vector<ProxyPair> &a,
        ProxyPair const &p,
        ProxyPair const &q,
        double e,
        double dpa,
        double e_dpa = 0.02
    ) {
        std::vector<ProxyPair>::iterator idx = a.end();
        double dd_min = 1.E+10;
        //double dpa_min = e_dpa;

        for (auto i = a.begin(); i != a.end(); ++i) {
            double dd = std::fabs(i->distance - p.distance);
    #if 1
            if (dd < e &&
                absDeltaAngle(p.pa, i->pa - dpa) < e_dpa &&
                dd < dd_min &&
                (i->first == q.first)) {
                dd_min = dd;
                idx = i;
            }
    #else
            if (dd < e &&
                absDeltaAngle(p.pa, i->pa - dpa) < dpa_min) {
                dpa_min = std::fabs(p.pa - i->pa - dpa);
                idx = i;
            }
    #endif
        }

        return idx;
    }

    void transform(
        int order,
        boost::shared_array<double> const & coeff,
        double x,
        double y,
        double *xn,
        double *yn
    ) {
        int ncoeff = (order + 1) * (order + 2) / 2;
        *xn = 0.0;
        *yn = 0.0;
        int n = 0;
        for (int i = 0; i <= order; i++) {
            for (int k = 0; k <= i; k++) {
                int j = i - k;
                *xn += coeff[n] * pow(x, j) * pow(y, k);
                *yn += coeff[n+ncoeff] * pow(x, j) * pow(y, k);
                n++;
            }
        }
    }

    boost::shared_array<double> polyfit(
        int order,
        ProxyVector const &img,
        ProxyVector const &posRefCat
    ) {
        int ncoeff = (order + 1) * (order + 2) / 2;
        boost::scoped_array<int> xorder(new int[ncoeff]);
        boost::scoped_array<int> yorder(new int[ncoeff]);

        int n = 0;
        for (int i = 0; i <= order; i++) {
            for (int k = 0; k <= i; k++) {
                int j = i - k;
                xorder[n] = j;
                yorder[n] = k;
                n++;
            }
        }

        boost::scoped_array<int> flag(new int[img.size()]);
        for (size_t k = 0; k < img.size(); k++) {
            flag[k] = 1;
        }

        boost::scoped_array<double> a_data(new double[ncoeff*ncoeff]);
        boost::scoped_array<double> b_data(new double[ncoeff]);
        boost::scoped_array<double> c_data(new double[ncoeff]);

        boost::shared_array<double> coeff(new double[ncoeff*2]);

        for (int loop = 0; loop < 1; loop++) {
            for (int i = 0; i < ncoeff; i++) {
                for (int j = 0; j < ncoeff; j++) {
                    a_data[i*ncoeff+j] = 0.0;
                    for (size_t k = 0; k < img.size(); k++) {
                        if (flag[k] == 1) {
                            a_data[i*ncoeff+j] += pow(img[k].getX(), xorder[i]) * 
                                pow(img[k].getY(), yorder[i]) * 
                                pow(img[k].getX(), xorder[j]) * 
                                pow(img[k].getY(), yorder[j]);
                        }
                    }
                }
                b_data[i] = c_data[i] = 0.0;
                for (unsigned int k = 0; k < img.size(); k++) {
                    if (flag[k] == 1) {
                        b_data[i] += pow(img[k].getX(), xorder[i]) * 
                            pow(img[k].getY(), yorder[i]) * 
                            posRefCat[k].getX();
                        c_data[i] += pow(img[k].getX(), xorder[i]) * 
                            pow(img[k].getY(), yorder[i]) * 
                            posRefCat[k].getY();
                    }
                }
            }

            gsl_matrix_view a = gsl_matrix_view_array(a_data.get(), ncoeff, ncoeff);
            gsl_vector_view b = gsl_vector_view_array(b_data.get(), ncoeff);
            gsl_vector_view c = gsl_vector_view_array(c_data.get(), ncoeff);

            boost::shared_ptr<gsl_vector> x(gsl_vector_alloc(ncoeff), gsl_vector_free);
            boost::shared_ptr<gsl_vector> y(gsl_vector_alloc(ncoeff), gsl_vector_free);

            int s;

            boost::shared_ptr<gsl_permutation> p(gsl_permutation_alloc(ncoeff), gsl_permutation_free);

            gsl_linalg_LU_decomp(&a.matrix, p.get(), &s);
            gsl_linalg_LU_solve(&a.matrix, p.get(), &b.vector, x.get());
            gsl_linalg_LU_solve(&a.matrix, p.get(), &c.vector, y.get());

            for (int i = 0; i < ncoeff; i++) {
                coeff[i] = x->data[i];
                coeff[i+ncoeff] = y->data[i];
            }

            double S, Sx, Sy, Sxx, Syy;
            S = Sx = Sy = Sxx = Syy = 0.0;
            for (size_t k = 0; k < img.size(); k++) {
                if (flag[k] == 1) {
                    double x0 = img[k].getX();
                    double y0 = img[k].getY();
                    double x1, y1;
                    transform(order, coeff, x0, y0, &x1, &y1);
                    S   += 1.;
                    Sx  += (x1 - posRefCat[k].getX());
                    Sxx += (x1 - posRefCat[k].getX()) * (x1 - posRefCat[k].getX());
                    Sy  += (y1 - posRefCat[k].getY());
                    Syy += (y1 - posRefCat[k].getY()) * (y1 - posRefCat[k].getY());
                }
            }
            double x_sig = std::sqrt((Sxx - Sx * Sx / S) / S);
            double y_sig = std::sqrt((Syy - Sy * Sy / S) / S);
            //std::cout << x_sig << " " << y_sig << std::endl;
    
            for (size_t k = 0; k < img.size(); k++) {
                double x0 = img[k].getX();
                double y0 = img[k].getY();
                double x1, y1;
                transform(order, coeff, x0, y0, &x1, &y1);
                if (std::fabs(x1-posRefCat[k].getX()) > 2. * x_sig ||
                    std::fabs(y1-posRefCat[k].getY()) > 2. * y_sig) {
                    flag[k] = 0;
                }
            }

        }

        return coeff;
    }

    /**
    Return the reference object nearest the given source

    @param[in] posRefCat  list of reference object ProxyRecords;
        read: x, y
    @param[in] x  source pixel position in x
    @param[in] y  source pixel position in y
    @param[in] matchingAllowancePix  maximum allowed distance between reference object and source (pixels);
    @return a pair of:
    - posRefCat iterator; if no match then posRefCat.end()
    - distance between reference object and source (in pixels); approx. matchingAllowancePix if no match
    */
    std::pair<ProxyVector::const_iterator, double> searchNearestPoint(
        ProxyVector const &posRefCat,
        double x,
        double y,
        double matchingAllowancePix
    ) {
        auto minDistSq = matchingAllowancePix*matchingAllowancePix;
        auto foundPtr = posRefCat.end();
        for (auto posRefPtr = posRefCat.begin(); posRefPtr != posRefCat.end(); ++posRefPtr) {
            auto const dx = posRefPtr->getX() - x;
            auto const dy = posRefPtr->getY() - y;
            auto const distSq = dx*dx + dy*dy;
            if (distSq < minDistSq) {
                foundPtr = posRefPtr;
                minDistSq = distSq;
            }
        }
        return std::make_pair(foundPtr, std::sqrt(minDistSq));
    }

    /**
    Find nearest reference object to a source and add the match to a list

    Callers should match sources from brightest to faintest, because once
    a reference object has been used in a match it will not be accepted again
    as a match for another source (and the later source will have no match).

    It is easy to change this code to handle duplicate reference objects by keeping
    the better match (the one that has a smaller distance between source and reference
    object). However, that raises the danger of false matches in the situation that
    the user provides a source catalog with many faint stars, some of which may be spurious.

    @param[in,out] proxyPairList  list of reference object, source matches
    @param[in] coeff  array of TAN WCS coefficients (I think)
    @param[in] posRefCat  list of position reference object proxies, sorted by decreasing brightness
    @param[in] source  source to match
    @param[in] matchingAllowancePix  maximum allowed distance between reference object and source (pixels)
    */
    void addNearestMatch(
        MultiIndexedProxyPairList &proxyPairList,
        boost::shared_array<double> coeff,
        ProxyVector const & posRefCat,
        RecordProxy const &source,
        double matchingAllowancePix
    ) {
        double x1, y1;
        auto x0 = source.getX();
        auto y0 = source.getY();
        transform(1, coeff, x0, y0, &x1, &y1);
        auto refObjDist = searchNearestPoint(posRefCat, x1, y1, matchingAllowancePix);
        if (refObjDist.first == posRefCat.end()) {
            // no reference object sufficiently close; do nothing
            return;
        }
        auto existingMatch = proxyPairList.get<refIdTag>().find(refObjDist.first->record->getId()); 
        if (existingMatch == proxyPairList.get<refIdTag>().end()) {
            // reference object not used before; add new entry
            auto proxyPair = ProxyPair(*refObjDist.first, source);
            proxyPairList.get<refIdTag>().insert(proxyPair);
            return;
#if 0
        } else {
            // This code keeps the closest match. It is disabled because it appears to be unnecessary
            // and may be undesirable, since the image may have very faint sources which may give
            // false matches. Note that the calling algorithm checks for matches starting with the
            // brightest source and works down in source brightness.
            if (existingMatch->distance <= refObjDist.second) {
                // reference object used before and old source is closer; do nothing
                return;
            } else {
                // reference object used before, but new source is closer; delete old entry and add new
                proxyPairList.get<refIdTag>().erase(existingMatch);
                auto proxyPair = ProxyPair(*refObjDist.first, source);
                proxyPairList.get<refIdTag>().insert(proxyPair);
                return;
            }
#endif
        }
    }

    /**
    Perform a verification pass on a possible match

    @param[in] coeff  array of TAN WCS coefficients?
    @param[in] posRefCat  list of position reference object proxies, sorted by decreasing brightness
    @param[in] sourceCat  list of sources proxies, sorted by decreasing brightness
    @param[in] matchingAllowancePix  maximum allowed distance between reference object and source (pixels)
    */
    afwTable::ReferenceMatchVector FinalVerify(
        boost::shared_array<double> coeff,
        ProxyVector const & posRefCat,
        ProxyVector const & sourceCat,
        double matchingAllowancePix,
        bool verbose
    ) {
        MultiIndexedProxyPairList proxyPairList;

        for (auto sourcePtr = sourceCat.begin(); sourcePtr != sourceCat.end(); ++sourcePtr) {
            addNearestMatch(proxyPairList, coeff, posRefCat, *sourcePtr, matchingAllowancePix);
        }
        int order = 1;
        if (proxyPairList.size() > 5) {
            for (auto j = 0; j < 100; j++) {
                auto prevNumMatches = proxyPairList.size();
                auto srcMat = ProxyVector();
                auto catMat = ProxyVector();
                srcMat.reserve(proxyPairList.size());
                catMat.reserve(proxyPairList.size());
                for (auto matchPtr = proxyPairList.get<refIdTag>().begin();
                    matchPtr != proxyPairList.get<refIdTag>().end(); ++matchPtr) {
                    catMat.push_back(matchPtr->first);
                    srcMat.push_back(matchPtr->second);
                }
                coeff = polyfit(order, srcMat, catMat);
                proxyPairList.clear();

                for (ProxyVector::const_iterator sourcePtr = sourceCat.begin();
                    sourcePtr != sourceCat.end(); ++sourcePtr) {
                    addNearestMatch(proxyPairList, coeff, posRefCat, *sourcePtr, matchingAllowancePix);
                }
                if (proxyPairList.size() == prevNumMatches) {
                    break;
                }
            }
        }
        auto matPair = afwTable::ReferenceMatchVector();
        matPair.reserve(proxyPairList.size());
        for (auto proxyPairIter = proxyPairList.get<insertionOrderTag>().begin();
            proxyPairIter != proxyPairList.get<insertionOrderTag>().end(); ++proxyPairIter) {
            matPair.push_back(afwTable::ReferenceMatch(
                proxyPairIter->first.record,
                boost::static_pointer_cast<afwTable::SourceRecord>(proxyPairIter->second.record),
                proxyPairIter->distance
            ));
        }

        return matPair;
    }

} // anonymous namespace

namespace lsst {
namespace meas {
namespace astrom {

    void MatchOptimisticBControl::validate() const {
        if (refFluxField.empty()) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "refFluxField must be specified");
        }
        if (sourceFluxField.empty()) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "sourceFluxField must be specified");
        }
        if (numBrightStars <= 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "numBrightStars must be positive");
        }
        if (minMatchedPairs < 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "minMatchedPairs must not be negative");
        }
        if (matchingAllowancePix <= 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "matchingAllowancePix must be positive");
        }
        if (maxOffsetPix <= 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "maxOffsetPix must be positive");
        }
        if (maxRotationDeg <= 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "maxRotationRad must be positive");
        }
        if (allowedNonperpDeg <= 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "allowedNonperpDeg must be positive");
        }
        if (numPointsForShape <= 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "numPointsForShape must be positive");
        }
        if (maxDeterminant <= 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "maxDeterminant must be positive");
        }
    }

    ProxyVector makeProxies(afwTable::SourceCatalog const & sourceCat) {
        ProxyVector r;
        r.reserve(sourceCat.size());
        for (afwTable::SourceCatalog::const_iterator sourcePtr = sourceCat.begin();
            sourcePtr != sourceCat.end(); ++sourcePtr) {
            r.push_back(RecordProxy(sourcePtr, sourcePtr->getCentroid()));
        }
        return r;
    }

    ProxyVector makeProxies(afwTable::SimpleCatalog const & posRefCat) {
        auto centroidKey = posRefCat.getSchema().find<afwTable::Point<double>>("centroid").key;
        auto hasCentroidKey = posRefCat.getSchema().find<afwTable::Flag>("hasCentroid").key;
        ProxyVector r;
        r.reserve(posRefCat.size());
        for (afwTable::SimpleCatalog::const_iterator posRefPtr = posRefCat.begin();
            posRefPtr != posRefCat.end(); ++posRefPtr) {
            if (!posRefPtr->get(hasCentroidKey)) {
                throw LSST_EXCEPT(pexExcept::InvalidParameterError, "unknown centroid");
            }
            r.push_back(RecordProxy(posRefPtr, posRefPtr->get(centroidKey)));
        }
        return r;
    }

    afwTable::ReferenceMatchVector matchOptimisticB(
        afwTable::SimpleCatalog const &posRefCat,
        afwTable::SourceCatalog const &sourceCat,
        MatchOptimisticBControl const &control,
        int posRefBegInd,
        bool verbose
    ) {
        control.validate();
        if (posRefBegInd < 0) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "posRefBegInd < 0");
        }
        if (posRefBegInd >= posRefCat.size()) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "posRefBegInd too big");
        }
        double const maxRotationRad = afw::geom::degToRad(control.maxRotationDeg);
        double const allowedNonperpRad = afw::geom::degToRad(control.allowedNonperpDeg);

        ProxyVector posRefProxyCat = makeProxies(posRefCat);
        ProxyVector sourceProxyCat = makeProxies(sourceCat);

        // sourceSubCat contains at most the numBrightStars brightest sources, sorted by decreasing flux
        ProxyVector sourceSubCat = selectPoint(
            sourceProxyCat,
            sourceCat.getSchema().find<double>(control.sourceFluxField).key,
            control.numBrightStars);

        // posRefSubCat skips the initial posRefBegInd brightest reference objects and contains
        // at most the next len(sourceSubCat) + 25 brightest reference objects, sorted by decreasing flux
        ProxyVector posRefSubCat = selectPoint(
            posRefProxyCat,
            posRefCat.getSchema().find<double>(control.refFluxField).key,
            sourceSubCat.size()+25,
            posRefBegInd);
        if (verbose) {
            std::cout << "Catalog sizes: " << sourceSubCat.size() << " " << posRefSubCat.size() << std::endl;
        }

        // Construct a list of pairs of position reference stars sorted by increasing separation
        std::vector<ProxyPair> posRefPairList;
        size_t const posRefCatSubSize = posRefSubCat.size();
        for (size_t i = 0; i < posRefCatSubSize-1; i++) {
            for (size_t j = i+1; j < posRefCatSubSize; j++) {
                posRefPairList.push_back(ProxyPair(posRefSubCat[i], posRefSubCat[j]));
            }
        }
        std::sort(posRefPairList.begin(), posRefPairList.end(), cmpPair);

        // Construct a list of pairs of sources sorted by increasing separation
        std::vector<ProxyPair> sourcePairList;
        size_t const sourceSubCatSize = sourceSubCat.size();
        for (size_t i = 0; i < sourceSubCatSize-1; i++) {
            for (size_t j = i+1; j < sourceSubCatSize; j++) {
                sourcePairList.push_back(ProxyPair(sourceSubCat[i], sourceSubCat[j]));
            }
        }
        std::sort(sourcePairList.begin(), sourcePairList.end(), cmpPair);

        afwTable::ReferenceMatchVector matPair;
        afwTable::ReferenceMatchVector matPairSave;
        std::vector<afwTable::ReferenceMatchVector> matPairCand;

        size_t const fullShapeSize = control.numPointsForShape - 1;   // Max size of shape array
        for (size_t ii = 0; ii < sourcePairList.size(); ii++) {
            ProxyPair p = sourcePairList[ii];

            std::vector<ProxyPair> q = searchPair(posRefPairList, p,
                control.matchingAllowancePix, maxRotationRad);

            // If candidate pairs are found
            if (q.size() != 0) {

                std::vector<ProxyPair> srcMatPair;
                std::vector<ProxyPair> catMatPair;

                // Go through candidate pairs
                for (size_t l = 0; l < q.size(); l++) {
                    // sign matters, so don't use deltaAng; no need to wrap because
                    // the result is used with deltaAng later
                    double dpa = p.pa - q[l].pa;

                    srcMatPair.clear();
                    catMatPair.clear();

                    srcMatPair.push_back(p);
                    catMatPair.push_back(q[l]);

                    if (verbose) {
                        std::cout << "p dist: " << p.distance << " pa: " << p.pa << std::endl;
                        std::cout << "q dist: " << q[l].distance << " pa: " << q[l].pa << std::endl;
                    }

                    for (size_t k = 0; k < sourceSubCat.size(); k++) {
                        if (p.first == sourceSubCat[k] || p.second == sourceSubCat[k]) continue;

                        ProxyPair pp(p.first, sourceSubCat[k]);
                
                        std::vector<ProxyPair>::iterator r = searchPair3(posRefPairList, pp, q[l],
                            control.matchingAllowancePix, dpa, maxRotationRad);
                        if (r != posRefPairList.end()) {
                            srcMatPair.push_back(pp);
                            catMatPair.push_back(*r);
                            if (verbose) {
                                std::cout << "  p dist: " << pp.distance << " pa: " << pp.pa << std::endl;
                                std::cout << "  r dist: " << (*r).distance << " pa: " << (*r).pa << std::endl;
                            }
                            if (srcMatPair.size() == fullShapeSize) {
                                break;
                            }
                        }
                    }

                    bool goodMatch = false;
                    if (srcMatPair.size() == fullShapeSize) {
                        goodMatch = true;
                        for (size_t k = 1; k < catMatPair.size(); k++) {
                            if (catMatPair[0].first != catMatPair[k].first) {
                                goodMatch = false;
                                break;
                            }
                        }
                    }

                    if (goodMatch && srcMatPair.size() == fullShapeSize) {

                        ProxyVector srcMat;
                        ProxyVector catMat;

                        srcMat.push_back(srcMatPair[0].first);
                        catMat.push_back(catMatPair[0].first);
                        for (size_t k = 0; k < srcMatPair.size(); k++) {
                            srcMat.push_back(srcMatPair[k].second);
                            catMat.push_back(catMatPair[k].second);
                        }

                        boost::shared_array<double> coeff = polyfit(1, srcMat, catMat);

                        if (verbose) {
                            for (size_t k = 0; k < srcMat.size(); k++) {
                                std::cout << "circle(" << srcMat[k].getX() << ","
                                          << srcMat[k].getY() << ",10) # color=green" << std::endl;
                                std::cout << "circle(" << catMat[k].getX() << ","
                                          << catMat[k].getY() << ",10) # color=red" << std::endl;
                                std::cout << "line(" << srcMat[0].getX() << "," << srcMat[0].getY() << ","
                                          << srcMat[k].getX() << "," << srcMat[k].getY()
                                          << ") # line=0 0 color=green" << std::endl;
                                std::cout << "line(" << catMat[0].getX() << "," << catMat[0].getY() << ","
                                          << catMat[k].getX() << "," << catMat[k].getY()
                                          << ") # line=0 0 color=red" << std::endl;
                            }
                        }

                        double a = coeff[1];
                        double b = coeff[2];
                        double c = coeff[4];
                        double d = coeff[5];
                        double theta = std::acos((a*b+c*d)/(std::sqrt(a*a+c*c)*std::sqrt(b*b+d*d))) /
                            M_PI * 180.0;
                        if (verbose) {
                            std::cout << "Linear fit from match:" << std::endl;
                            std::cout << coeff[0] << " " << coeff[1] << " " << coeff[2] << std::endl;
                            std::cout << coeff[3] << " " << coeff[4] << " " << coeff[5] << std::endl;
                            std::cout << coeff[1] * coeff[5] - coeff[2] * coeff[4] - 1. << std::endl;
                            std::cout << theta << std::endl;
                        }
                        if (std::fabs(coeff[1] * coeff[5] - coeff[2] * coeff[4] - 1.) 
                                > control.maxDeterminant ||
                            std::fabs(theta - 90.) > allowedNonperpRad ||
                            std::fabs(coeff[0]) > control.maxOffsetPix ||
                            std::fabs(coeff[3]) > control.maxOffsetPix) {
                            if (verbose) {
                                std::cout << "Bad; continuing" << std::endl;
                            }
                            continue;
                        } else {
                            double x0, y0, x1, y1;
                            int num = 0;
                            srcMat.clear();
                            catMat.clear();
                            for (size_t i = 0; i < sourceSubCat.size(); i++) {
                                x0 = sourceSubCat[i].getX();
                                y0 = sourceSubCat[i].getY();
                                transform(1, coeff, x0, y0, &x1, &y1);
                                auto refObjDist = searchNearestPoint(posRefSubCat, x1, y1,
                                    control.matchingAllowancePix);
                                if (refObjDist.first != posRefSubCat.end()) {
                                    num++;
                                    srcMat.push_back(sourceSubCat[i]);
                                    catMat.push_back(*refObjDist.first);
                                    if (verbose) {
                                        std::cout << "Match: " << x0 << "," << y0 << " --> "
                                            << x1 << "," << y1 << " <==> "
                                            << refObjDist.first->getX() << "," << refObjDist.first->getY()
                                            << std::endl;
                                    }
                                }
                            }
                            if (num <= control.numPointsForShape) {
                                // Can get matrix = 0,0,0,0; everything matches a single catalog object
                                if (verbose) {
                                    std::cout << "Insufficient initial matches; continuing" << std::endl;
                                }
                                continue;
                            }
                            coeff = polyfit(1, srcMat, catMat);
                            if (verbose) {
                                std::cout << "Coefficients from initial matching:" << std::endl;
                                for (size_t i = 0; i < 6; ++i) {
                                    std::cout << coeff[i] << " ";
                                }
                                std::cout << std::endl;
                            }

                            matPair = FinalVerify(coeff, posRefProxyCat, sourceProxyCat,
                                control.matchingAllowancePix, verbose);
                            if (verbose) {
                                std::cout << "Number of matches: " << matPair.size() << " vs " <<
                                    control.minMatchedPairs << std::endl;
                            }
                            if (matPair.size() <= control.minMatchedPairs) {
                                if (verbose) {
                                    std::cout << "Insufficient final matches; continuing" << std::endl;
                                }
                                if (matPair.size() > matPairSave.size()) {
                                    matPairSave = matPair;
                                }
                                continue;
                            } else {
                                if (verbose) {
                                    std::cout << "Finish" << std::endl;
                                }
                                matPairCand.push_back(matPair);
                                if (matPairCand.size() == 3) {
                                    goto END;
                                }
                            }
                        }
                    }
                }
            }
        }

     END:
        if (matPairCand.size() == 0) {
            return matPairSave;
        } else {
            size_t nmatch = matPairCand[0].size();
            afwTable::ReferenceMatchVector matPairRet = matPairCand[0];
            for (size_t i = 1; i < matPairCand.size(); i++) {
                if (matPairCand[i].size() > nmatch) {
                    nmatch = matPairCand[i].size();
                    matPairRet = matPairCand[i];
                }
            }
            return matPairRet;
        }
    }

}}} // namespace lsst::afw::astrom
