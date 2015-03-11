// -*- lsst-c++ -*-
#if !defined(LSST_MEAS_ASTROM_MATCHOPTIMISTICB_H)
#define LSST_MEAS_ASTROM_MATCHOPTIMISTICB_H

#include <cmath>
#include <string>
#include <vector>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/Match.h"

namespace lsst {
namespace meas {
namespace astrom {

    /**
    A wrapper around a SimpleRecord or SourceRecord that allows us to record a pixel position
    in a way that is independent of the record type.
    */
    struct RecordProxy {
        PTR(lsst::afw::table::SimpleRecord) record;
        lsst::afw::geom::Point2D position;

        double getX() const { return position.getX(); }
        double getY() const { return position.getY(); }

        /**
        Support implicit conversion to SimpleRecord (discarding the pixel position)
        */
        operator PTR(lsst::afw::table::SimpleRecord)() const { return record; }

        bool operator==(RecordProxy const & other) const { return record == other.record; }
        bool operator!=(RecordProxy const & other) const { return record != other.record; }

        /**
        Construct a RecordProxy

        @param[in] record  the record to wrap
        @param[in] position  pixel position; this is likely to be an initial guess while fitting a WCS
        */
        RecordProxy(
            PTR(lsst::afw::table::SimpleRecord) record,
            lsst::afw::geom::Point2D const & position
        ) : record(record), position(position) {}

        explicit RecordProxy() {}  // default constructor needed so we can call ProxyVector::resize()
    };

    typedef std::vector<RecordProxy> ProxyVector;

    ProxyVector makeProxies(lsst::afw::table::SourceCatalog const & sourceCat);

    ProxyVector makeProxies(lsst::afw::table::SimpleCatalog const & posRefCat);

    struct ProxyPair {
        RecordProxy first;
        RecordProxy second;
        double distance;
        double pa;

        ProxyPair(RecordProxy const & s1, RecordProxy const & s2) : first(s1), second(s2) {
            double x1 = first.position.getX();
            double y1 = first.position.getY();
            double x2 = second.position.getX();
            double y2 = second.position.getY();
            distance = std::hypot(x2-x1, y2-y1);
            pa = std::atan2(y2-y1, x2-x1);
        }
    };

    struct MatchOptimisticBControl {
        LSST_CONTROL_FIELD(refFluxField, std::string, "name of flux field in reference catalog");
        LSST_CONTROL_FIELD(sourceFluxField, std::string, "name of flux field in source catalog");
        LSST_CONTROL_FIELD(numBrightStars, int, "maximum number of bright reference stars to use");
        LSST_CONTROL_FIELD(minMatchedPairs, int, "minimum number of matches");
        LSST_CONTROL_FIELD(matchingAllowancePix, double, "?");
        LSST_CONTROL_FIELD(maxOffsetPix, double, "maximum translation allowed (pixels)");
        LSST_CONTROL_FIELD(maxRotationRad, double, "maximum allowed rotation (rad)");
        LSST_CONTROL_FIELD(angleDiffFrom90, double, "? (deg)");
        LSST_CONTROL_FIELD(numPointsForShape, int, "number of points in a matching shape");
        LSST_CONTROL_FIELD(maxDeterminant, double, "?");

        MatchOptimisticBControl() :
            refFluxField("r_flux"),
            sourceFluxField("slot_ApFlux_flux"),
            numBrightStars(100),
            minMatchedPairs(50),
            matchingAllowancePix(10.0),
            maxOffsetPix(300),
            maxRotationRad(0.02),
            angleDiffFrom90(0.05),
            numPointsForShape(6),
            maxDeterminant(0.02)
        {
            validate();
        }

        void validate() const;

        ~MatchOptimisticBControl() {};
    };
    
    /**
    Match sources to stars in a position reference catalog using optimistic pattern matching B

    Optimistic Pattern Matching is described in V. Tabur 2007, PASA, 24, 189-198
    "Fast Algorithms for Matching CCD Images to a Stellar Catalogue"

    @param[in] posRefCat  catalog of position reference stars
        fields that are used:
        coord
        centroid
        hasCentroid
        control.refFluxField  flux for desired filter
    @param[in] sourceCat  catalog of detected sources
    @param[in] wcs  estimated WCS
    @param[in] control  control object
    @param[in] posRefBegInd  index of first start to use in posRefCat
    @param[in] verbose  true to print diagnostic information to std::cout
    */
    lsst::afw::table::ReferenceMatchVector matchOptimisticB(
        lsst::afw::table::SimpleCatalog const &posRefCat,
        lsst::afw::table::SourceCatalog const &sourceCat,
        MatchOptimisticBControl const &control,
        int posRefBegInd=0,
        bool verbose = false
    );

}}}

#endif
