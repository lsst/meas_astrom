// -*- lsst-C++ -*-

// Astrometry.net include files...
extern "C" {
#include "astrometry/solver.h"
#include "astrometry/index.h"
#include "astrometry/multiindex.h"
#include "astrometry/starkd.h"
#include "astrometry/fitsioutils.h"
#include "astrometry/fitstable.h"
#include "astrometry/log.h"
#include "astrometry/tic.h"
#include "astrometry/healpix.h"

#undef ATTRIB_FORMAT
#undef FALSE
#undef TRUE

#undef logdebug
#undef debug
}

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "boost/format.hpp"

#include "lsst/utils/pybind11.h"
#include "lsst/meas/astrom/detail/utils.h"
#include "lsst/base.h"
#include "lsst/log/Log.h"
#include "lsst/daf/persistence.h"
#include "lsst/daf/base/Persistable.h"
#include "lsst/daf/base/PropertyList.h"
#include "lsst/afw/table.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/coord.h"

namespace lsst {
namespace meas {
namespace astrom {

/*
 * A thin C++ wrapper around astrometry_net's multiindex_t struct.
 *
 * This provide memory management and a few methods used by LSST.
 */
class MultiIndex {

public:

    /**
     * Construct a MultiIndex from an astrometry.net multi-index file
     */
    MultiIndex(std::string const & filepath);

    ~MultiIndex() {}

    /**
     * Get the specified index
     */
    index_t * operator[](int i) const {
        return multiindex_get(_multiindex.get(), i);
    }

    /**
     * Add an index read from a file
     */
    void addIndex(std::string const & filepath, bool metadataOnly);

    /**
     * Is this multi-index in range of the specified cone?
     */
    int isWithinRange(double ra, double dec, double radius_deg);

    /**
     * Unload the indices
     */
    void unload() {
        multiindex_unload(_multiindex.get());
    }

    std::string getName() const {
        return _multiindex->fits->filename;
    }

    /**
     * Get the number of indices
     */
    int getLength() const {
        return multiindex_n(_multiindex.get());
    }

    /**
     * Reload the indices
     */
    void reload();

private:
    
    struct _Deleter {
        void operator()(multiindex_t* m) {
            multiindex_free(m);
        }
    };

    std::unique_ptr<multiindex_t, _Deleter> _multiindex;
};


/**
 * A thin C++ wrapper around astrometry.net's solver_t struct.
 *
 * This provide memory management and methods used by LSST.
 */
class Solver {

public:

    explicit Solver();

    ~Solver();

    /**
    Load reference objects in a region of the sky described by a center coordinate and a radius

    @param[in] inds  list of star kd-trees from astrometry.net
    @param[in] ctrCoord  center of search region
    @param[in] radius  search radius
    @param[in] idCol  name of ID column in astrometry.net data
    @param[in] filterNameList  names of filters in astrometry.net data
    @param[in] magColList  names of magnitude columns in astrometry.net data
    @param[in] magErrColList  names of magnitude uncertainty (sigma) columns in astrometry.net data
    @param[in] starGalCol  name of "starGal" column (true if object is a star) in astrometry.net data
    @param[in] varCol  name of "var" column (true if brightness is variable) in astrometry.net data
    @param[in] uniqueIds  if true then only return unique IDs (the first of each seen)

    Returned schema:
    - id
    - coord: sky position (an lsst::afw::coord::IcrsCoord)
    - centroid: centroid on some exposure, if relevant (an lsst::afw::geom::Point2D);
        returned value is not set
    - hasCentroid: if true then centroid has been set; returned value is false
    - <filterName>_flux: flux in the specified filter (double)
    - <filterName>_fluxSigma: flux uncertainty in the specified filter (double)
    - resolved (if starGalCol specified): true if object is not resolved
    - variable (if varCol specified): true if brightness is variable
    - photometric: true if not resolved (or starGalCol blank) and not variable (or varCol blank);
        note that if starGalCol and varCol both blank then all objects are claimed to be photometric
    */
    lsst::afw::table::SimpleCatalog getCatalog(
        std::vector<index_t*> inds,
        lsst::afw::coord::Coord const &ctrCoord,
        lsst::afw::geom::Angle const &radius,
        const char* idCol,
        std::vector<std::string> const& filterNameList,
        std::vector<std::string> const& magColList,
        std::vector<std::string> const& magErrColList,
        const char* starGalCol,
        const char* varCol,
        bool uniqueIds=true);

    std::shared_ptr<lsst::daf::base::PropertyList> getSolveStats() const;

    std::shared_ptr<lsst::afw::image::Wcs> getWcs();

    bool didSolve() const {
        return solver_did_solve(_solver.get());
    }

    void run(double cpulimit);

    std::pair<double, double> getQuadSizeRangeArcsec() const {
        double qlo,qhi;
        solver_get_quad_size_range_arcsec(_solver.get(), &qlo, &qhi);
        return std::make_pair(qlo, qhi);
    }

    /**
     * Add indices to the solver
     *
     * The indices are bare pointers whose memory is managed by the caller.
     * Typically the indices are owned by a MultiIndex object owned by the caller.
     */
    void addIndices(std::vector<index_t*> inds);

    /**
     * Set parity to flipped (if true) or normal (if false)
     */
    void setParity(bool flipped) {
        _solver->parity = flipped ? PARITY_FLIP : PARITY_NORMAL;
    }

    void setMatchThreshold(double threshold) {
        solver_set_keep_logodds(_solver.get(), threshold);
    }

    void setPixelScaleRange(double low, double high) {
        _solver->funits_lower = low;
        _solver->funits_upper = high;
    }

    void setRaDecRadius(double ra, double dec, double radius_deg) {
        solver_set_radec(_solver.get(), ra, dec, radius_deg);
    }

    void setImageSize(int width, int height);

    void setMaxStars(int maxStars) {
        _solver->endobj = maxStars;
    }

    void setStars(lsst::afw::table::SourceCatalog const & srcs, int x0, int y0);

private:

    struct _Deleter {
        void operator()(solver_t* m) {
            solver_free(m);
        }
    };

    std::unique_ptr<solver_t, _Deleter> _solver;
};

/**
 * Calculate the distance from coordinates to a healpix
 *
 * Note that this assumes that the astrometry.net catalog reference system is ICRS.
 */
lsst::afw::geom::Angle healpixDistance(int hp, int nside, lsst::afw::coord::Coord const& coord);

}}}  // namespace lsst::meas::astrom