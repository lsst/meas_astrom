#if !defined(LSST_MEAS_ASTROM_UTILS_H)
#define LSST_MEAS_ASTROM_UTILS_H 1

#include <string>
#include <vector>

#include "lsst/afw/table/Source.h"

namespace lsst {
namespace meas {
namespace astrom {
namespace detail {

    typedef struct {
        std::string name;
        std::string magcol;
        std::string magerrcol;

        bool hasErr() const {
            return (magerrcol.size() > 0);
        }
    } mag_column_t;


/// RAII manager for astrometry.net indices
///
/// Ensures index files are closed when done, to prevent "Too many open files" errors.
/// See also #2292, #2879.
struct IndexManager {
    index_t* index;
    IndexManager(index_t* ind) : index(ind) {}
    ~IndexManager() {
        // Change once astrometry.net-0.40+ is in...
        /*
          if (index_close_fds(ind)) {
          throw LSST_EXCEPT(lsst::pex::exceptions::IoError,
          "Failed to index_close_fds() an astrometry_net_data file");
          }
        */
        if (index->quads && index->quads->fb && index->quads->fb->fid) {
            _close(index->quads->fb->fid);
        }
        if (index->codekd && index->codekd->tree && index->codekd->tree->io) {
            _close(index->codekd->tree->io);
        }
        if (index->starkd && index->starkd->tree && index->starkd->tree->io) {
            _close(index->starkd->tree->io);
        }
    }
    void _close(FILE* & fid) {
        if (fid) {
            if (fclose(fid)) {
                std::cerr << "Error closing an astrometry_net_data quadfile" << std::endl;
            }
            fid = NULL;
        }
    }
    void _close(void* io) {
        _close(reinterpret_cast<kdtree_fits_t*>(io)->fid);
    }
};



/*
 * Implementation for index_t::getCatalog method
 */
lsst::afw::table::SimpleCatalog
getCatalogImpl(std::vector<index_t*> inds,
               double ra, double dec, double radius,
               const char* idcol,
               std::vector<mag_column_t> const& magcols,
               const char* stargalcol,
               const char* varcol,
               bool unique_ids);
}}}}
#endif
