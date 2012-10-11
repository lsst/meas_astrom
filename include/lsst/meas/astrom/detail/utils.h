#if !defined(LSST_MEAS_ASTROM_UTILS_H)
#define LSST_MEAS_ASTROM_UTILS_H 1

#include <string>
#include <vector>

#include "lsst/afw/table/Source.h"

extern "C" {
    typedef struct index_s index_t;     // from astrometry.net's index.h
}

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

/*
 * Implementation for index_s::getCatalog method
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
