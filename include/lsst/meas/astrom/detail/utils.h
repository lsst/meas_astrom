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
/*
 * Implementation for index_s::getCatalog method
 */
lsst::afw::table::SimpleCatalog
getCatalogImpl(std::vector<index_t*> inds,
	       double ra, double dec, double radius,
	       const char* idcol,
	       std::vector<std::string> const & magcolVec,
	       std::vector<std::string> const & magerrcolVec,
	       const char* stargalcol,
	       const char* varcol,
	       bool unique_ids);
}}}}
#endif
