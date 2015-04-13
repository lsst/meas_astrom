// -*- lsst-c++ -*-
#if !defined(LSST_MEAS_ASTROM_FITTANWCS_H)
#define LSST_MEAS_ASTROM_FITTANWCS_H

#include "lsst/afw/image.h"
#include "lsst/afw/table.h"

namespace lsst {
namespace meas {
namespace astrom {

    /**
    Fit a pure TAN-WCS

    @param[in] matches  reference object/source matches, such as is returned by matchOptimisticB
        (an lsst::afw::table::ReferenceMatchVector)
    @param[in] verbose  if true print diagnostic messages to std::cout
    */
    PTR(lsst::afw::image::Wcs) fitTanWcs(
        lsst::afw::table::ReferenceMatchVector const &matches,
        bool verbose = false
    );

}}}

#endif
