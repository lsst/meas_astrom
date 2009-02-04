#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/meas/astrom/centroid/SdssCentroid.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace astrom { namespace centroid {

template<typename ImageT> SdssCentroider<ImageT>* SdssCentroider<ImageT>::_instance = 0;

template<typename ImageT>
Centroid SdssCentroider<ImageT>::doApply(ImageT const& image, int x, int y, double background) const {
    typename ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1, -1) + im(-1, 0) + im(-1,  1) +
         im( 0, -1) + im( 0, 0) + im( 0,  1) +
         im( 1, -1) + im( 1, 0) + im( 1,  1)) - 9*background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::UnderflowErrorException,
                          (boost::format("Object at (%d, %d) has no counts") % x % y).str());
    }

    double const sum_x =
        -im(-1, -1) + im(-1,  1) +
        -im( 0, -1) + im( 0,  1) +
        -im( 1, -1) + im( 1,  1);
    double const sum_y =
        -(im(-1, -1) + im(-1, 0) + im(-1,  1)) +
          im( 1, -1) + im( 1, 0) + im( 1,  1);

    return Centroid(x + sum_x/sum, y + sum_y/sum);
}

//
// Explicit instantiations
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
                template class SdssCentroider<lsst::afw::image::Image<IMAGE_T> >;
                
MAKE_CENTROIDERS(float)

// \endcond

}}}}
