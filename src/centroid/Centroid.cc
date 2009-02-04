#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw.h"

#include "lsst/meas/astrom/centroid/Centroid.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

/*
 * Include concrete implementations
 */
#include "lsst/meas/astrom/centroid/NaiveCentroid.h"

namespace lsst { namespace meas { namespace astrom { namespace centroid {

/************************************************************************************************************/

template<typename ImageT>
Centroid Centroider<ImageT>::apply(ImageT const& image, int x, int y, double background) const {
    if (x < 1 || x > image.getWidth() - 2 || y < 1 || y > image.getHeight() - 2) {
        throw LSST_EXCEPT(pexExceptions::RangeErrorException,
                          (boost::format("Object at (%d, %d) is too close to the edge of the frame") % x % y).str());
    }
    pexLogging::TTrace<8>("meas.astrom.centroid", "Centroiding object at (%d, %d)", x, y);

    return doApply(image, x, y, background);
}

/************************************************************************************************************/
/**
 * @brief A factory function to return a Centroider of the specified type
 *
 * The Centroider has a method (apply) that can be used to return a Centroid
 */
template<typename ImageT>
Centroider<ImageT>* make_Centroider(centroidType type) {
    switch (type) {
      case NAIVE:
        return NaiveCentroider<ImageT>::getInstance();
      default:
        throw LSST_EXCEPT(pexExceptions::NotFoundException, 
                          (boost::format("Centroider of type %d is not implemented") % type).str());
    }
    // NOTREACHED
}

//
// Explicit instantiations
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
                template Centroid Centroider<IMAGE_T>::apply(IMAGE_T const&, int, int, double) const; \
                template Centroider<IMAGE_T>* make_Centroider<IMAGE_T>(centroidType);

                
MAKE_CENTROIDERS(lsst::afw::image::Image<float>)

// \endcond
                
}}}}
