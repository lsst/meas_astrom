#if !defined(LSST_MEAS_ASTROM_NAIVECENTROID_H)
#define LSST_MEAS_ASTROM_NAIVECENTROID_H 1
/**
 * @file
 */
#include "lsst/meas/astrom/centroid/Centroid.h"

namespace lsst { namespace meas { namespace astrom { namespace centroid {

/**
 * @brief Return the centroid as a simple unweighted first moment of the 3x3 region around a pixel
 */
template<typename ImageT>
class NaiveCentroider : public Centroider<ImageT> {
public:
    static Centroider<ImageT>* getInstance() {
        if (_instance == NULL) {
            _instance = new NaiveCentroider;
        }
        return _instance;
    }
private:
    Centroid doApply(ImageT const& image, int x, int y, double background) const;

    static NaiveCentroider* _instance;
};
}}}}
#endif
