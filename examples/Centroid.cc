//
// Demonstrate use of Centroider
//
#include <iostream>
#include "lsst/afw.h"
#include "lsst/meas/astrom/centroid/Centroid.h"

using namespace std;
namespace centroid = lsst::meas::astrom::centroid;
namespace image = lsst::afw::image;

template<typename ImageT>
void computeCentroid(centroid::Centroider<ImageT> const* cc) {
    centroid::Centroid cen = cc->apply(ImageT(100, 100), 10, 20);

    cout << "(x, y) = " << cen.getX() << ", " << cen.getY() << endl;
}

int main() {
    typedef image::Image<float> ImageT;
    centroid::Centroider<ImageT> *nc = centroid::make_Centroider<ImageT>(centroid::SDSS);

    computeCentroid(nc);
}
