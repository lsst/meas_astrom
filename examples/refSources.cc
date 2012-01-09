
extern "C" {
#include <stdlib.h>
#include <math.h>
}

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>
#include "boost/shared_ptr.hpp"

#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
//#include "lsst/meas/algorithms/Measure.h"

using namespace std;
namespace except = lsst::pex::exceptions;
namespace pexLog = lsst::pex::logging;
namespace pexPolicy = lsst::pex::policy;
namespace afwCoord = lsst::afw::coord;
namespace afwGeom = lsst::afw::geom;
namespace afwImg = lsst::afw::image;
namespace afwDet = lsst::afw::detection;
namespace det = lsst::afw::detection;
namespace math = lsst::afw::math;

namespace measAstrom = lsst::meas::astrom;

int main() {
    const char* candata = getenv("ASTROMETRY_NET_DATA_DIR");
    if (!candata) {
        printf("astrometry_net_data must be setup.\n");
        exit(-1);
    }
    string andata = string(candata);
    cout << "ANDATA: " << andata << endl;
    andata += "/metadata.paf";
    //pexPolicy::Policy pol(andata);
    measAstrom::net::GlobalAstrometrySolution gas(andata);

    measAstrom::net::ReferenceSources ref;

    ref = gas.getCatalogue(1.62002199864 * afwGeom::degrees, -3.84209636408 * afwGeom::degrees, 627.541610333 * afwGeom::arcseconds,
                           "g", "id", 101217104);
    cout << "Get refs: " << ref.refsources.size() << endl;

    vector<bool> stargal = gas.getTagAlongBool(ref.intrefsources, "starnotgal");
    cout << "Got stargals: " << stargal.size() << endl;

    for (unsigned int i = 0; i < stargal.size(); i++) {
        cout << (stargal[i] ? "T" : "F") << " ";
    }
    cout << endl;

    vector<bool> vars = gas.getTagAlongBool(ref.intrefsources, "variable");
    cout << "Got variability: " << vars.size() << endl;

    for (unsigned int i = 0; i < vars.size(); i++) {
        cout << (vars[i] ? "T" : "F") << " ";
    }
    cout << endl;

    vector<double> magerrs = gas.getTagAlongDouble(ref.intrefsources, "g_err");
    cout << "Got mag errs: " << magerrs.size() << endl;

    for (unsigned int i = 0; i < magerrs.size(); i++) {
        cout << magerrs[i] << " ";
    }
    cout << endl;

    for (unsigned int i = 0; i < ref.refsources.size(); i++) {
        ref.refsources[i]->setPsfFluxErr(magerrs[i] * ref.refsources[i]->getPsfFlux() * -log(10)/2.5);
        bool isstar = stargal[i];
        if (vars[i])
            isstar = false;
        //
        int starflag = 0x2000;

        ref.refsources[i]->setFlagForDetection(ref.refsources[i]->getFlagForDetection() | starflag);
    }

        
}
