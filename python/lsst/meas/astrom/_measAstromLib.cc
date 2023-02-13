#include "pybind11/pybind11.h"
#include "lsst/cpputils/python.h"

namespace py = pybind11;
using namespace pybind11::literals;
using lsst::cpputils::python::WrapperCollection;

namespace lsst {
namespace meas {
namespace astrom {
namespace sip {
void wrapCreateWcsWithSip(WrapperCollection &wrappers);
void wrapLeastSqFitter1d(WrapperCollection &wrappers);
void wrapLeastSqFitter2d(WrapperCollection &wrappers);
void wrapMatchSrcToCatalogue(WrapperCollection &wrappers);
}  // namespace sip

void wrapPolynomialTransform(WrapperCollection &wrappers);
void wrapScaledPolynomialTransformFitter(WrapperCollection &wrappers);
void wrapSipTransform(WrapperCollection &wrappers);
void wrapMatchOptimisticB(WrapperCollection &wrappers);
void wrapMakeMatchStatistics(WrapperCollection &wrappers);
void wrapPessimisticPatternMatcherUtils(WrapperCollection &wrappers);

PYBIND11_MODULE(_measAstromLib, mod) {
    WrapperCollection wrappers(mod, "lsst.meas.astrom");

    wrappers.addInheritanceDependency("lsst.afw.math");

    sip::wrapLeastSqFitter1d(wrappers);
    sip::wrapCreateWcsWithSip(wrappers);
    sip::wrapMatchSrcToCatalogue(wrappers);
    sip::wrapLeastSqFitter2d(wrappers);

    wrapPolynomialTransform(wrappers);
    wrapScaledPolynomialTransformFitter(wrappers);
    wrapSipTransform(wrappers);
    wrapMatchOptimisticB(wrappers);
    wrapMakeMatchStatistics(wrappers);
    wrapPessimisticPatternMatcherUtils(wrappers);
    wrappers.finish();
};

}  // namespace astrom
}  // namespace meas
}  // namespace lsst