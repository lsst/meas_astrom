# -*- python -*-
#
# Setup our environment
#
import glob, os.path, sys
import lsst.SConsUtils as scons


# Stolen from ap/SConstruct.
boostInt64IsLongCheckSrc = """
    #include "boost/cstdint.hpp"
    #include "boost/static_assert.hpp"
    #include "boost/type_traits/is_same.hpp"

    int main() {
        BOOST_STATIC_ASSERT((boost::is_same<long, boost::int64_t>::value));
        return 0;
    }
    """
def CustomCompileCheck(context, message, source, extension = '.c'):
    context.Message(message)
    result = context.TryCompile(source, extension)
    context.Result(result)
    return result



env = scons.makeEnv(
    "meas_astrom",
    r"$HeadURL$",
    [
        ["boost", "boost/version.hpp", "boost_system:C++"],
        ["boost", "boost/version.hpp", "boost_filesystem:C++"],
        ["boost", "boost/regex.hpp", "boost_regex:C++"],
        ["boost", "boost/test/unit_test.hpp", "boost_unit_test_framework:C++"],
        ["python", "Python.h"],
        ["base", "lsst/base.h"],
        #["numpy", "Python.h numpy/arrayobject.h"], # see numpy workaround below
        ["m", "math.h", "m", "sqrt"],
        ["cfitsio", "fitsio.h", "cfitsio", "ffopen"],
        ["wcslib", "wcslib/wcs.h", "wcs"],
        ["xpa", "xpa.h", "xpa", "XPAPuts"],
        ["minuit2", "Minuit2/FCNBase.h", "Minuit2:C++"],
        ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
        ["utils", "lsst/utils/Utils.h", "utils:C++"],
        ["daf_base", "lsst/daf/base.h", "daf_base:C++"],
        ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
        ["security", "lsst/security/Security.h", "security:C++"],
        ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
        ["daf_persistence", "lsst/daf/persistence/Persistence.h", "daf_persistence:C++"],
        ["daf_data", "lsst/daf/data/LsstBase.h", "daf_data:C++"],
        ["eigen", "Eigen/Core.h"],
        ["ndarray", "lsst/ndarray/Array.h"],
        ["afw", "lsst/afw/image/MaskedImage.h", "afw"],
        ["astrometry_net", "solver.h", "pthread backend"], 
    ],
)

env.libs["meas_astrom"] +=  env.getlibs("daf_base daf_data daf_persistence pex_logging pex_exceptions " + \
    "pex_policy security minuit2 ndarray afw boost utils wcslib astrometry_net")
if True:
    #
    # Workaround SConsUtils failure to find numpy .h files. Fixed in sconsUtils >= 3.3.2
    #
    import numpy
    env.Append(CCFLAGS = ["-I", numpy.get_include()])

if not env.CleanFlagIsSet():
    conf = Configure(env, custom_tests = {'CustomCompileCheck' : CustomCompileCheck, })
                                          
    # Without some help, SWIG disagrees with boost on the actual type of int64_t
    if conf.CustomCompileCheck('Checking whether boost::int64_t is long ... ',
                               boostInt64IsLongCheckSrc, extension='.cc'):
        conf.env["MYSWIGFLAGS"] = '-DSWIGWORDSIZE64'
        env = conf.Finish()

#
# Build/install things
#
for d in (
    ".",
    "doc",
    "examples",
    "lib",
    "python/lsst/meas/astrom/net",
    "python/lsst/meas/astrom/sip",
    "tests",
):
    if d != ".":
        try:
            SConscript(os.path.join(d, "SConscript"))
        except Exception, e:
            print >> sys.stderr, "In processing file %s:" % (os.path.join(d, "SConscript"))
            print >> sys.stderr, e
    Clean(d, Glob(os.path.join(d, "*~")))
    Clean(d, Glob(os.path.join(d, "*.pyc")))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", [
    env.InstallAs(os.path.join(env['prefix'], "doc", "doxygen"), os.path.join("doc", "htmlDir")),
    env.Install(env['prefix'], "examples"),
    env.Install(env['prefix'], "include"),
    env.Install(env['prefix'], "lib"),
    env.Install(env['prefix'], "python"),
    env.InstallEups(env['prefix'] + "/ups"),
])

scons.CleanTree(r"*~ core *.so *.os *.o")
#
# Build TAGS files
#
files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
Wcs determination
""")
