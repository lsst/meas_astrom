# -*- python -*-
#
# Setup our environment
#
import glob, os.path, sys
import lsst.SConsUtils as scons

env = scons.makeEnv(
    "meas_astrom",
    r"$HeadURL$",
    [
        ["boost", "boost/version.hpp", "boost_system:C++"],
        ["boost", "boost/version.hpp", "boost_filesystem:C++"],
        ["boost", "boost/regex.hpp", "boost_regex:C++"],
        ["boost", "boost/test/unit_test.hpp", "boost_unit_test_framework:C++"],
        ["python", "Python.h"],
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
        ["afw", "lsst/afw/image/MaskedImage.h", "afw"],
        ["astrometry_net", "solver.h", "pthread backend"], 
    ],
)

env.libs["meas_astrom"] +=  env.getlibs("daf_base daf_data daf_persistence pex_logging pex_exceptions " + \
    "pex_policy security minuit2 afw boost utils wcslib astrometry_net")
if True:
    #
    # Workaround SConsUtils failure to find numpy .h files. Fixed in sconsUtils >= 3.3.2
    #
    import numpy
    env.Append(CCFLAGS = ["-I", numpy.get_include()])
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
            print >> sys.stderr, traceback.format_exc()
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
