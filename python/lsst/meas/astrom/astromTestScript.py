
import numpy as np
import subprocess

shifts = np.linspace(1, 299.99, 10)
rots = np.linspace(-5.99, 5.99, 10)
scales = np.logspace(0, np.log10(1.01), 10)
shears = np.logspace(-5, -3, 10)
print(shifts)
print(rots)
print(scales)
print(shears)

for shift in shifts:
    job = subprocess.run(
                    "processCcd.py /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/output/ "
                    f"--output shift{shift:.2f}rot0scale0shear0 "
                    "--calib /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/calibingested/ "
                    "-C config/procCcdCalibrate.py --id visit=411420^419802^411371 "
                    f"-c calibrate.astrometry.shiftSize={shift:.2f} -j6 "
                    "-c calibrate.astrometry.matcher.maxRotationDeg=5.999 "
                    "-c calibrate.astrometry.outputFile=/project/morriscb/src/ap_verify_ci_hits2015/wcsTweak",
                    shell=True)
for rot in rots:
    job = subprocess.run(
                    "processCcd.py /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/output/ "
                    f"--output shift0rot{rot:.2f}scale0shear0 "
                    "--calib /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/calibingested/ "
                    "-C config/procCcdCalibrate.py --id visit=411420^419802^411371 "
                    f"-c calibrate.astrometry.rotsize={rot:.2f} -j6 "
                    "-c calibrate.astrometry.matcher.maxRotationDeg=5.999 "
                    "-c calibrate.astrometry.outputFile=/project/morriscb/src/ap_verify_ci_hits2015/wcsTweak",
                    shell=True)
for scale in scales:
    job = subprocess.run(
                    "processCcd.py /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/output/ "
                    f"--output shift0rot0scale{scale:.5f}shear0 "
                    "--calib /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/calibingested/ "
                    "-C config/procCcdCalibrate.py --id visit=411420^419802^411371 "
                    f"-c calibrate.astrometry.affineXScale={scale:.5f} -j6 "
                    "-c calibrate.astrometry.matcher.maxRotationDeg=5.999 "
                    "-c calibrate.astrometry.outputFile=/project/morriscb/src/ap_verify_ci_hits2015/wcsTweak",
                    shell=True)
for shear in shears:
    job = subprocess.run(
                    "processCcd.py /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/output/ "
                    f"--output shift0rot0scale0shear{shear:.5f} "
                    "--calib /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/calibingested/ "
                    "-C config/procCcdCalibrate.py --id visit=411420^419802^411371 "
                    f"-c calibrate.astrometry.affineXShear={shear:.5f} -j6 "
                    "-c calibrate.astrometry.matcher.maxRotationDeg=5.999 "
                    "-c calibrate.astrometry.outputFile=/project/morriscb/src/ap_verify_ci_hits2015/wcsTweak",
                    shell=True)
