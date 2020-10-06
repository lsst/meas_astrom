
import numpy as np
import subprocess

shifts = np.linspace(1, 300, 2)
rots = np.linspace(0, 6, 2)
scales = np.logspace(-5, -3, 2)
shears = np.logspace(-5, -3, 2)

for shift in shifts:
    job = subprocess.run(
                    "processCcd.py /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/output/ "
                    f"--output shift{shift:.2f}rot0scale0shear0 "
                    "--calib /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/calibingested/ "
                    "-C config/procCcdCalibrate.py --id visit=411420^419802^411371 "
                    f"-c calibrate.astrometry.shiftSize={shift:.2f} "
                    "-c calibrate.astrometry.matcher.maxRotationDeg=6.0 "
                    "-c calibrate.astrometry.outputFile=/project/morriscb/src/ap_verify_ci_hits2015/wcsTweak",
                    shell=True)
    job.wait()
for rot in rots:
    job = subprocess.run(
                    "processCcd.py /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/output/ "
                    f"--output shift0rot{rot:.2f}scale0shear0 "
                    "--calib /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/calibingested/ "
                    "-C config/procCcdCalibrate.py --id visit=411420^419802^411371 "
                    f"-c calibrate.astrometry.rotsize={rot:.2f} "
                    "-c calibrate.astrometry.matcher.maxRotationDeg=6.0 "
                    "-c calibrate.astrometry.outputFile=/project/morriscb/src/ap_verify_ci_hits2015/wcsTweak",
                    shell=True)
    job.wait()
for scale in scales:
    job = subprocess.run(
                    "processCcd.py /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/output/ "
                    f"--output shift0rot0scale{scale:.5f}shear0 "
                    "--calib /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/calibingested/ "
                    "-C config/procCcdCalibrate.py --id visit=411420^419802^411371 "
                    f"-c calibrate.astrometry.affineXScale={scale:.5f} "
                    "-c calibrate.astrometry.matcher.maxRotationDeg=6.0 "
                    "-c calibrate.astrometry.outputFile=/project/morriscb/src/ap_verify_ci_hits2015/wcsTweak",
                    shell=True)
    job.wait()
for shear in shears:
    job = subprocess.run(
                    "processCcd.py /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/output/ "
                    f"--output shift0rot0scale0shear{shear:.5f} "
                    "--calib /project/morriscb/src/ap_verify_ci_hits2015/defaultRun/calibingested/ "
                    "-C config/procCcdCalibrate.py --id visit=411420^419802^411371 "
                    f"-c calibrate.astrometry.affineXShear={shear:.5f} "
                    "-c calibrate.astrometry.matcher.maxRotationDeg=6.0 "
                    "-c calibrate.astrometry.outputFile=/project/morriscb/src/ap_verify_ci_hits2015/wcsTweak",
                    shell=True)
    job.wait()
