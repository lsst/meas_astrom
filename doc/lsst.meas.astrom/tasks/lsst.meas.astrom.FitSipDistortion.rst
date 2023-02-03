.. lsst-task-topic:: lsst.meas.astrom.FitSipDistortionTask

####################
FitSipDistortionTask
####################

.. Summary paragraph (a few sentences)
.. The aim is to say what the task is for

``FitSipDistortionTask`` is a drop-in replacement for
:lsst-task:`lsst.meas.astrom.FitTanSipWcsTask`.  It is built on fundamentally
stronger fitting algorithms, but has received significantly less testing.
    
Like :lsst-task:`lsst.meas.astrom.FitTanSipWcsTask`, this task is most easily
used as the wcsFitter component of
:lsst-task:`lsst.meas.astrom.AstrometryTask`; it can be enabled in a config
file via e.g.

.. code-block:: py

    from lsst.meas.astrom import FitSipDistortionTask
    config.(...).astometry.wcsFitter.retarget(FitSipDistortionTask)

.. _lsst.meas.astrom.FitSipDistortionTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``FitSipDistortionTask`` involves three steps:

- We set the CRVAL and CRPIX reference points to the mean positions of
  the matches, while holding the CD matrix fixed to the value passed in
  to the run() method.  This work is done by the makeInitialWcs method.i
- We fit the SIP "reverse transform" (the AP and BP polynomials that map
  "intermediate world coordinates" to pixels).  This happens iteratively;
  while fitting for the polynomial coefficients given a set of matches is
  a linear operation that can be done without iteration, outlier
  rejection using sigma-clipping and estimation of the intrinsic scatter
  are not. By fitting the reverse transform first, we can do outlier
  rejection in pixel coordinates, where we can better handle the source
  measurement uncertainties that contribute to the overall scatter.  This
  fit results in a
  :cpp:class:`lsst::meas::astrom::ScaledPolynomialTransform`, which is
  somewhat more general than the SIP reverse transform in that it allows
  an affine transform both before and after the polynomial.  This is
  somewhat more numerically stable than the SIP form, which applies only
  a linear transform (with no offset) before the polynomial and only a
  shift afterwards.  We only convert to SIP form once the fitting is
  complete.  This conversion is exact (though it may be subject to
  significant round-off error) as long as we do not attempt to null the
  low-order SIP polynomial terms (we do not).
- Once the SIP reverse transform has been fit, we use it to populate a
  grid of points that we use as the data points for fitting its inverse,
  the SIP forward transform.  Because our "data" here is artificial,
  there is no need for outlier rejection or uncertainty handling.  We
  again fit a general scaled polynomial, and only convert to SIP form
  when the fit is complete.

.. _lsst.meas.astrom.FitSipDistortionTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.astrom.FitSipDistortionTask

.. _lsst.meas.astrom.FitSipDistortionTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.astrom.FitSipDistortionTask

.. _lsst.meas.astrom.FitSipDistortionTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.astrom.FitSipDistortionTask

.. _lsst.meas.astrom.FitSipDistortionTask-debug:

Debugging
=========

Enabling DEBUG-level logging on this task will report the number of
outliers rejected and the current estimate of intrinsic scatter at each
iteration.

FitSipDistortionTask also supports the following :ref:`lsstDebug` variables to
control diagnostic displays:

- FitSipDistortionTask.display: if True, enable display diagnostics.
- FitSipDistortionTask.frame: frame to which the display will be sent
- FitSipDistortionTask.pause: whether to pause (by dropping into pdb)
  between iterations (default is True).  If False, multiple frames
  will be used, starting at the given number.

The diagnostic display displays the image (or an empty image if
exposure=None) overlaid with the positions of sources and reference
objects will be shown for every iteration in the reverse transform fit.
The legend for the overlay is:

Red X
    Reference sources transformed without SIP distortion terms; this
    uses a TAN WCS whose CRPIX, CRVAL and CD matrix are the same
    as those in the TAN-SIP WCS being fit.  These are not expected to
    line up with sources unless distortion is small.

Magenta X
    Same as Red X, but for matches that were rejected as outliers.

Red O
    Reference sources using the current best-fit TAN-SIP WCS.  These
    are connected to the corresponding non-distorted WCS position by
    a red line, and should be a much better fit to source positions
    than the Red Xs.

Magenta O
    Same as Red O, but for matches that were rejected as outliers.

Green Ellipse
    Source positions and their error ellipses, including the current
    estimate of the intrinsic scatter.

Cyan Ellipse
    Same as Green Ellipse, but for matches that were rejected as outliers.

Reference to parameters:
See :lsst-task:`lsst.pipe.base.Task`; FitSipDistortionTask does not add any
additional constructor parameters.
