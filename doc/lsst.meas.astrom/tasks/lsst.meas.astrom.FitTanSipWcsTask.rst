.. lsst-task-topic:: lsst.meas.astrom.FitTanSipWcsTask

################
FitTanSipWcsTask
################

.. Summary paragraph (a few sentences)
.. The aim is to say what the task is for

``FitTanSipWcsTask`` Fit a TAN-SIP WCS given a list of reference object/source
matches.

.. _lsst.meas.astrom.FitTanSipWcsTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

Measure the distortions in an image plane and express them a SIP polynomials.

Given a list of matching sources between a catalog and an image,
and a linear Wcs that describes the mapping from pixel space in the image
and ra/dec space in the catalog, calculate discrepancies between the two
and compute SIP distortion polynomials to describe the discrepancy.

SIP polynomials are defined in Shupe at al. (2005) ASPC 347 491.

Note that the SIP standard insists (although it is only mentioned obliquely 
between Eqns 3 and 4) that the lowest three terms in the distortion
polynomials be zero (A00, A10, A01, B00, etc.). To achieve this, we need to
adjust the values of CD and CRPIX from the input wcs. This may not be the
behavior you expect.

.. _lsst.meas.astrom.FitTanSipWcsTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.astrom.FitTanSipWcsTask

.. _lsst.meas.astrom.FitTanSipWcsTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.astrom.FitTanSipWcsTask

.. _lsst.meas.astrom.FitTanSipWcsTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.astrom.FitTanSipWcsTask

.. _lsst.meas.astrom.FitTanSipWcsTask-examples:

Examples
========

.. Add a brief example here.
.. If there are multiple examples
.. (such as one from a command-line context and another that uses the Python API)
.. you can separate each example into a different subsection for clarity.

See :lsst-task:`lsst.pipe.tasks.photoCal.PhotoCalTask`
.. note:: Pipe task will require conversion before this link is usable.

.. _lsst.meas.astrom.FitTanSipWcsTask-debug:

Debugging
=========

FitTanSipWcsTask does not support any debug variables.
