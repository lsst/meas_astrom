.. lsst-task-topic:: lsst.meas.astrom.DirectMatchTask

###############
DirectMatchTask
###############

.. Summary paragraph (a few sentences)
.. The aim is to say what the task is for

``DirectMatchTask`` implements a brute force nearest neighbor match between
a input source catalog and a reference catalog.

The matching permits no rotation or scaling, but uses the existing sky
positions in the source catalog. This is often useful for QA, as it allows
validating the pipeline astrometry and photometry against the reference
catalog.

Note that this DirectMatchTask is not currently suitable for use within the
AstrometryTask, as it has a different interface and serves a different purpose.

.. _lsst.meas.astrom.DirectMatchTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``DirectMatchTask`` runs this sequence of operations:

- Finds an on sky circle covering the input source catalog.
- Loads reference objects within the circle.
- Performs a simple nearest neighbor assignment within a specified tolerance.

.. _lsst.meas.astrom.DirectMatchTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.astrom.DirectMatchTask

.. _lsst.meas.astrom.DirectMatchTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.astrom.DirectMatchTask

.. _lsst.meas.astrom.DirectMatchTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.astrom.DirectMatchTask

.. _lsst.meas.astrom.DirectMatchTask-examples:

Examples
========

.. Add a brief example here.
.. If there are multiple examples
.. (such as one from a command-line context and another that uses the Python API)
.. you can separate each example into a different subsection for clarity.

.. code-block:: py

    config = DirectMatchConfig()
    task = DirectMatchTask(butler=butler, config=config)
    matchResults = task.run(catalog)
