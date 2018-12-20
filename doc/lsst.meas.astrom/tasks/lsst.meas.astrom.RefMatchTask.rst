.. lsst-task-topic:: lsst.meas.astrom.RefMatchTask

##############
RefMatchTask
##############

.. Summary paragraph (a few sentences)
.. The aim is to say what the task is for

``RefMatchTask`` matches a source catalog with objects from a reference
catalog.

.. _lsst.meas.astrom.RefMatchTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``RefMatchTask`` runs this sequence of operations:

- Find position reference stars that overlap the exposure.
- Match input source catalog to the position reference stars.
- Return statistics on quality of match.

.. _lsst.meas.astrom.RefMatchTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.astrom.RefMatchTask

.. _lsst.meas.astrom.RefMatchTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.astrom.RefMatchTask

.. _lsst.meas.astrom.RefMatchTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.astrom.RefMatchTask

.. _lsst.meas.astrom.RefMatchTask-debug:
