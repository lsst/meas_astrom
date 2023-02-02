.. lsst-task-topic:: lsst.meas.astrom.AstrometryTask

##############
AstrometryTask
##############

.. Summary paragraph (a few sentences)
.. The aim is to say what the task is for

``AstrometryTask`` matches a source catalog with objects from a reference
catalog and attempts to fit a WCS.

.. _lsst.meas.astrom.AstrometryTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``AstrometryTask`` runs this sequence of operations:

- Find position reference stars that overlap the exposure.
- Match input source catalog to the position reference stars.
- Fit a WCS based on the matches.

.. _lsst.meas.astrom.AstrometryTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.astrom.AstrometryTask

.. _lsst.meas.astrom.AstrometryTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.astrom.AstrometryTask

.. _lsst.meas.astrom.AstrometryTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.astrom.AstrometryTask

.. _lsst.meas.astrom.AstrometryTask-debug:

Debugging
=========

The `lsst.pipe.base.cmdLineTask.CmdLineTask` command line task interface supports a flag -d to import debug.py from your PYTHONPATH; see :ref:`lsstDebug` for more about debug.py files.

The available variables in AstrometryTask are

display (bool)
    If True display information at three stages: after finding reference objects, after matching sources to reference objects, and after fitting the WCS; defaults to False
frame (int)
    frame to use to display the reference objects; the next two frames are used to display the match list and the results of the final WCS; defaults to 0

To investigate the meas_astrom_astrometry_Debug, put something like

.. code-block:: py

    import lsstDebug
    def DebugInfo(name):
        debug = lsstDebug.getInfo(name) # N.b. lsstDebug.Info(name) would call us recursively
        if name == "lsst.meas.astrom.astrometry":
            debug.display = True

        return debug

    lsstDebug.Info = DebugInfo

into your debug.py file and run this task with the ``--debug`` flag.
