.. lsst-task-topic:: lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask

####################
MatchOptimisticBTask
####################

.. Summary paragraph (a few sentences)
.. The aim is to say what the task is for

``MatchOptimisticBTask`` matches sources to reference objects. This is often done
as a preliminary step to fitting an astrometric or photometric solution.

Optimistic Pattern Matching is described in [Tabur2007]_

.. _lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``MatchOptimisticBTask`` runs this sequence of operations:

- Flags sources with bad centroids and low signal to noise and remove them from
  the matching.
- Match the usable sources with an input reference catalog using the V. Tabur
  2007 algorithm.
- Further remove sources detected on the edge of the image and those that are
  saturated.
- Return these sources matched to the references.

.. _lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask

.. _lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask

.. _lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask

.. _lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask-examples:

Examples
========

.. Add a brief example here.
.. If there are multiple examples
.. (such as one from a command-line context and another that uses the Python API)
.. you can separate each example into a different subsection for clarity.

MatchOptimisticBTask is a subtask of AstrometryTask, which is called by
PhotoCalTask.

See :lsst-task:`lsst.pipe.tasks.photoCal.PhotoCalTask`
.. note:: Pipe task will require conversion before this link is useable.

.. _lsst.meas.astrom.matchOptimisticB.matchOptimisticBContinued.MatchOptimisticBTask-debug:

Debugging
=========

The `lsst.pipe.base.cmdLineTask.CmdLineTask` command line task interface supports a flag -d to import debug.py from your PYTHONPATH; see `lsstDebug` for more about debug.py files.

The available variables in MatchOptimisticB are

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

into your debug.py file and run this task with the --debug flag.
