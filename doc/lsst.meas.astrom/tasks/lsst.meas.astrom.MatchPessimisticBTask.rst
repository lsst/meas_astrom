.. lsst-task-topic:: lsst.meas.astrom.MatchPessimisticBTask

#####################
MatchPessimisticBTask
#####################

.. Summary paragraph (a few sentences)
.. The aim is to say what the task is for

``MatchPessimisticBTask`` matches sources to reference objects. This is often
done as a preliminary step to fitting an astrometric or photometric solution.

The algorithm is based on a more "Pessimistic" version of the Optimistic
Pattern Matcher B as described in `DMTN-013 <http://ls.st/DMTN-031>`_.

Optimistic Pattern Matching is described in [Tabur2007]_

.. [Tabur2007] Fast algorithms for matching CCD images to a stellar catalogue*
               `arxiv:0710.3618 <https://arxiv.org/abs/0710.3618>`_

.. _lsst.meas.astrom.MatchPessimisticBTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``MatchPessimisticBTask`` runs this sequence of operations:

- Flags sources with bad centroids and low signal to noise and remove them from
  the matching.
- Match the usable sources with an input reference catalog using the updated
  V. Tabur 2007 algorithm.
- Further remove sources detected on the edge of the image and those that are
  saturated.
- Return these sources matched to the references.

.. _lsst.meas.astrom.MatchPessimisticBTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.astrom.MatchPessimisticBTask

.. _lsst.meas.astrom.MatchPessimisticBTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.meas.astrom.MatchPessimisticBTask

.. _lsst.meas.astrom.MatchPessimisticBTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.astrom.MatchPessimisticBTask

.. _lsst.meas.astrom.MatchPessimisticBTask-examples:

Examples
========

.. Add a brief example here.
.. If there are multiple examples
.. (such as one from a command-line context and another that uses the Python API)
.. you can separate each example into a different subsection for clarity.

MatchPessimisticBTask is a subtask of AstrometryTask, which is called by
PhotoCalTask.

See :lsst-task:`lsst.pipe.tasks.photoCal.PhotoCalTask`
.. note:: Pipe task will require conversion before this link is usable.

.. _lsst.meas.astrom.MatchPessimisticBTask-debug:

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
