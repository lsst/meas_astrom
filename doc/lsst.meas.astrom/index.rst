.. py:currentmodule:: lsst.meas.astrom

.. _lsst.meas.astrom:

################
lsst.meas.astrom
################

Tasks and methods for finding an astrometric solution, including:

- `~lsst.meas.astrom.AstrometryTask` compute an astrometric solution for a set of sources detected
    on an Exposure

Default subtasks used by `~lsst.meas.astrom.AstrometryTask`:

- `~lsst.meas.astrom.MatchOptimisticBTask` match sources to reference objects
- `~lsst.meas.astrom.FitTanSipWcsTask` fit a TAN-SIP WCS given a list of matches of sources
    and reference objects


.. .. _lsst.meas.astrom-using:

.. Using lsst.meas.astrom
.. ==================

.. toctree linking to topics related to using the module's APIs.

.. .. toctree::
..    :maxdepth: 1

.. _lsst.meas.astrom-contributing:

Contributing
============

``lsst.meas.astrom`` is developed at https://github.com/lsst/meas_astrom.
You can find Jira issues for this module under `here <https://jira.lsstcorp.org/issues/?jql=text%20~%20%22meas_astrom%22>`_.

.. If there are topics related to developing this module (rather than using it), link to this from a toctree placed here.

.. .. toctree::
..    :maxdepth: 1

.. _lsst.meas.astrom-command-line-taskref:

Task reference
==============

.. _lsst.meas.astrom-command-line-tasks:

Command-line tasks
------------------

.. lsst-cmdlinetasks::
   :root: lsst.meas.astrom

.. _lsst.meas.astrom-tasks:

Tasks
-----

.. lsst-tasks::
   :root: lsst.meas.astrom
   :toctree: tasks

.. _lsst.meas.astrom-configs:

Configurations
--------------

.. lsst-configs::
   :root: lsst.meas.astrom
   :toctree: configs

.. _lsst.meas.astrom-pyapi:

Python API reference
====================

.. automodapi:: lsst.meas.astrom
   :no-main-docstr:
   :no-inheritance-diagram:
