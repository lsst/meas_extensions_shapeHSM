.. py:currentmodule:: lsst.meas.extensions.shapeHSM

.. _lsst.meas.extensions.shapeHSM:

#############################
lsst.meas.extensions.shapeHSM
#############################

The ``lsst.meas.extensions.shapeHSM`` module provides algorithms for HSM shape measurement.
The algorithm was initially described in `Hirata & Seljak (2003) <https://ui.adsabs.harvard.edu/abs/2003MNRAS.343..459H/abstract>`_, and was modified later in `Mandelbaum et al. (2005) <https://ui.adsabs.harvard.edu/abs/2005MNRAS.361.1287M/abstract>`_.
HSM is named after the primary authors: Christopher Hirata, Uros Seljak, and Rachel Mandelbaum.
Their implementation of this algorithm lives within `GalSim <https://github.com/GalSim-developers/GalSim>`_, and this package interacts with the Python layer of GalSim to make the measurements.

Using lsst.meas.extensions.shapeHSM
===================================

.. toctree linking to topics related to using the module's APIs.

.. toctree::
   :maxdepth: 1

   hsm-moments
   hsm-shape

.. _lsst.meas.extensions.shapeHSM-contributing:

Contributing
============

``lsst.meas.extensions.shapeHSM`` is developed at https://github.com/lsst/meas_extensions_shapeHSM.
You can find Jira issues for this module under the `meas_extensions_shapeHSM <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20meas_extensions_shapeHSM>`_ component.

.. _lsst.meas.extensions.shapeHSM-pyapi:

Python API reference
====================

.. automodapi:: lsst.meas.extensions.shapeHSM
   :no-main-docstr:
