.. py:currentmodule:: lsst.gauss2d.fit

.. _lsst.gauss2d.fit:

###############
lsst.gauss2d.fit
###############

gauss2dfit is a submodule for gauss2d that implements a 2D Gaussian mixture
model class, along with its constituent parts, including parameters and data.
The Model class can evaluate the likelihood and gradients thereof. gauss2dfit
does not yet provide optimizers, although support for GSL fitters is planned.
Users should turn to `MultiProFit <https://github.com/lsst-dm/multiprofit>`_
for access to Python (scipy) optimizers.

Gaussian mixture approximations to the Sersic profile are provided. Use of
the GSL library is strongly recommended as it enables nonlinear interpolation
of the profile weights, which is necessary to compute accurate derivatives of
the likelihood and model for the Sersic index parameter.

.. _lsst.gauss2d.fit-using:

Using lsst.gauss2d.fit
==============================

Example usage can be found in the unit tests and also in dependent packages,
particularly `MultiProFit <https://github.com/lsst-dm/multiprofit>`_.

.. toctree::
   :maxdepth: 2

.. _lsst.gauss2d.fit-contributing:

Contributing
============

``lsst.gauss2d.fit`` is developed at https://github.com/lsst-dm/gauss2d_fit.
You can find Jira issues for this module under the
`gauss2d_fit <https://rubinobs.atlassian.net/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20gauss2d_fit>`_
component.

.. If there are topics related to developing this module (rather than using it), link to this from a toctree placed here.

.. .. toctree::
..    :maxdepth: 2

.. _lsst.gauss2d.fit-pyapi:

Python API reference
====================

``lsst.gauss2d.fit`` has Python bindings for classes using numpy-based single
and double precision arrays. Support for GSL arrays is forthcoming with
`DM-38617 <https://jira.lsstcorp.org/browse/DM-38617>`_.

.. automodapi:: lsst.gauss2d.fit
   :no-main-docstr:
   :no-inheritance-diagram:
