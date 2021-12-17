=======
RavenPy
=======

.. image:: https://img.shields.io/pypi/v/ravenpy.svg
    :target: https://pypi.python.org/pypi/ravenpy

.. image:: https://anaconda.org/conda-forge/ravenpy/badges/installer/conda.svg
    :target: https://conda.anaconda.org/conda-forge

.. image:: https://github.com/CSHS-CWRA/RavenPy/actions/workflows/main.yml/badge.svg
    :target: https://github.com/CSHS-CWRA/RavenPy/actions/workflows/main.yml
    :alt: Build status

.. image:: https://readthedocs.org/projects/ravenpy/badge/?version=latest
    :target: https://ravenpy.readthedocs.io/en/latest/?version=latest
    :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/CSHS-CWRA/RavenPy/badge.svg?branch=master
    :target: https://coveralls.io/github/CSHS-CWRA/RavenPy?branch=master
    :alt: Coveralls

A Python wrapper to setup and run the hydrologic modelling framework Raven_.

* Free software: MIT license
* Documentation: https://ravenpy.readthedocs.io

Features
--------

* Download and compile Raven with `pip`
* Configure, run and parse Raven outputs from Python
* Parallel simulations over parameters, models or watersheds
* Utility command to create grid weight files
* Extract physiographic information about watersheds
* Algorithms to estimate model parameters from ungauged watersheds
* Exposes outputs (flow, storage) as `xarray.DataArray` objects

Install
-------

Please see the detailed `installation docs`_.

Credits
-------

RavenPy's development has been funded by CANARIE_.

This package was created with Cookiecutter_ and the `Ouranosinc/cookiecutter-pypackage`_ project template.

.. _`installation docs`: https://ravenpy.readthedocs.io/en/latest/installation.html
.. _Raven: http://raven.uwaterloo.ca
.. _Cookiecutter: https://github.com/audreyfeldroy/cookiecutter-pypackage
.. _`Ouranosinc/cookiecutter-pypackage`: https://github.com/Ouranosinc/cookiecutter-pypackage
.. _`CANARIE`: https://www.canarie.ca
