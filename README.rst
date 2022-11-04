=======
RavenPy
=======

.. image:: https://img.shields.io/pypi/v/ravenpy.svg
    :target: https://pypi.python.org/pypi/ravenpy
    :alt: PyPI

.. image:: https://anaconda.org/conda-forge/ravenpy/badges/version.svg
    :target: https://anaconda.org/conda-forge/ravenpy
    :alt: Conda-Forge

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

`RavenPy` is a Python wrapper for Raven_, accompanied by utility functions that facilitate model configuration, calibration, and evaluation.

Raven_ is an hydrological modeling framework that lets hydrologists build hydrological models by combining different hydrological processes together. It can also be used to emulate a variety of existing lumped and distributed models. Model structure, parameters, initial conditions and forcing files are configured in text files, which Raven parses to build and run hydrological simulations. A detailed description about modeling capability of Raven_ can be found in the `docs`_.

`RavenPy` provides a Python interface to Raven_, automating the creation of configuration files and allowing the model to be launched from Python. Results, or errors, are automatically parsed and exposed within the programming environment. This facilitates the launch of parallel simulations, multi-model prediction ensembles, sensitivity analyses and other experiments involving a large number of model runs.

The code is currently undergoing a number of changes to support semi-distributed watersheds, so expect some API changes over the next versions.

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

Acknowledgements
----------------

RavenPy's development has been funded by CANARIE_ and Ouranos_ and would be not be possible without the help of Juliane Mai and James Craig.

This package was created with Cookiecutter_ and the `Ouranosinc/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyfeldroy/cookiecutter-pypackage
.. _Raven: http://raven.uwaterloo.ca
.. _`CANARIE`: https://www.canarie.ca
.. _`Ouranos`: http://ouranos.ca
.. _`Ouranosinc/cookiecutter-pypackage`: https://github.com/Ouranosinc/cookiecutter-pypackage
.. _`docs`: http://raven.uwaterloo.ca/files/v3.0.1/RavenManual_v3.0.1.pdf
.. _`installation docs`: https://ravenpy.readthedocs.io/en/latest/installation.html
