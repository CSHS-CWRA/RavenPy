=======
RavenPy
=======

|pypi| |conda| |license| |build| |docs| |coveralls|

A Python wrapper to setup and run the hydrologic modelling framework Raven_.

* Free software: MIT license
* Documentation: https://ravenpy.readthedocs.io

`RavenPy` is a Python wrapper for Raven_, accompanied by utility functions that facilitate model configuration, calibration, and evaluation.

Raven_ is an hydrological modeling framework that lets hydrologists build hydrological models by combining different hydrological processes together. It can also be used to emulate a variety of existing lumped and distributed models. Model structure, parameters, initial conditions and forcing files are configured in text files, which Raven parses to build and run hydrological simulations. A detailed description about modeling capability of Raven_ can be found in the `docs`_.

`RavenPy` provides a Python interface to Raven_, automating the creation of configuration files and allowing the model to be launched from Python. Results, or errors, are automatically parsed and exposed within the programming environment. This facilitates the launch of parallel simulations, multi-model prediction ensembles, sensitivity analyses and other experiments involving a large number of model runs.

Note that version 0.12 includes major changes compared to the previous 0.11 release, and breaks backward compatibility. The benefits of these changes are a much more intuitive interface for configuring and running the model.

Features
--------

* Configure, run and parse Raven outputs from Python
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

.. _CANARIE: https://www.canarie.ca
.. _Cookiecutter: https://github.com/audreyfeldroy/cookiecutter-pypackage
.. _Ouranos: https://www.ouranos.ca
.. _Ouranosinc/cookiecutter-pypackage: https://github.com/Ouranosinc/cookiecutter-pypackage
.. _Raven: http://raven.uwaterloo.ca
.. _docs: https://www.civil.uwaterloo.ca/raven/files/v3.7/RavenManual_v3.7.pdf
.. _installation docs: https://ravenpy.readthedocs.io/en/latest/installation.html
.. _raven-hydro: https://github.com/Ouranosinc/raven-hydro


.. |pypi| image:: https://img.shields.io/pypi/v/ravenpy.svg
        :target: https://pypi.python.org/pypi/ravenpy
        :alt: PyPI

.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/ravenpy
        :target: https://anaconda.org/conda-forge/ravenpy
        :alt: Conda-Forge

.. |license| image:: https://img.shields.io/github/license/CSHS-CWRA/RavenPy.svg
        :target: https://github.com/CSHS-CWRA/RavenPy/blob/master/LICENSE
        :alt: License

.. |build| image:: https://github.com/CSHS-CWRA/RavenPy/actions/workflows/main.yml/badge.svg
        :target: https://github.com/CSHS-CWRA/RavenPy/actions/workflows/main.yml
        :alt: Build status

.. |docs| image:: https://readthedocs.org/projects/ravenpy/badge/?version=latest
        :target: https://ravenpy.readthedocs.io
        :alt: Documentation Status

.. |coveralls| image:: https://coveralls.io/repos/github/CSHS-CWRA/RavenPy/badge.svg?branch=master
        :target: https://coveralls.io/github/CSHS-CWRA/RavenPy?branch=master
        :alt: Coveralls
