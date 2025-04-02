==============
RavenPy |logo|
==============

+----------------------------+-----------------------------------------------------+
| Versions                   | |pypi| |versions|                                   |
+----------------------------+-----------------------------------------------------+
| Documentation and Support  | |docs|                                              |
+----------------------------+-----------------------------------------------------+
| Open Source                | |license| |ossf-score|                              |
+----------------------------+-----------------------------------------------------+
| Coding Standards           | |black| |isort| |ruff| |ossf-bp| |pre-commit|       |
+----------------------------+-----------------------------------------------------+
| Development Status         | |status| |build| |coveralls|                        |
+----------------------------+-----------------------------------------------------+


A Python wrapper for configuring and running the hydrologic modelling framework Raven_.

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
.. _Cookiecutter: https://github.com/cookiecutter/cookiecutter
.. _Ouranos: https://www.ouranos.ca
.. _Ouranosinc/cookiecutter-pypackage: https://github.com/Ouranosinc/cookiecutter-pypackage
.. _Raven: https://raven.uwaterloo.ca
.. _docs: https://raven.uwaterloo.ca/files/v3.8/RavenManual_v3.8.pdf
.. _installation docs: https://ravenpy.readthedocs.io/en/latest/installation.html
.. _raven-hydro: https://github.com/Ouranosinc/raven-hydro


.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/psf/black
        :alt: Python Black

.. |build| image:: https://github.com/CSHS-CWRA/RavenPy/actions/workflows/main.yml/badge.svg
        :target: https://github.com/CSHS-CWRA/RavenPy/actions
        :alt: Build Status

.. |coveralls| image:: https://coveralls.io/repos/github/CSHS-CWRA/RavenPy/badge.svg
        :target: https://coveralls.io/github/CSHS-CWRA/RavenPy
        :alt: Coveralls

.. |docs| image:: https://readthedocs.org/projects/ravenpy/badge/?version=latest
        :target: https://ravenpy.readthedocs.io/en/latest
        :alt: Documentation Status

.. |isort| image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
        :target: https://pycqa.github.io/isort/
        :alt: Isort

.. |license| image:: https://img.shields.io/github/license/CSHS-CWRA/RavenPy.svg
        :target: https://github.com/CSHS-CWRA/RavenPy/blob/master/LICENSE
        :alt: License

.. |logo| image:: https://raw.githubusercontent.com/CSHS-CWRA/RavenPy/master/docs/_static/_images/logos/ravenpy-logo-small.png
        :target: https://github.com/CSHS-CWRA/RavenPy
        :alt: RavenPy Logo

.. |ossf-bp| image:: https://bestpractices.coreinfrastructure.org/projects/10064/badge
        :target: https://bestpractices.coreinfrastructure.org/projects/10064
        :alt: Open Source Security Foundation Best Practices

.. |ossf-score| image:: https://api.securityscorecards.dev/projects/github.com/CSHS-CWRA/RavenPy/badge
        :target: https://securityscorecards.dev/viewer/?uri=github.com/CSHS-CWRA/RavenPy
        :alt: OpenSSF Scorecard

.. |pre-commit| image:: https://results.pre-commit.ci/badge/github/CSHS-CWRA/RavenPy/master.svg
        :target: https://results.pre-commit.ci/latest/github/CSHS-CWRA/RavenPy/master
        :alt: pre-commit.ci status

.. |pypi| image:: https://img.shields.io/pypi/v/RavenPy.svg
        :target: https://pypi.python.org/pypi/RavenPy
        :alt: PyPI

.. |ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
        :target: https://github.com/astral-sh/ruff
        :alt: Ruff

.. |status| image:: https://www.repostatus.org/badges/latest/active.svg
        :target: https://www.repostatus.org/#active
        :alt: Project Status: Active - The project has reached a stable, usable state and is being actively developed.

.. |versions| image:: https://img.shields.io/pypi/pyversions/RavenPy.svg
        :target: https://pypi.python.org/pypi/RavenPy
        :alt: Supported Python Versions
