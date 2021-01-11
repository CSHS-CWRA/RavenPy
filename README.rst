=======
RavenPy
=======


.. image:: https://img.shields.io/pypi/v/ravenpy.svg
        :target: https://pypi.python.org/pypi/ravenpy

..
    .. image:: https://img.shields.io/travis/CSHS-CWRA/ravenpy.svg
            :target: https://travis-ci.com/CSHS-CWRA/ravenpy

.. image:: https://readthedocs.org/projects/ravenpy/badge/?version=latest
        :target: https://ravenpy.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

..
    .. image:: https://pyup.io/repos/github/CSHS-CWRA/ravenpy/shield.svg
            :target: https://pyup.io/repos/github/CSHS-CWRA/ravenpy/
            :alt: Updates


A Python wrapper to setup and run the hydrologic modelling framework Raven_.


* Free software: MIT license
* Documentation: https://ravenpy.readthedocs.io.


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

RavenPy can only be installed in a virtual environment, see `installation docs`_ .

.. code-block:: bash

    sudo apt install libnetcdf-dev libspatialindex-dev libgdal-dev
    python -m venv ./venv
    source ./venv/bin/activate
    pip install . --verbose --install-option="--with-raven" --install-option="--with-testdata"


Grid Weight Generation Script
-----------------------------

To run the grid weight generation script (written and maintained by `Julie Mai <https://github.com/julemai/GridWeightsGenerator>`_):

.. code-block:: bash

    ravenpy generate-grid-weights /path/to/HRU/shapefile /path/to/NC/file --var-names lon lat --output raven


Credits
-------

RavenPy's development has been funded by CANARIE_.

This package was created with Cookiecutter_ and the `Ouranosinc/cookiecutter-pypackage`_ project template.


.. _`installation docs`: https://ravenpy.readthedocs.io/installation
.. _Raven: http://raven.uwaterloo.ca
.. _Cookiecutter: https://github.com/audreyfeldroy/cookiecutter-pypackage
.. _`Ouranosinc/cookiecutter-pypackage`: https://github.com/Ouranosinc/cookiecutter-pypackage
.. _`CANARIE`: https://www.canarie.ca
