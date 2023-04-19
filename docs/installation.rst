============
Installation
============

Full Installation (Anaconda)
----------------------------

For many reasons, we recommend using a `Conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
to work with the full RavenPy installation. This implementation is able to manage the harder-to-install GIS dependencies,like `GDAL`.
Begin by creating an environment:

.. code-block:: console

   $ conda create -c conda-forge --name ravenpy

The newly created environment must then be activated:

.. code-block:: console

   $ conda activate ravenpy

RavenPy can then be installed directly via its `conda-forge` package by running:

.. code-block:: console

   (ravenpy) $ conda install -c conda-forge ravenpy

This approach installs both the `Raven <http://raven.uwaterloo.ca>`_ binary directly to your environment `PATH`,
as well as installs all the necessary Python and C libraries supporting GIS functionalities.


Custom Installation (Python/Pip)
--------------------------------

.. warning::
   As of April 2023, this installation method does not install the `raven-hydro` hydrological model. You must install the `raven-hydro` model and add it to your `bin` or append it to your `$PATH` manually. This may change in the future

If you wish to install RavenPy and its C-libraries manually, compiling the `Raven` binary for your system,
you can install the entire system directly, placing them in the `bin` folder of your environment.
In order to perform this from Ubuntu/Debian:

.. code-block:: console

   $ sudo apt-get install gcc libnetcdf-dev gdal proj geos geopandas

Then, from your python environment, run:

.. code-block:: console

   $ pip install ravenpy[gis]

If desired, the core functions of `RavenPy` can be installed without its GIS functionalities as well.
This implementation of RavenPy is much lighter on dependencies and can be installed easily with `pip`,
without the need for `conda` or `virtualenv`.

.. code-block:: console

   $ pip install ravenpy

But then you will be in charge of providing ``raven``  binaries on your PATH, or setting values for the
``RAVENPY_RAVEN_BINARY_PATH`` environment variable (as an absolute path) in the
terminal/command prompt/shell used at runtime.

Development Installation (from sources)
---------------------------------------

.. warning::
   As of April 2023, this installation method does not install the `raven-hydro` hydrological model. You must install the `raven-hydro` model and add it to your `bin` or append it to your `$PATH` manually. This may change in the future

The sources for RavenPy can be obtained from the GitHub repo:

.. code-block:: console

    $ git clone git://github.com/CSHS-CWRA/ravenpy

You can then create and activate your `Conda environment
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
by doing:

.. code-block:: console

   $ cd /path/to/ravenpy
   $ conda env create -f environment.yml
   $ conda activate ravenpy

You can then install RavenPy with:

.. code-block:: console

   # for the python dependencies
   (ravenpy) $ pip install --editable ".[dev]"

Install the pre-commit hook (to make sure that any code you contribute is properly formatted):

.. code-block:: console

   (ravenpy-env) $ pre-commit install

If everything was properly installed the test suite should run successfully:

.. code-block:: console

   (ravenpy-env) $ pytest tests
