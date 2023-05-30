============
Installation
============

Anaconda Python Installation
----------------------------

For many reasons, we recommend using a `Conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
to work with the full RavenPy installation. This implementation is able to manage the harder-to-install GIS dependencies, like `GDAL`.

Begin by creating an environment:

.. code-block:: console

   $ conda create -c conda-forge --name ravenpy

The newly created environment must then be activated:

.. code-block:: console

   $ conda activate ravenpy

RavenPy can then be installed directly via its `conda-forge` package by running:

.. code-block:: console

   (ravenpy) $ conda install -c conda-forge ravenpy

This approach installs the `Raven <http://raven.uwaterloo.ca>`_ binary directly to your environment `PATH`,
as well as installs all the necessary Python and C libraries supporting GIS functionalities.

Python Installation (pip)
-------------------------

.. warning::

   In order to compile the Raven model (provided by the `raven-hydro` package, a C++ compiler (`GCC`, `Clang`, `MSVC`, etc.) and either `GNU Make` (Linux/macOS) or `Ninja` (Windows) must be exposed on the `$PATH`.

.. warning::

   The Raven model also requires that NetCDF4 libraries are installed on the system, exposed on the `$PATH`, and discoverable using the `FindNetCDF.cmake` helper script bundled with `raven-hydro`.

   On Linux, this can be provided by the `libnetcdf-dev` system library; On macOS by the `netcdf` homebrew package; And on Windows by using UNIDATA's `pre-built binaries <https://docs.unidata.ucar.edu/netcdf-c/current/winbin.html>`_.

In order to perform this from Ubuntu/Debian:

.. code-block:: console

   $ sudo apt-get install gcc libnetcdf-dev gdal proj geos

Then, from your python environment, run:

.. code-block:: console

   $ pip install ravenpy[gis]

If desired, the core functions of `RavenPy` can be installed without its GIS functionalities as well. This implementation of RavenPy is much lighter on dependencies and can be installed easily with `pip`, without the need for `conda` or `virtualenv`.

.. code-block:: console

   $ pip install ravenpy

Using A Custom Raven Model Binary
---------------------------------

If you wish to install the `Raven` model, either compiling the `Raven` binary from sources for your system or installing the pre-built binary offered by UWaterloo, we encourage you to consult the `Raven` documentation (`Raven Downloads <https://www.civil.uwaterloo.ca/raven/Downloads.html>`_).

Once downloaded/compiled, the binary can be pointed to manually (as an absolute path) by setting the environment variable ``RAVENPY_RAVEN_BINARY_PATH`` in the terminal/command prompt/shell used at runtime.

.. code-block:: console

   $ export RAVENPY_RAVEN_BINARY_PATH=/path/to/my/custom/raven

Development Installation (from sources)
---------------------------------------

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
   (ravenpy) $ pip install --editable ".[dev,gis]"

Install the pre-commit hook (to make sure that any code you contribute is properly formatted):

.. code-block:: console

   (ravenpy-env) $ pre-commit install

If everything was properly installed the test suite should run successfully:

.. code-block:: console

   (ravenpy-env) $ pytest tests
