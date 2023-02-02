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

.. warning::

   `OSTRICH <http://www.civil.uwaterloo.ca/envmodelling/Ostrich.html>`_ is another dependency needed for some key functions
   in RavenPy but is not currently available for Windows conda users. We suggest performing either a manual installation of OSTRICH
   or compiling it from source code: https://github.com/usbr/ostrich

For Unix/Linux users, OSTRICH can be installed from conda by running:

.. code-block:: console

   (ravenpy) $ conda install -c conda-forge ostrich

Custom Installation (Python/Pip)
--------------------------------

.. warning::
   The following instructions will only work on POSIX-like systems (Unix/Linux; not supported on Windows).

If you wish to install RavenPy and its C-libraries manually, compiling the `Raven` and `Ostrich binaries for your system,
you can install the entire system directly, placing them in the `bin` folder of your environment.
In order to perform this from Ubuntu/Debian:

.. code-block:: console

   $ sudo apt-get install gcc libnetcdf-dev gdal proj geos geopandas

Then, from your python environment, run:

.. code-block:: console

   $ pip install ravenpy[gis]
   $ pip install ravenpy[gis] --verbose --install-option="--with-binaries"

.. warning::

  It is imperative that the Python dependencies are pre-installed before running the `--with-binaries`
  option; This install step will fail otherwise.

If desired, the core functions of `RavenPy` can be installed without its GIS functionalities as well.
This implementation of RavenPy is much lighter on dependencies and can be installed easily with `pip`,
without the need for `conda` or `virtualenv`.

The only libraries required for RavenPy in this approach are a C++ compiler and the NetCDF4 development libraries.

.. code-block:: console

   $ sudo apt-get install gcc libnetcdf-dev

.. code-block:: console

   $ pip install ravenpy
   $ pip install ravenpy --verbose --install-option="--with-binaries"

.. warning::

  It is imperative that the Python dependencies are pre-installed before running the `--with-binaries`
  option; This install step will fail otherwise.

If for any reason you prefer to install without the binaries, from a fresh python environment, run the following:

.. code-block:: console

   (ravenpy-env) $ pip install ravenpy[gis]

But then you will be in charge of providing either ``raven`` and ``ostrich`` binaries on your PATH, or setting values for
``RAVENPY_RAVEN_BINARY_PATH`` and ``RAVENPY_OSTRICH_BINARY_PATH`` environment variables (both as absolute paths) in the
terminal/command prompt/shell used at runtime.

.. note::

  The `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ implementation also works well, but the
  GIS system libraries it depends on (specifically `GDAL` and `GEOS`) can be more difficult to configure.

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
   (ravenpy) $ pip install --editable ".[dev]"

.. warning::

   The following command will only work on POSIX-like systems (Unix/Linux; not supported on Windows).

.. code-block:: console

   # for the Raven and OSTRICH binaries
   (ravenpy) $ pip install --editable "." --install-option="--with-binaries"

Install the pre-commit hook (to make sure that any code you contribute is properly formatted):

.. code-block:: console

   (ravenpy-env) $ pre-commit install

If everything was properly installed the test suite should run successfully:

.. code-block:: console

   (ravenpy-env) $ pytest tests
