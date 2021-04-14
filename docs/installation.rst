.. highlight:: shell

============
Installation
============

Full Installation (Anaconda)
----------------------------

For many reasons, we recommend using a `Conda environment
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
to work with the full RavenPy installation. This implementation is able to manage
the harder to install GIS dependencies, like `GDAL`.

.. code-block:: console

   $ conda create -c conda-forge --name ravenpy-env click clisops matplotlib pip proj statsmodels xarray xclim xskillscore

The newly created environment must then be activated:

.. code-block:: console

   $ conda activate ravenpy-env

To install the remaining Python dependencies, run:

.. code-block::

   (ravenpy-env) $ pip install ravenpy[gis]

`RavenPy` relies on the `Raven <http://raven.uwaterloo.ca>`_ and `OSTRICH
<http://www.civil.uwaterloo.ca/envmodelling/Ostrich.html>`_ binaries, which can be conveniently
downloaded, compiled, and made available on your `PATH`. There are currently two options available:

* To directly download, compile, and place them in the `bin` folder of your environment, use this command:

.. code-block:: console

   (ravenpy-env) $ pip install ravenpy --verbose --install-option="--with-binaries"

.. warning::

  It is imperative that the Python dependencies are pre-installed before running the `--with-binaries`
  option; This install step will fail otherwise.

* Alternatively, the Raven and Ostrich binaries can be installed directly from `conda` with the following command:

.. code-block::

  (ravenpy-env) $ conda install -c zeitsperre raven ostrich

If successful, this should install ``raven`` and ``ostrich`` binaries in the ``bin``
folder of your environment, which should then be already available in your
path, and that you can verify by doing:

.. code-block:: console

   (ravenpy-env) $ which raven
   (ravenpy-env) $ which ostrich

If for any reason you prefer to install without the binaries, you can
simply omit the option:

.. code-block:: console

   (ravenpy-env) $ pip install ravenpy[gis]

But then you will be in charge of providing either ``raven`` and ``ostrich`` binaries on your PATH,
or values for ``RAVENPY_RAVEN_BINARY_PATH`` and ``RAVENPY_OSTRICH_BINARY_PATH`` environment
variables (both as absolute paths) at runtime.

.. note::

  The `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ implementation also works well, but the
  GIS system libraries it depends on (specifically `GDAL` and `GEOS`) can be more difficult to configure.

Light Installation
------------------

If desired, the core functions of `RavenPy` can be installed without its GIS functionalities as well.
This implementation of RavenPy is much lighter on dependencies and can be installed easily with `pip`,
without the need for `conda` or `virtualenv`.

.. code-block:: console

  $ pip install ravenpy
  $ pip install ravenpy --verbose --install-option="--with-binaries"

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
   $ conda activate ravenpy-env

You can then install RavenPy with:

.. code-block:: console

   # for the python dependencies
   (ravenpy-env) $ pip install --editable ".[dev]"
   # for the Raven and OSTRICH binaries
   (ravenpy-env) $ pip install --editable "." --install-option="--with-binaries"

Then clone the Raven Test Data repo somewhere on your disk:

.. code-block:: console

   (ravenpy-env) $ git clone git@github.com:Ouranosinc/raven-testdata.git

Install the pre-commit hook (to make sure that any code you contribute is properly formatted):

.. code-block:: console

   (ravenpy-env) $ pre-commit install

If everything was properly installed the test suite should run successfully:

.. code-block:: console

   (ravenpy-env) $ export RAVENPY_TESTDATA_PATH=/path/to/raven-testdata
   (ravenpy-env) $ pytest tests

Or set the conda environment variable permanently:

.. code-block:: console

   (ravenpy-env) $ conda env config vars set RAVENPY_TESTDATA_PATH=/path/to/raven-testdata

then deactivate and reactivate the environment.
