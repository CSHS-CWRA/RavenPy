============
Installation
============

Anaconda Python Installation
----------------------------

For many reasons, we recommend using a `Conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ to work with the full RavenPy installation. This implementation is able to manage the harder-to-install GIS dependencies, like `GDAL`.
Furthermore, due to the complexity of some packages, the default dependency solver can take a long time to resolve the environment.
If `mamba` is not already your default solver, consider running the following commands in order to speed up the process:

    .. code-block:: console

        conda install -n base conda-libmamba-solver
        conda config --set solver libmamba

Begin by creating an environment:

    .. code-block:: console

        conda create -c conda-forge --name ravenpy

The newly created environment must then be activated:

    .. code-block:: console

        conda activate ravenpy

RavenPy can then be installed directly via its `conda-forge` package by running:

    .. code-block:: console

        (ravenpy) $ conda install -c conda-forge ravenpy

This approach installs the `Raven <https://raven.uwaterloo.ca>`_ binary directly to your environment `PATH`, as well as installs all the necessary Python and C libraries supporting GIS functionalities.

Python Installation (pip)
-------------------------

.. warning::

   In order to compile the Raven model (provided by the `raven-hydro` package, a C++ compiler (`GCC`, `Clang`, `MSVC`, etc.) and either `GNU Make` (Linux/macOS) or `Ninja` (Windows) must be exposed on the `$PATH`.

.. warning::

   The Raven model also requires that NetCDF4 libraries are installed on the system, exposed on the `$PATH`, and discoverable using the `FindNetCDF.cmake` helper script bundled with `raven-hydro`.

   On Linux, this can be provided by the `libnetcdf-dev` system library; On macOS by the `netcdf` homebrew package; And on Windows by using UNIDATA's `pre-built binaries <https://docs.unidata.ucar.edu/netcdf-c/current/winbin.html>`_.

In order to perform this from Ubuntu/Debian:

    .. code-block:: console

        sudo apt-get install gcc libnetcdf-dev gdal proj geos

Then, from your python environment, run:

    .. code-block:: console

        python -m pip install ravenpy[gis,raven-hydro]

If desired, the core functions of `RavenPy` can be installed without its GIS functionalities as well. This implementation of RavenPy is much lighter on dependencies and can be installed easily with `pip`, without the need for `conda` or `virtualenv`.

    .. code-block:: console

        python -m pip install ravenpy[raven-hydro]

Finally, if you wish to provide your own `Raven` binary, you can install `RavenPy` without installing the `raven-hydro` package:

    .. code-block:: console

        python -m pip install ravenpy

Using A Custom Raven Model Binary
---------------------------------

If you wish to install the `Raven` model, either compiling the `Raven` binary from sources for your system or installing the pre-built binary offered by UWaterloo, we encourage you to consult the `Raven` documentation (`Raven Downloads <https://www.civil.uwaterloo.ca/raven/Downloads.html>`_).

Once downloaded/compiled, the binary can be pointed to manually (as an absolute path) by setting the environment variable ``RAVENPY_RAVEN_BINARY_PATH`` in the terminal/command prompt/shell used at runtime.

    .. code-block:: console

        export RAVENPY_RAVEN_BINARY_PATH=/path/to/my/custom/raven

Customizing remote service datasets
-----------------------------------

A number of functions and tests within `RavenPy` are dependent on remote services (THREDDS, GeoServer) for providing climate datasets, hydrological boundaries, and other data. These services are provided by `Ouranos <https://www.ouranos.ca>`_ through the `PAVICS <https://pavics.ouranos.ca>`_ project and may be subject to change in the future.

If for some reason you wish to use alternate services, you can set the following environment variables to point to your own instances of THREDDS and GeoServer:

    .. code-block:: console

        export RAVENPY_THREDDS_URL=https://my.domain.org/thredds
        export RAVENPY_GEOSERVER_URL=https://my.domain.org/geoserver

Development Installation (from sources)
---------------------------------------

The sources for `RavenPy` can be obtained from the GitHub repo:

#. Download the source code from the `Github repo`_ using one of the following methods:

    * Clone the public repository:

        .. code-block:: console

            git clone git://github.com/CSHS-CWRA/ravenpy

    * Download the `tarball <https://github.com/CSHS-CWRA/RavenPy/tarball/master>`_:

        .. code-block:: console

            curl -OJL https://github.com/CSHS-CWRA/RavenPy/tarball/master

#. Once you have a copy of the source, you can install it with:

    .. code-block:: console

        conda env create -f environment-dev.yml
        conda activate ravenpy-dev
        (ravenpy-dev) make dev

    If you are on Windows, replace the ``make dev`` command with the following:

    .. code-block:: console

        (ravenpy-dev) python -m pip install -e .[dev]

    Even if you do not intend to contribute to `RavenPy`, we favor using `environment-dev.yml` over `environment.yml` because it includes additional packages that are used to run all the examples provided in the documentation.
    If for some reason you wish to install the `PyPI` version of `RavenPy` into an existing Anaconda environment (*not recommended if requirements are not met*), only run the last command above.

#. When new changes are made to the `Github repo`_, if using a clone, you can update your local copy using the following commands from the root of the repository:

    .. code-block:: console

        git fetch
        git checkout main
        git pull origin main
        conda env update -n ravenpy-dev -f environment-dev.yml
        conda activate ravenpy-dev
        (ravenpy-dev) make dev

    These commands should work most of the time, but if big changes are made to the repository, you might need to remove the environment and create it again.

#. If everything was properly installed the test suite should run successfully:

    .. code-block:: console

        (ravenpy) python -m pytest tests

.. _Github repo: https://github.com/CSHS-CWRA/RavenPy
