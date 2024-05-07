=======
History
=======

0.13.1 (2024-05-07)
-------------------
* Fixed buggy `CustomOutput` command. (PR #360)
* Make sure config and output paths are absolute. (PR #360)

0.13.0 (2023-01-10)
-------------------
* Fixed problem with scalar elevation in netCDF files parsed with `nc_specs`. (issue #279, PR #323)
* Added notebook on sensitivity analysis. (PR #320)
* Updated Notebooks 03 and 04. (PR #319)
* Upgrade to `pydantic` v2.0. (PR #326)
* Pin `cf-xarray` for Python3.8. (PR #325)
* Fix `Coveralls` Workflows. (PR #328)
* Fix notebook execution. (PR #329)
* Refactor and simplify testing data fetching. (PR #332)

Breaking changes
^^^^^^^^^^^^^^^^
* Update to `pydantic` v2.0. (PR #326)
* Added `h5netcdf` as a core dependency to provide a stabler backend for `xarray.open_dataset`. (PR #332)
* Switched from `autodoc_pydantic` to `autodoc-pydantic` for `pydantic` v2.0+ support in documentation. (PR #326)

Internal changes
^^^^^^^^^^^^^^^^
* Removed some redundant `pytest` fixtures for running `emulators` tests.
* `"session"`-scoped `pytest` fixtures used for hindcasting/forecasting are now always yielded and copied to new objects within tests.

0.12.3 (2023-10-02)
-------------------
* `RavenPy` now uses `platformdirs` to write `raven_testing` to the user's cache directory. Dynamic paths are now used to cache data dependent on the user's operating system. Developers can now safely delete the `.raven_testing_data` folder in their home directory without affecting the functionality of `RavenPy`.
* Updated `raven-hydro` to v0.2.4 to address CMake build issues.

Breaking changes
^^^^^^^^^^^^^^^^
* In tests, set `xclim`'s missing value option to ``skip``. As of `xclim` v0.45, missing value checks are applied to the ``fit`` indicator, meaning that parameters will be set to `None` if missing values are found in the fitted time series. Wrap calls to ``fit`` with ``xclim.set_options(check_missing="skip")`` to reproduce the previous behavior of xclim.
* The `_determine_upstream_ids` function under `ravenpy.utilities.geoserver` has been removed as it was a duplicate of `ravenpy.utilities.geo.determine_upstream_ids`. The latter function is now used in its place.

Internal changes
^^^^^^^^^^^^^^^^
* Added a GitHub Actions workflow to remove obsolete GitHub Workflow cache files.
* `RavenPy` now accepts a `RAVENPY_THREDDS_URL` for setting the URL globally to the THREDDS-hosted climate data service. Defaults to `https://pavics.ouranos.ca/twitcher/ows/proxy/thredds`.
* `RavenPy` processes and tests that depend on remote GeoServer calls now allow for optional server URL and file location targets. The server URL can be set globally with the following environment variable:
    * `RAVENPY_GEOSERVER_URL`: URL to the GeoServer-hosted vector/raster data. Defaults to `https://pavics.ouranos.ca/geoserver`. This environment variable was previously called `GEO_URL` but was renamed to narrow its scope to `RavenPy`.
        * `GEO_URL` is still supported for backward compatibility but may eventually be removed in a future release.
* `RavenPy` has temporarily pinned `xarray` below v2023.9.0 due to incompatibilities with `xclim` v0.45.0`.

0.12.2 (2023-07-04)
-------------------
This release is primarily a bugfix to address issues arising from dependencies.

Breaking changes
^^^^^^^^^^^^^^^^
* `raven-hydro` version has been bumped from v0.2.1 to v0.2.3. This version provides better support for builds on Windows and MacOS.
* Due to major breaking changes, `pydantic` has been pinned below v2.0 until changes can be made to adapt to their new API.
* `numpy` has been pinned below v1.25.0 to ensure compatibility with `numba`.

Internal changes
^^^^^^^^^^^^^^^^
* ``test_geoserver::test_select_hybas_ar_domain_point`` is now temporarily skipped when testing on MacOS due to a mysterious domain identification error.

0.12.1 (2023-06-01)
-------------------
This release is largely a bugfix to better stabilize performance and enhance the documentation.

* Avoid repeatedly calling `xr.open_dataset` in `OutputReader`'s `hydrograph` and `storage` properties. This seems to cause kernel failures in Jupyter notebooks.

Internal changes
^^^^^^^^^^^^^^^^
* Hyperlinks to documented functions now points to entries in the `User API` section.
* Docstrings are now more conformant to numpy-docstring conventions and formatting errors raised from badly-formatted pydantic-style docstrings have been addressed.
* In order to prevent timeout and excessive memory usage, Jupyter notebooks have been adjusted to no longer run on ReadTheDocs. All notebooks have been updated to the latest RavenPy and remain tested against RavenPy externally.
* Documentation built on ReadTheDocs is now set to `fail_on_warning`.

0.12.0 (2023-05-25)
-------------------
This release includes major breaking changes. It completely overhauls how models are defined, and how to run
simulations, and any code relying on the previous release will most likely break. Please check the documentation
to see how to use the new improved interface.

Breaking changes
^^^^^^^^^^^^^^^^
* The entire model configuration and simulation interface (see PR #269).
* The Raven model executable is now updated to v3.7.
* Added support for Ensemble Kalman Filter using RavenC.
* Now employing the `spotpy` package for model calibration instead of `ostrich`.
* BasinMaker importer assumes `SubBasin=HRU` in order to work with files downloaded from the BasinMaker web site.
* Ravenpy now employs a new method for installing the Raven model using the `raven-hydro <https://github.com/Ouranosinc/raven-hydro>`_ python package  (based on `scikit-build-core`) (see PR #278).
* Replaced `setup.py`, `requirements.txt`, and `Manifest.in` for `PEP 517 <https://peps.python.org/pep-0517>`_ compliance (`pyproject.toml`) using the flit backend (see PR #278).
* Dealt with an import-based error that occurred due to the sequence in which modules are loaded at import (attempting to call ravenpy before it is installed).
* Updated pre-commit hooks to include formatters and checkers for TOML files.
* The build recipes no longer build on each other, so when installing the `dev` or `docs` recipe, you must also install the gis recipe.
* Updated the GeoServer API calls to work with the GeoPandas v0.13.0.

0.11.0 (2023-02-16)
-------------------

* Update RavenC executable to v3.6.
* Update xclim library to v0.40.0.
* Update fiona library to v1.9.
* Address some failures that can be caused when attempting to run CLI commands without the proper GIS dependencies installed.
* Addressed warnings raised in conda-forge compilation due to badly-configured MANIFEST.in.
* Update installation documentation to reflect most recent changes.

0.10.0 (2022-12-21)
-------------------

* Update Raven executable to 3.5. Due to a bug in RavenC, simulations storing reservoir information to netCDF will fail. We expect this to be resolved in the next release. Note that we only test RavenPy with one Raven version. There is no guarantee it will work with other versions.
* Relax geo test to avoid failures occurring due to GDAL 3.6.
* Pin numpy below 1.24 (see https://github.com/numba/numba/issues/8615)

0.9.0 (2022-11-16)
------------------

Breaking changes
^^^^^^^^^^^^^^^^
* HRUState's signature has changed. Instead of passing variables as keyword arguments (e.g. `soil0=10.`), it now expects a `state` dictionary keyed by variables' Raven name (e.g. `{"SOIL[0]": 10}). This change makes `rvc` files easier to read, and avoids Raven warnings regarding 'initial conditions for state variables not in model'.
* `nc_index` renamed to `meteo_idx` to enable the specification of distinct indices for observed streamflow using `hydro_idx`. `nc_index` remains supported for backward compatibility.
* The distributed python testing library, `pytest-xdist` is now a testing and development requirement.
* `xarray` has been pinned below "2022.11.0" due to incompatibility with `climpred=="2.2.0"`.

New features
^^^^^^^^^^^^
* Add support for hydrometric gauge data distinct from meteorological input data. Configuration parameter `hydro_idx` identifies the gauge station index, while `meteo_idx` (previously `nc_index`) stands for the meteo station index.
* Add support for multiple gauge observations. If a list of `hydro_idx` is provided, it must be accompanied with a list of corresponding subbasin identifiers (`gauged_sb_ids`) of the same length.
* Automatically infer scale and offset `:LinearTransform` parameters from netCDF file metadata, so that input data units are automatically converted to Raven-compliant units whenever possible.
* Add support for the command `:RedirectToFile`. Tested for grid weights only.
* Add support for the command `:WriteForcingFunctions`.
* Add support for the command `:CustomOutput`.
* Multiple other new RavenCommand objects added, but not integrated in the configuration, including `:SoilParameterList`, `:VegetationParameterList` and `:LandUseParameterList`.
* Multichoice options (e.g. calendars) moved from RV classes to `config.options`, but aliases created for backward compatibility.
* Patch directory traversal vulnerability (`CVE-2007-4559 <https://github.com/advisories/GHSA-gw9q-c7gh-j9vm>`_).
* A local copy of the raven-testdata with environment variable (`RAVENPY_TESTDATA_PATH`) set to that location is now no longer needed in order to run the testing suite. Test data is fetched automatically and now stored at `~/.raven_testing_data`.
* RavenPy now leverages `pytest-xdist` to distribute tests among Python workers and significantly speed up the testing suite, depending on number of available CPUs. File access within the testing suite has also been completely rewritten for thread safety.
    - On pytest launch with "`--numprocesses` > 0", testing data will be fetched automatically from `Ouranosinc/raven-testdata` by one worker, blocking others until this step is complete. Spawned pytest workers will then copy the testing data to their respective temporary directories before beginning testing.
* To aid with development and debugging purposes, two new environment variables and pytest fixtures are now available:
    - In order to skip the data collection step: `export SKIP_TEST_DATA=true`
    - In order to target a specific branch of `Ouranosinc/raven-testdata` for data retrieval: `export MAIN_TESTDATA_BRANCH="my_branch"`
    - In order to fetch testing data using the user-set raven-testdata branch, pytest fixtures for `get_file` and `get_local_testdata` are now available for convenience

0.8.1 (2022-10-26)
------------------

* Undo change related to `suppress_output`, as it breaks multiple tests in raven. New `Raven._execute` method runs models but does not parse results.

0.8.0
-----

Breaking changes
^^^^^^^^^^^^^^^^
* Parallel parameters must be provided explicitly using the `parallel` argument when calling emulators.
* Multiple `nc_index` values generate multiple *gauges*, instead of being parallelized.
* Python3.7 is no longer supported.
* Documentation now uses sphinx-apidoc at build-time to generate API pages.

* Add ``generate-hrus-from-routing-product`` script.
* Do not write RV zip file and merge outputs when `suppress_output` is True. Zipping rv files during multiple calibration runs leads to a non-linear performance slow-down.
* Fixed issues with coverage reporting via tox and GitHub Actions
* Add partial support for `:RedirectToFile` command, tested with GridWeights only.

0.7.8
-----

* Added functionalities in Data Assimilation utils and simplified tests.
* Removed pin on setuptools.
* Fixed issues related to symlinks, working directory, and output filenames.
* Fixed issues related to GDAL version handling in conda-forge.
* Updated jupyter notebooks.

0.7.7
-----

* Updated internal shapely calls to remove deprecated ``.to_wkt()`` methods.

0.7.6
-----

* Automate release pipeline to PyPI using GitHub CI actions.
* Added coverage monitoring GitHub CI action.
* Various documentation adjustments.
* Various metadata adjustments.
* Pinned owslib to 0.24.1 and above.
* Circumvented a bug in GitHub CI that was causing tests to fail at collection stage.

0.7.5
-----

* Update test so that it works with xclim 0.29.

0.7.4
-----

* Pinned climpred below v2.1.6.

0.7.3
-----

* Pinned xclim below v0.29.

0.7.2
-----

* Update cruft.
* Subclass ``derived_parameters`` in Ostrich emulators to avoid having to pass ``params``.

0.7.0
-----

* Add support for V2.1 of the Routing Product in ``ravenpy.extractors.routing_product``.
* Add ``collect-subbasins-upstream-of-gauge`` CLI script.
* Modify WFS request functions to use spatial filtering (``Intersects``) supplied by OWSLib.

0.6.0
-----

* Add support for EvaluationPeriod commands. Note that as a result of this, the model's ``diagnostics`` property contains one list per key, instead of a single scalar. Also note that for calibration, Ostrich will use the first period and the first evaluation metric.
* Add ``SACSMA``, ``CANADIANSHIELD`` and ``HYPR`` model emulators.

0.5.2
-----

* Simplify RVC configuration logic.
* Add ``ravenpy.utilities.testdata.file_md5_checksum`` (previously in ``xarray.tutorial``).

0.5.1
-----

* Some adjustments and bugfixes needed for RavenWPS.
* Refactoring of some internal logic in ``ravenpy.config.rvs.RVT``.
* Improvements to typing with the help of mypy.

0.5.0
-----

* Refactoring of the RV config subsystem:

  * The config is fully encapsulated into its own class: ``ravenpy.config.rvs.Config``.
  * The emulator RV templates are inline in their emulator classes.

* The emulators have their own submodule: ``ravenpy.models.emulators``.
* The "importers" have been renamed to "extractors" and they have their own submodule: ``ravenpy.extractors``.

0.4.2
-----

* Update to RavenC revision 318 to fix OPeNDAP access for StationForcing commands.
* Fix grid_weights set to None by default.
* Pass nc_index to ObservationData command.
* Expose more cleanly RavenC errors and warnings.

0.4.1
-----

* Add notebook about hindcast verification skill.
* Add notebook about routing capability.
* Modify geoserver functions to have them return GeoJSON instead of GML.
* Collect upstream watershed aggregation logic.
* Fix RVC bug.

0.4.0
-----

This is an interim version making one step toward semi-distributed modeling support.
Model configuration is still in flux and will be significantly modified with 0.5.
The major change in this version is that model configuration supports passing multiple HRU objects,
instead of simply passing area, latitude, longitude and elevation for a single HRU.

* GR4JCN emulator now supports routing mode.
* Add BLENDED model emulator.
* DAP links for forcing files are now supported.
* Added support for ``tox``-based localized installation and testing with python-pip.
* Now supporting Python 3.7, 3.8, and 3.9.
* Build testing for ``pip`` and ``conda``-based builds with GitHub CI.

0.3.1
-----

* Update external dependencies (Raven, OSTRICH) to facilitate Conda packaging.

0.3.0
-----

* Migration and refactoring of GIS and IO utilities (``utils.py``, ``utilities/gis.py``) from RavenWPS to RavenPy.
* RavenPy can now be installed from PyPI without GIS dependencies (limited functionality).
* Hydro routing product is now supported from ``geoserver.py`` (a notebook has been added to demonstrate the new functions).
* New script ``ravenpy aggregate-forcings-to-hrus`` to aggregate NetCDF files and compute updated grid weights.
* Add the basis for a new routing emulator option (WIP).
* Add climpred verification capabilities.

0.2.3
-----

* Regionalisation data is now part of the package.
* Fix tests that were not using testdata properly.
* Add tests for external dataset access.
* ``utilities.testdata.get_local_testdata`` now raises an exception when it finds no dataset corresponding to the user pattern.

0.2.2
-----

* Set wcs.getCoverage timeout to 120 seconds.
* Fix ``Raven.parse_results`` logic when no flow observations are present and no diagnostic file is created.
* Fix ECCC test where input was cached and shadowed forecast input data.

0.2.1
-----

* Fix xarray caching bug in regionalization.

0.2.0
-----

* Refactoring of ``ravenpy.utilities.testdata`` functions.
* Bump xclim to 0.23.

0.1.7
-----

* Fix xarray caching bug affecting climatological ESP forecasts (#33).
* Fix deprecation issue with Fiona.

0.1.6 (2021-01-15)
------------------

* Correct installer bugs.

0.1.5 (2021-01-14)
------------------

* Release with docs.

0.1.0 (2020-12-20)
------------------

* First release on PyPI.
