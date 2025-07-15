=========
Changelog
=========

v0.19.0 (2025-07-16)
--------------------

New features
^^^^^^^^^^^^
* Added `parsers.parse_rv` to extract a Command value from an RV file. (PR #503)
* New module `ravenpy.testing` has been added to provide utility functions and support for testing and testing data management. (PR #513)

Breaking changes
^^^^^^^^^^^^^^^^
* `ravenpy` now requires `pooch>=1.8.0` for downloading and caching remote testing data. (PR #513)
* `ravenpy.utilities.testdata` has been refactored to new module `ravenpy.testing`. The `publish_release_notes` function is now located in `ravenpy.utilities.publishing`. (PR #513)
* The `ravenpy.testing.utils` module now provides a `yangtze()` class for fetching and caching the `raven-testdata` testing data. A convenience function (`get_file`) replaces the previous `get_local_testdata`. (PR #513)
* The `ravenpy.testing.utils.open_dataset` function no longer supports OPeNDAP URLs or local file paths. Instead, it uses the `yangtze()` class to fetch datasets from the testing data repository or the local cache. Users should now use `xarray.open_dataset()` directly for OPeNDAP URLs or local files. (PR #513)

Bug fixes
^^^^^^^^^
* Fixed bug affecting `GriddedForcing`, where `station_idx` in the call to ``nc_specs`` was set to ``1`` instead of ``None``. (PR #501)
* Fixed bug in the `run` method, where `overwrite=False` was not being respected. (PR #503)
* Pin `scipy` below v1.16.0 due to a breaking change that affects `statsmodels` below v0.14.4. (PR #521)

Internal changes
^^^^^^^^^^^^^^^^
* `ravenpy` now requires `xclim>=0.57.0` and `xsdba` (v0.4.0+). (PR #511)
* The `tests` folder no longer contains an `__init__.py` file and is no longer treated as a package. `pytest` fixtures from `emulators.py` are now directly imported into `conftest.py` for use in tests, and existing `pytest` fixtures have been modified to use the new `yangtze()` class for fetching testing data. (PR #513)

v0.18.2 (2025-05-05)
--------------------

New features
^^^^^^^^^^^^
* Added `RelativeHumidityMethod` to RVIs. (PR #490).
* Add `parse` methods to `LineCommand`, `SubBasins`, `HRUs`, `Reservoir`, `SubBasinGroup`, `ChannelProfile`. (PR #492).
* Tweak the `GridWeightExtractor` to support datasets from the Canadian River and Lake Hydrofabric. Allows setting the `routing_id_field` to `__INDEX__` in order to match HRU IDs. (PR #492).
* Add support for `Recharge` process. (PR #492).

v0.18.1 (2025-04-15)
--------------------

New features
^^^^^^^^^^^^
* `ravenpy` no longer installs `raven-hydro` by default. The Raven model executable can now be provided by explicitly setting the `RAVENPY_RAVEN_BINARY_PATH` environment variable. (PR #486).

Bug fixes
^^^^^^^^^
* Fixed a bug in ``ravenpy.utilities.regionalization.multiple_linear_regression`` that was calling a class method incorrectly. (PR #484).

Internal changes
^^^^^^^^^^^^^^^^
* `pydap` has been pinned below v3.5.5 temporarily until `xarray` offers support for it. (PR #486).
* More than 7500 DeprecationWarnings emitted during the testing suite have been addressed. Minimum supported `pydantic` has been raised to v2.11. (PR #487).
* Regenerated the notebook outputs using newer version of `xclim`. (PR #484).

v0.18.0 (2025-04-03)
--------------------

New features
^^^^^^^^^^^^
* `ravenpy` now supports Python3.13. (PR #459)
* Updated `raven-hydro` to v0.4.0 (`RavenHydroFramework` v4.0.1). (PR #459)
* Updated `xclim` to v0.54.0, `pint` to v0.24.4, and `numpy` to v1.24.0 (no longer pinned below v2.0). (PR #459)
* `ravenpy` is now registered with the Open Source Security Foundation (OSSF) Best Practices initiative (`RavenPy OpenSSF-BP Status <https://www.bestpractices.dev/en/projects/10064>`_). (PR #464)
* `ravenpy` now enables new EvaluationMetrics commands in the model configuration. Other features from `RavenHydroFramework` will be included in newer releases. (PR #476)

Bug fixes
^^^^^^^^^
* Fix bug in _MonthlyRecord class definition crashing the pydantic-autodoc serialization. (PR #458)
* Fixed a small API bug in the `Comparing_hindcasts_and_ESP_forecasts.ipynb` notebook. (PR #463)
* The `Raven` model previously always reported version "3.7", regardless of the installed `Raven` version. It now uses `raven-hydro`'s `__raven_version__` attribute. (PR #464)

Internal changes
^^^^^^^^^^^^^^^^
* Updated the cookiecutter template to the latest commit: (PR #454)
    * GitHub Actions and Python dependencies have been updated.
    * New `pre-commit` hooks for `vulture` (find dead code) and `codespell` (spelling errors).
    * Removed several `type: ignore` statements.
    * Spelling errors in documentation have been addressed.
* GitHub Workflows now test `ravenpy` using macOS as well as Python3.13. (PR #459)
* Several small deprecation and usage warnings as well as a few variable typing issues have been addressed. (PR #464)
* Updated the license to reflect current year. (PR #476)
* Documentation version now supports showing hyphens in the version number. (PR #476)
* Call signatures and docstrings of functions have been modified to be more precise for the expected variable type. (PR #476)

v0.17.0 (2025-01-27)
--------------------

* Updated the cookiecutter template to the latest commit and synchronized dependencies between PyPI and Anaconda recipes. (PR #427)
* Updated `ts_fit_graph` logic for `matplotlib` >= 3.10.0 compatibility. (PR #434)
* Temporarily pinned `pygments` below v2.19 due to a breaking change affecting `sphinx-codeautolink`. (PR #434)
* Adopted a new RavenPy logo for the documentation. (PR #428)
* Documentation Updates: (PR #436)
    * Cleaner imports, removed some unneeded library imports.
    * Typo and grammar fixes.
    * Updated the Python, Anaconda, and Ubuntu versions used to generate the documentation.
* Small import fixes and minor code cleanup (`ravenpy.extractors`). (PR #436)
* Adjusted pins for `intake`, `intake-esm` and `zarr` to ensure notebooks run correctly. (PR #440)
* Added a Security Policy (`SECURITY.md`) to the repository. (PR #441)
* Updated the cookiecutter template to the latest commit: (PR #444)
    * Workflows now use Checkout with `persist-credentials: false`.
    * CodeQL workflow has been updated.
    * `pre-commit` hooks for `vulture` (finding dead code) and `zizmor` (finding security vulnerabilities) have been added.

v0.16.1 (2024-12-05)
--------------------

* Improved the HBV-EC emulator by reformatting some of the content (cosmetic changes), adding the spatial interpolation and the soil parameters for the TOPSOIL explicitly. (PR #410)
* Add support for `pymbolic` > 2022.2 (PR #420)
* Convert `SymConfig` to a dict to silence pydantic deprecation warnings (PR #420)
* Drop support for Python 3.9 (PR #420)

v0.16.0 (2024-10-18)
--------------------

* Set base required `geopandas` to v1.0. (PR #394)
* Removed the pin on `pyogrio` (set by `geopandas` now). (PR #394)
* Removed the `requests` dependency (now using `urllib`/`urllib3`). (PR #394)

Internal changes
^^^^^^^^^^^^^^^^
* The cookiecutter template has been updated to the latest commit: (PR #386)
    * `ravenpy` now uses a `src`-layout for the package.
    * `HISTORY.rst` has been renamed to `CHANGELOG.rst`.
    * `ruff` checks have replaced most of the `flake8` checks.
    * `ravenpy` now has a `CODE_OF_CONDUCT.md` file.
    * Many `numpydoc`-style docstrings have been adjusted for consistency.
* Added `setuptools` to the `gis` build recipe to ensure that the `gdal` bindings are built successfully. (PR #400)
* Modified the sub-basin and channel profile extraction functions to correctly set the river length to zero and set default values for reach attributes in sub-basins with no channel routing (i.e., sub-basins with lakes or headwater basins). (issue #354, PR #401)
* Improved the HBV-EC emulator by adding parameter information (name, definition, and Raven default values), fixed the variable name for the adiabatic temperature lapse rate, and added an alias for rain snow fraction to match other emulators. (PR #404 and #408)
* Modified the `sphinx` configuration to better support SVG and to remove incompatible elements from the PDF build. (PR #407)

v0.15.0 (2024-06-20)
--------------------

* Pinned `pint` below version 0.24 due to a breaking change in their API. (PR #375)
* Pinned `numpy` below v2.0.0 due to a breaking change in their API. (PR #378)
* Update `raven-hydro` to v0.3.1 and `RavenHydroFramework` to v3.8.1. (PR #378)
* Fixed bug in `Config.duplicate` dating from the switch to Pydantic V2 in 0.13 (PR #367)

Internal changes
^^^^^^^^^^^^^^^^
* Synchronize several dependencies between `pyproject.toml`, `environment*.yml`, and `tox.ini`. (PR #378)
* Drop the code formatting conventions for Python3.8, extend to Python3.11 and Python3.12. (PR #378)
* Addresses a bunch of small warnings in the pytest output. (PR #378)

v0.14.1 (2024-05-07)
--------------------

* Upgraded `owslib` to `>=0.29.1`. (PR #358)
* All operations that open NetCDF files or DAP links accept an `engine` argument. The default for all of these is `h5netcdf`. (PR #358)
* Added `pydap` as an alternate backend for opening DAP links. (PR #358)
* Fixed buggy CustomOutput command. (PR #360)
* Make sure config and output paths are absolute. (PR #360)

Internal changes
^^^^^^^^^^^^^^^^
* Added some development dependencies that were missing to the `environment.yml`. (PR #358)
* `test_climpred_hindcast_verif` is now skipped for Python3.10 builds. It seems to only fail on the particular version of Python. When examining the dependencies, other than the Python version (and ABI version), there are no differences in the environments between Python3.10 and Python3.11. Possibly an issue with `climpred`. (PR #358)
* Temporarily disabled tests for macOS on GitHub due to architecture changes. (PR #358)
* Pinned `pyogrio` below v0.8.0 until `geopandas` supports it. (PR #363)
* Updated linting dependencies to the latest versions. (PR #363)

v0.14.0 (2024-03-13)
--------------------

* Add support for new processes and methods added in Raven v3.8. (PR #335)
* Add Interpolation command options. (PR #338)
* Let VegetationClass records contain symbolic expressions. (PR #338)
* Add support for custom RV subclasses. (PR #338)
* Use HRU_ID (if available) instead of SubId in BasinMaker reservoirs extraction logic. (PR #338)
* Added support for Python 3.12 and dropped support for Python3.8. (PR #341, PR #343)
* Added support for `raven-hydro` v0.3.0 and `RavenHydroFramework` to v3.8. (PR #341, PR #351)
* `ravenpy` now requires `xclim` >= v0.48.2, `xarray` >= v2023.11.0, and `pandas` >= 2.2.0. (PR #341)
* Now automatically filters HRUs based on the ``hru_type``. (issue #340, PR #334)

Internal changes
^^^^^^^^^^^^^^^^
* Updated GitHub publishing workflows to use Trusted Publisher for TestPyPI/PyPI releases. (PR #341)
* Added Dependabot to keep dependencies up-to-date. (PR #342)
* Now using step-security/harden-runner action to harden GitHub Actions runners. (PR #341)
* Adjusted GitHub Workflows to test against Python 3.9, 3.10, 3.11, and 3.12. (PR #341, PR #343)
* Updated the build-system requirements when testing with `tox` to use newer `setuptools` and `wheel` versions when building `gdal`. (PR #341)

v0.13.0 (2024-01-10)
--------------------

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

v0.12.3 (2023-10-02)
--------------------

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

v0.12.2 (2023-07-04)
--------------------

This release is primarily a bugfix to address issues arising from dependencies.

Breaking changes
^^^^^^^^^^^^^^^^
* `raven-hydro` version has been bumped from v0.2.1 to v0.2.3. This version provides better support for builds on Windows and MacOS.
* Due to major breaking changes, `pydantic` has been pinned below v2.0 until changes can be made to adapt to their new API.
* `numpy` has been pinned below v1.25.0 to ensure compatibility with `numba`.

Internal changes
^^^^^^^^^^^^^^^^
* ``test_geoserver::test_select_hybas_ar_domain_point`` is now temporarily skipped when testing on MacOS due to a mysterious domain identification error.

v0.12.1 (2023-06-01)
--------------------

This release is largely a bugfix to better stabilize performance and enhance the documentation.

* Avoid repeatedly calling `xr.open_dataset` in `OutputReader`'s `hydrograph` and `storage` properties. This seems to cause kernel failures in Jupyter notebooks.

Internal changes
^^^^^^^^^^^^^^^^
* Hyperlinks to documented functions now points to entries in the `User API` section.
* Docstrings are now more conformant to numpy-docstring conventions and formatting errors raised from badly-formatted pydantic-style docstrings have been addressed.
* In order to prevent timeout and excessive memory usage, Jupyter notebooks have been adjusted to no longer run on ReadTheDocs. All notebooks have been updated to the latest RavenPy and remain tested against RavenPy externally.
* Documentation built on ReadTheDocs is now set to `fail_on_warning`.

v0.12.0 (2023-05-25)
--------------------

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

v0.11.0 (2023-02-16)
--------------------

* Update RavenC executable to v3.6.
* Update xclim library to v0.40.0.
* Update fiona library to v1.9.
* Address some failures that can be caused when attempting to run CLI commands without the proper GIS dependencies installed.
* Addressed warnings raised in conda-forge compilation due to badly-configured MANIFEST.in.
* Update installation documentation to reflect most recent changes.

v0.10.0 (2022-12-21)
--------------------

* Update Raven executable to 3.5. Due to a bug in RavenC, simulations storing reservoir information to netCDF will fail. We expect this to be resolved in the next release. Note that we only test RavenPy with one Raven version. There is no guarantee it will work with other versions.
* Relax geo test to avoid failures occurring due to GDAL 3.6.
* Pin numpy below 1.24 (see https://github.com/numba/numba/issues/8615)

v0.9.0 (2022-11-16)
-------------------

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

v0.8.1 (2022-10-26)
-------------------

* Undo change related to `suppress_output`, as it breaks multiple tests in raven. New `Raven._execute` method runs models but does not parse results.

v0.8.0 (2022-10-04)
-------------------

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

v0.7.8 (2022-01-14)
-------------------

* Added functionalities in Data Assimilation utils and simplified tests.
* Removed pin on setuptools.
* Fixed issues related to symlinks, working directory, and output filenames.
* Fixed issues related to GDAL version handling in conda-forge.
* Updated jupyter notebooks.

0.7.7 (2021-12-21)
------------------

* Updated internal shapely calls to remove deprecated ``.to_wkt()`` methods.

0.7.6 (2021-12-20)
------------------

* Automate release pipeline to PyPI using GitHub CI actions.
* Added coverage monitoring GitHub CI action.
* Various documentation adjustments.
* Various metadata adjustments.
* Pinned owslib to 0.24.1 and above.
* Circumvented a bug in GitHub CI that was causing tests to fail at collection stage.

v0.7.5 (2021-09-10)
-------------------

* Update test so that it works with xclim 0.29.

v0.7.4 (2021-09-02)
-------------------

* Pinned climpred below v2.1.6.

v0.7.3 (2021-08-31)
-------------------

* Pinned xclim below v0.29.

v0.7.2 (2021-08-31)
-------------------

* Update cruft.
* Subclass ``derived_parameters`` in Ostrich emulators to avoid having to pass ``params``.

v0.7.0 (2021-07-27)
-------------------

* Add support for V2.1 of the Routing Product in ``ravenpy.extractors.routing_product``.
* Add ``collect-subbasins-upstream-of-gauge`` CLI script.
* Modify WFS request functions to use spatial filtering (``Intersects``) supplied by OWSLib.

v0.6.0 (2021-07-14)
-------------------

* Add support for EvaluationPeriod commands. Note that as a result of this, the model's ``diagnostics`` property contains one list per key, instead of a single scalar. Also note that for calibration, Ostrich will use the first period and the first evaluation metric.
* Add ``SACSMA``, ``CANADIANSHIELD`` and ``HYPR`` model emulators.

v0.5.2 (2021-05-25)
-------------------

* Simplify RVC configuration logic.
* Add ``ravenpy.utilities.testdata.file_md5_checksum`` (previously in ``xarray.tutorial``).

v0.5.1 (2021-05-12)
-------------------

* Some adjustments and bugfixes needed for RavenWPS.
* Refactoring of some internal logic in ``ravenpy.config.rvs.RVT``.
* Improvements to typing with the help of mypy.

v0.5.0 (2021-04-30)
-------------------

* Refactoring of the RV config subsystem:

  * The config is fully encapsulated into its own class: ``ravenpy.config.rvs.Config``.
  * The emulator RV templates are inline in their emulator classes.

* The emulators have their own submodule: ``ravenpy.models.emulators``.
* The "importers" have been renamed to "extractors" and they have their own submodule: ``ravenpy.extractors``.

v0.4.2 (2021-04-14)
-------------------

* Update to RavenC revision 318 to fix OPeNDAP access for StationForcing commands.
* Fix grid_weights set to None by default.
* Pass nc_index to ObservationData command.
* Expose more cleanly RavenC errors and warnings.

v0.4.1 (2021-04-13)
-------------------

* Add notebook about hindcast verification skill.
* Add notebook about routing capability.
* Modify geoserver functions to have them return GeoJSON instead of GML.
* Collect upstream watershed aggregation logic.
* Fix RVC bug.

v0.4.0 (2021-04-09)
-------------------

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

v0.3.1 (2021-04-06)
-------------------

* Update external dependencies (Raven, OSTRICH) to facilitate Conda packaging.

v0.3.0 (2021-03-11)
-------------------

* Migration and refactoring of GIS and IO utilities (``utils.py``, ``utilities/gis.py``) from RavenWPS to RavenPy.
* RavenPy can now be installed from PyPI without GIS dependencies (limited functionality).
* Hydro routing product is now supported from ``geoserver.py`` (a notebook has been added to demonstrate the new functions).
* New script ``ravenpy aggregate-forcings-to-hrus`` to aggregate NetCDF files and compute updated grid weights.
* Add the basis for a new routing emulator option (WIP).
* Add climpred verification capabilities.

v0.2.3 (2021-02-01)
-------------------

* Regionalisation data is now part of the package.
* Fix tests that were not using testdata properly.
* Add tests for external dataset access.
* ``utilities.testdata.get_local_testdata`` now raises an exception when it finds no dataset corresponding to the user pattern.

v0.2.2 (2021-01-29)
-------------------

* Set wcs.getCoverage timeout to 120 seconds.
* Fix ``Raven.parse_results`` logic when no flow observations are present and no diagnostic file is created.
* Fix ECCC test where input was cached and shadowed forecast input data.

v0.2.1 (2021-01-28)
-------------------

* Fix xarray caching bug in regionalization.

v0.2.0 (2021-01-26)
-------------------

* Refactoring of ``ravenpy.utilities.testdata`` functions.
* Bump xclim to 0.23.

v0.1.7 (2021-01-19)
-------------------

* Fix xarray caching bug affecting climatological ESP forecasts (#33).
* Fix deprecation issue with Fiona.

v0.1.6 (2021-01-15)
-------------------

* Correct installer bugs.

v0.1.5 (2021-01-14)
-------------------

* Release with docs.

v0.1.0 (2020-12-20)
-------------------

* First release on PyPI.
