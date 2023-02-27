=======
History
=======

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
