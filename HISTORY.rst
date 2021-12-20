=======
History
=======

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

* Update cruft
* Subclass `derived_parameters` in Ostrich emulators to avoid having to pass `params`.

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
