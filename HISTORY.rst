=======
History
=======

0.4.0
-----

* GR4JCN emulator now supports routing mode.
* Add BLENDED emulator.

0.3.1
-----

* Update external dependencies (Raven, OSTRICH) to facilitate Conda packaging.

0.3.0
-----

* Migration and refactoring of GIS and IO utilities (`utils.py`, `utilities/gis.py`) from RavenWPS to RavenPy.
* RavenPy can now be installed from PyPI without GIS dependencies (limited functionality).
* Hydro routing product is now supported from `geoserver.py` (a notebook has been added to demonstrate the new functions).
* New script `ravenpy aggregate-forcings-to-hrus` to aggregate NetCDF files and compute updated grid weights.
* Add the basis for a new routing emulator option (WIP).
* Add climpred verification capabilities.

0.2.3
-----

* Regionalisation data is now part of the package.
* Fix tests that were not using testdata properly.
* Add tests for external dataset access.
* `utilities.testdata.get_local_testdata` now raises an exception when it finds no dataset corresponding to the user pattern.

0.2.2
-----

* Set wcs.getCoverage timeout to 120 seconds.
* Fix `Raven.parse_results` logic when no flow observations are present and no diagnostic file is created.
* Fix ECCC test where input was cached and shadowed forecast input data.

0.2.1
-----

* Fix xarray caching bug in regionalization.

0.2.0
-----

* Refactoring of `ravenpy.utilities.testdata` functions.
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
