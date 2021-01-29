=======
History
=======

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
