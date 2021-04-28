import datetime as dt

from ravenpy.models import HMETS
from ravenpy.utilities.testdata import get_local_testdata


class TestRavenNRCAN:
    def test_simple(self):

        params = (
            9.5019,
            0.2774,
            6.3942,
            0.6884,
            1.2875,
            5.4134,
            2.3641,
            0.0973,
            0.0464,
            0.1998,
            0.0222,
            -1.0919,
            2.6851,
            0.3740,
            1.0000,
            0.4739,
            0.0114,
            0.0243,
            0.0069,
            310.7211,
            916.1947,
        )

        ts = get_local_testdata("nrcan/NRCAN_2006-2007_subset.nc")
        start_date = dt.datetime(2006, 1, 1)
        end_date = dt.datetime(2007, 12, 31)

        model = HMETS()
        model(
            ts=ts,
            params=params,
            start_date=start_date,
            end_date=end_date,
            run_name="test-hmets-NRCAN",
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            rain_snow_fraction="RAINSNOW_DINGMAN",
            tasmax={"scale": 1.0, "offset": -273.15},
            tasmin={"scale": 1.0, "offset": -273.15},
            pr={"scale": 86400, "offset": 0.0},
        )
