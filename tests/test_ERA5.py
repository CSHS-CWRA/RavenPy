import datetime as dt

from ravenpy.models import HMETS

SALMON_coords = (-123.3659, 54.4848)  # (lon, lat)

era5 = "era5/tas_pr_20180101-20180108.nc"


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


class TestRavenERA5:
    def test_simple(self, get_file):
        ts = get_file(era5)
        model = HMETS()
        model(
            ts=ts,
            params=params,
            start_date=dt.datetime(2018, 1, 1),
            end_date=dt.datetime(2018, 1, 3),
            run_name="test-hmets-era5",
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            rain_snow_fraction="RAINSNOW_DINGMAN",
            tas={"time_shift": -0.25},
            pr={"time_shift": -0.25},
        )
