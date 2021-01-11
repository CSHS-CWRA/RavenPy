import datetime as dt
import json
import tempfile

import xarray as xr

from ravenpy.models import HMETS
from ravenpy.utilities.testdata import get_test_data

# Get path to ncml file for NRCan data.
NRCAN_path = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/gridded_obs/nrcan_v2.ncml"

# Temporary path
filepath = tempfile.mkdtemp() + "/NRCAN_ts.nc"

# Get information for given catchment, could be passed in parameter to the function
ts = get_test_data(
    "raven-gr4j-cemaneige", "Salmon-River-Near-Prince-George_meteo_daily.nc"
)[0]
salmon = xr.open_dataset(ts)
lat = salmon.lat.values[0]
lon = salmon.lon.values[0]

# Start and end dates
start_date = dt.datetime(2006, 1, 1)
end_date = dt.datetime(2007, 12, 31)

# Get data for the covered period including a 1-degree bounding box for good measure. Eventually we will be
# able to take the catchment polygon as a mask and average points residing inside.

ds = (
    xr.open_dataset(NRCAN_path)
    .sel(
        lat=slice(lat + 1, lat - 1),
        lon=slice(lon - 1, lon + 1),
        time=slice(start_date, end_date),
    )
    .mean(dim={"lat", "lon"}, keep_attrs=True)
)
ds.to_netcdf(filepath)


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

        model = HMETS()
        model(
            ts=filepath,
            params=params,
            start_date=start_date,
            end_date=end_date,
            name="Salmon",
            run_name="test-hmets-NRCAN",
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            rain_snow_fraction="RAINSNOW_DINGMAN",
            tasmax={"linear_transform": (1.0, -273.15)},
            tasmin={"linear_transform": (1.0, -273.15)},
            pr={"linear_transform": (86400, 0.0)},
        )
