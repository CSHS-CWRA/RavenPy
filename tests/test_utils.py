import pytest
from packaging.version import Version
from pydap import __version__ as __pydap_version__

from ravenpy.config.utils import nc_specs


older_pydap = False
if Version(__pydap_version__) < Version("3.5.5"):
    older_pydap = True


def test_nc_specs(yangtze):
    f = yangtze.fetch("raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc")
    attrs = nc_specs(f, "PRECIP", station_idx=1, alt_names=("rain",))
    assert "file_name_nc" in attrs
    assert attrs["dim_names_nc"] == ("time",)

    # 2D with station dimension
    f = yangtze.fetch("raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc")
    attrs = nc_specs(f, "PRECIP", station_idx=1, alt_names=("rain",))
    assert attrs["dim_names_nc"] == (
        "region",
        "time",
    )

    # 3D - Since this file is not CF compliant, nc_specs cannot infer the correct dimension order
    f = yangtze.fetch("raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc")
    attrs = nc_specs(f, "PRECIP", station_idx=1, alt_names=("rain",))
    assert attrs["dim_names_nc"] == ("time", "lon", "lat")

    f = yangtze.fetch("cmip5/tas_Amon_CanESM2_rcp85_r1i1p1_200601-210012_subset.nc")
    attrs = nc_specs(f, "TEMP_AVE", station_idx=1, engine="netcdf4")
    assert attrs["dim_names_nc"] == ("lon", "lat", "time")


def test_nc_specs_bad(bad_netcdf):
    # Test with scalar elevation. Should normally have a station dimension, but not always the case.
    attrs_s = nc_specs(bad_netcdf, "PRECIP", station_idx=1, alt_names=("rain",))
    assert attrs_s["elevation"] == 1.0


@pytest.mark.online
def test_dap_specs():
    # Link to THREDDS Data Server netCDF testdata
    tds = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/testdata/raven"
    fn = f"{tds}/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"

    attrs = nc_specs(fn, "PRECIP", station_idx=1, alt_names=("rain",), engine="pydap")
    assert "units" in attrs
