import pytest

from ravenpy.config.utils import nc_specs


def test_nc_specs(get_local_testdata):
    f = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    attrs = nc_specs(f, "PRECIP", station_idx=1, alt_names=("rain",))
    assert "file_name_nc" in attrs


def test_nc_specs_bad(bad_netcdf):
    # Test with scalar elevation. Should normally have a station dimension, but not always the case.
    attrs_s = nc_specs(bad_netcdf, "PRECIP", station_idx=1, alt_names=("rain",))
    assert attrs_s["elevation"] == 1.0


@pytest.mark.online
def test_dap_specs():
    # Link to THREDDS Data Server netCDF testdata
    TDS = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/testdata/raven"
    fn = f"{TDS}/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"

    attrs = nc_specs(fn, "PRECIP", station_idx=1, alt_names=("rain",), engine="pydap")
    assert "units" in attrs
