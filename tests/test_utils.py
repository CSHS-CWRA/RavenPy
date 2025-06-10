import pytest
from packaging.version import Version
from pydap import __version__ as __pydap_version__

from ravenpy.config.utils import nc_specs

older_pydap = False
if Version(__pydap_version__) < Version("3.5.5"):
    older_pydap = True


def test_nc_specs(yangtze):
    f = yangtze.fetch(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    attrs = nc_specs(f, "PRECIP", station_idx=1, alt_names=("rain",))
    assert "file_name_nc" in attrs


def test_nc_specs_bad(bad_netcdf):
    # Test with scalar elevation. Should normally have a station dimension, but not always the case.
    attrs_s = nc_specs(bad_netcdf, "PRECIP", station_idx=1, alt_names=("rain",))
    assert attrs_s["elevation"] == 1.0


@pytest.mark.online
@pytest.mark.skipif(
    older_pydap, reason="pydap version 3.5.5 is required for this test", strict=False
)
def test_dap_specs():
    # Link to THREDDS Data Server netCDF testdata
    tds = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/testdata/raven"
    fn = f"{tds}/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"

    attrs = nc_specs(fn, "PRECIP", station_idx=1, alt_names=("rain",), engine="pydap")
    assert "units" in attrs
