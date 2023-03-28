import datetime as dt

import numpy as np
import pytest
import xarray as xr

from ravenpy import Emulator
from ravenpy.new_config import commands as rc
from ravenpy.new_config.rvs import Config

# Expected NSE for emulator configuration from the `config_rv` test fixture.
NSE = {
    "GR4JCN": -0.117301,
    "HMETS": -3.0132,
    "Mohyse": 0.194612,
    "HBVEC": 0.0186633,
    "CanadianShield": 0.39602,
    "HYPR": 0.685188,
    "SACSMA": -0.0382907,
    "Blended": -0.913785,
}


def test_rv(config_rv):
    """Test the model configuration can be written to disk."""
    name, path = config_rv
    assert (path / "test.rvi").exists()


def test_run(numeric_config, tmp_path):
    """Test that the emulator actually runs and returns the expected NSE."""
    name, conf = numeric_config
    assert conf.__config__.allow_mutation

    e = Emulator(config=conf, workdir=tmp_path)
    e.write_rv()
    out = e.run()

    d = out.diagnostics
    np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], NSE[name], 4)

    if name == "CanadianShield":
        pytest.skip("Missing solution due to SuppressOutput.")

    # Start new simulation with final state from initial run.
    new = e.resume()
    new.start_date = conf.end_date
    new.end_date = dt.datetime(2002, 1, 7)
    new.run_name = None

    e2 = Emulator(config=new, workdir=tmp_path / "resumed")
    out = e2.run()
    assert e2.modelname == "raven"
    assert isinstance(out.hydrograph, xr.Dataset)
    assert isinstance(out.storage, xr.Dataset)


@pytest.mark.skip("Need to find a clean way to freeze emulator config instance.")
def test_emulator_config_is_read_only(dummy_config, tmp_path):
    cls, P = dummy_config

    e = Emulator(config=cls(), workdir=tmp_path)

    # The emulator configuration should be read-only.
    with pytest.raises(TypeError):
        e.config.run_name = "Renamed"


def test_no_run_name(dummy_config, tmp_path):
    cls, P = dummy_config
    conf = cls(params=[0.5])

    paths = conf.write_rv(tmp_path, modelname="ham")
    assert "ham.rvi" in str(paths["rvi"])

    paths = conf.write_rv(tmp_path, overwrite=True)
    assert "raven.rvi" in str(paths["rvi"])


@pytest.mark.slow
@pytest.mark.online
@pytest.mark.xfail(raises=OSError, reason="pavics.ouranos.ca is unreachable")
def test_run_with_dap_link(minimal_emulator, tmp_path):
    """Test Raven with DAP link instead of local netCDF file."""
    # Link to THREDDS Data Server netCDF testdata
    TDS = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/testdata/raven"
    fn = f"{TDS}/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"

    alt_names = {
        "RAINFALL": "rain",
        "TEMP_MIN": "tmin",
        "TEMP_MAX": "tmax",
        "PET": "pet",
        "HYDROGRAPH": "qobs",
        "SNOWFALL": "snow",
    }

    conf = minimal_emulator
    conf.gauge = rc.Gauge.from_nc(fn, alt_names=alt_names)

    out = Emulator(conf, workdir=tmp_path).run()
