import datetime as dt

import numpy as np
import pytest

from ravenpy import Emulator
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
    e = Emulator(config=conf, workdir=tmp_path)
    e.write_rv()
    out = e.run()

    d = out.diagnostics
    np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], NSE[name], 4)

    if name == "CanadianShield":
        pytest.skip("Missing solution due to SuppressOutput.")

    new = e.resume()
    new.start_date = conf.end_date
    new.end_date = dt.datetime(2002, 1, 7)

    e2 = Emulator(config=new, workdir=tmp_path / "resumed")
    e2.run()


def test_emulator(dummy_config, tmp_path):
    cls, P = dummy_config

    e = Emulator(config=cls(), workdir=tmp_path)

    # The emulator configuration should be read-only.
    with pytest.raises(TypeError):
        e.config.run_name = "Renamed"
