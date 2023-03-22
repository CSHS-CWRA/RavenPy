import numpy as np

from ravenpy import Emulator

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


def test_run(numeric_config, tmpdir):
    """Test that the emulator actually runs and returns the expected NSE."""
    name, cls = numeric_config
    e = Emulator(config=cls, path=tmpdir)
    e.build()
    e.run()
    out = e.parse()
    d = out["diagnostics"]
    np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], NSE[name], 4)
