import numpy as np

from ravenpy import Emulator

NSE = {"GR4JCN": -0.117301, "HMETS": -3.0132, "Mohyse": 0.194612, "HBVEC": 0.0186633}


def test_rv(config_rv):
    """Test the model configuration can be written to disk."""
    assert (config_rv / "test.rvi").exists()


def test_run(numeric_config, tmpdir):
    """Check the emulator actually runs."""
    name = numeric_config.__class__.__name__
    e = Emulator(config=numeric_config, path=tmpdir)
    e.build()
    e.run()
    out = e.parse()
    d = out["diagnostics"]
    np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], NSE[name], 4)
