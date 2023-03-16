import numpy as np

from ravenpy import Emulator

NSE = {"GR4JCN": -0.117301, "HMETS": -3.0132, "Mohyse": 0.194612}


def test_rv(emulator_rv):
    """Test the model configuration can be written to disk."""
    assert (emulator_rv / "test.rvi").exists()


def test_run(emulator, tmpdir):
    """Check the emulator actually runs."""
    name = emulator.__class__.__name__
    e = Emulator(config=emulator, path=tmpdir)
    e.build()
    e.run()
    out = e.parse()
    d = out["diagnostics"]
    np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], NSE[name], 4)
