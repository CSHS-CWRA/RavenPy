import pytest
import spotpy

from ravenpy.utilities.calibration import SpotSetup

bounds = {
    "GR4JCN": dict(
        low=(0.01, -15.0, 10.0, 0.0, 1.0, 0.0),
        high=(2.5, 10.0, 700.0, 7.0, 30.0, 1.0),
    ),
    "HMETS": dict(
        low=(
            0.3,
            0.01,
            0.5,
            0.15,
            0.0,
            0.0,
            -2.0,
            0.01,
            0.0,
            0.01,
            0.005,
            -5.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.00001,
            0.0,
            0.00001,
            0.0,
            0.0,
        ),
        high=(
            20.0,
            5.0,
            13.0,
            1.5,
            20.0,
            20.0,
            3.0,
            0.2,
            0.1,
            0.3,
            0.1,
            2.0,
            5.0,
            1.0,
            3.0,
            1.0,
            0.02,
            0.1,
            0.01,
            0.5,
            2.0,
        ),
    ),
    "Mohyse": dict(
        low=(0.01, 0.01, 0.01, -5.00, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
        high=(20.0, 1.0, 20.0, 5.0, 0.5, 1.0, 1.0, 1.0, 15.0, 15.0),
    ),
    "HBVEC": dict(
        low=(
            -3.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.3,
            0.0,
            0.0,
            0.01,
            0.05,
            0.01,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.01,
            0.0,
            0.05,
            0.8,
            0.8,
        ),
        high=(
            3.0,
            8.0,
            8.0,
            0.1,
            1.0,
            1.0,
            7.0,
            100.0,
            1.0,
            0.1,
            6.0,
            5.0,
            5.0,
            0.2,
            1.0,
            30.0,
            3.0,
            2.0,
            1.0,
            1.5,
            1.5,
        ),
    ),
}


def test_spotpy_calibration(symbolic_config, tmpdir):
    name = symbolic_config.__class__.__name__
    if name not in bounds:
        pytest.skip("No bounds defined.")

    spot_setup = SpotSetup(
        config=symbolic_config,
        low=bounds[name]["low"],
        high=bounds[name]["high"],
        path=tmpdir,
    )

    sampler = spotpy.algorithms.dds(
        spot_setup, dbname="RAVEN_model_run", dbformat="ram", save_sim=False
    )

    sampler.sample(10, trials=1)
    d1 = sampler.status.objectivefunction_max
    sampler.sample(50, trials=1)
    d2 = sampler.status.objectivefunction_max
    assert d2 > d1
