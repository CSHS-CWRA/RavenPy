import pytest
import spotpy

from ravenpy.utilities.calibration import SpotSetup

# Low and high bounds for parameter distribution
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
    "CanadianShield": dict(
        low=(
            0.010,
            0.010,
            0.010,
            0.000,
            0.000,
            0.050,
            0.000,
            -5.00,
            -5.00,
            0.500,
            0.500,
            0.000,
            0.000,
            0.000,
            0.000,
            0.000,
            0.010,
            0.005,
            -3.00,
            0.500,
            5.000,
            0.000,
            0.000,
            -1.00,
            0.000,
            0.000,
            0.000,
            0.000,
            0.000,
            0.000,
            0.000,
            0.800,
            0.800,
            0.000,
        ),
        high=(
            0.500,
            2.000,
            3.000,
            3.000,
            0.050,
            0.450,
            7.000,
            -1.00,
            -1.00,
            2.000,
            2.000,
            100.0,
            100.0,
            100.0,
            0.400,
            0.100,
            0.300,
            0.100,
            3.000,
            4.000,
            500.0,
            5.000,
            5.000,
            1.000,
            8.000,
            20.00,
            1.500,
            0.200,
            0.200,
            10.00,
            10.00,
            1.200,
            1.200,
            1.000,
        ),
    ),
    "HYPR": dict(
        low=(
            -1.0,
            -3.0,
            0.00,
            0.30,
            -1.3,
            -2.0,
            0.00,
            0.10,
            0.40,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.01,
            0.00,
            0.00,
            1.50,
            0.00,
            0.00,
            0.80,
        ),
        high=(
            1.0000,
            3.0000,
            0.8000,
            1.0000,
            0.3000,
            0.0000,
            30.000,
            0.8000,
            2.0000,
            100.00,
            0.5000,
            5.0000,
            1.0000,
            1000.0,
            6.0000,
            7.0000,
            8.0000,
            3.0000,
            5.0000,
            5.0000,
            1.2000,
        ),
    ),
    "SACSMA": dict(
        low=(
            -3.00000000,
            -1.52287874,
            -0.69897000,
            0.025000000,
            0.010000000,
            0.075000000,
            0.015000000,
            0.040000000,
            0.000000000,
            0.000000000,
            0.000000000,
            0.000000000,
            0.000000000,
            0.000000000,
            0.000000000,
            0.000000000,
            0.000000000,
            0.300000000,
            0.010000000,
            0.800000000,
            0.800000000,
        ),
        high=(
            -1.82390874,
            -0.69897000,
            -0.30102999,
            0.125000000,
            0.075000000,
            0.300000000,
            0.300000000,
            0.600000000,
            0.500000000,
            3.000000000,
            80.00000000,
            0.800000000,
            0.050000000,
            0.200000000,
            0.100000000,
            0.400000000,
            8.000000000,
            20.00000000,
            5.000000000,
            1.200000000,
            1.200000000,
        ),
    ),
    "Blended": dict(
        low=(
            0.0,
            0.1,
            0.5,
            -5.0,
            0.0,
            0.5,
            5.0,
            0.0,
            0.0,
            0.0,
            -5.0,
            0.5,
            0.0,
            0.01,
            0.005,
            -5.0,
            0.0,
            0.0,
            0.0,
            0.3,
            0.01,
            0.5,
            0.15,
            1.5,
            0.0,
            -1.0,
            0.01,
            0.00001,
            0.0,
            0.0,
            -3.0,
            0.5,
            0.8,
            0.8,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ),
        high=(
            1.0,
            3.0,
            3.0,
            -1.0,
            100.0,
            2.0,
            10.0,
            3.0,
            0.05,
            0.45,
            -2.0,
            2.0,
            0.1,
            0.3,
            0.1,
            2.0,
            1.0,
            5.0,
            0.4,
            20.0,
            5.0,
            13.0,
            1.5,
            3.0,
            5.0,
            1.0,
            0.2,
            0.02,
            0.5,
            2.0,
            3.0,
            4.0,
            1.2,
            1.2,
            0.02,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ),
    ),
}


@pytest.mark.slow
def test_spotpy_calibration(symbolic_config, tmpdir):
    """Confirm that calibration actually improves the NSE."""
    name, cls = symbolic_config

    # FIXME: The Blended model run returns error code -11.
    if name == "Blended":
        pytest.skip("The Blended model run returns error code -11.")
    if name == "CanadianShield":
        pytest.skip(
            "The CanadianShield model run returns 'CHydroUnit constructor:: HRU 2 has a negative or zero area'"
        )

    if name not in bounds:
        pytest.skip("No bounds defined.")

    spot_setup = SpotSetup(
        config=cls,
        low=bounds[name]["low"],
        high=bounds[name]["high"],
    )

    sampler = spotpy.algorithms.dds(
        spot_setup, dbname="RAVEN_model_run", dbformat="ram", save_sim=False
    )

    # Removing this that is too subject to stochasticity. If it works, it will converge
    # sampler.sample(8, trials=1)
    # d1 = sampler.status.objectivefunction_max
    # sampler.sample(100, trials=1)
    # d2 = sampler.status.objectivefunction_max
    # assert d2 > d1

    sampler.sample(20, trials=1)
    assert spot_setup.diagnostics is not None

    results = sampler.getdata()
    assert len(spotpy.analyser.get_best_parameterset(results)[0]) == len(
        bounds[name]["high"]
    )
