import datetime as dt

import numpy as np
import pytest
import xarray as xr
from climpred import HindcastEnsemble

from ravenpy import Emulator, OutputReader
from ravenpy.utilities.forecasting import (
    climatology_esp,
    ensemble_prediction,
    hindcast_climatology_esp,
    warm_up,
)


def test_warm_up(minimal_emulator, tmp_path):
    conf = warm_up(minimal_emulator, duration=10, workdir=tmp_path)
    wup = OutputReader(path=tmp_path / "output")

    assert conf.start_date == minimal_emulator.start_date
    assert len(wup.hydrograph.q_sim.time) == 11
    assert conf.hru_state_variable_table != minimal_emulator.hru_state_variable_table


def test_climatology_esp(minimal_emulator, tmp_path):
    config = minimal_emulator.model_copy(deep=True)
    esp = climatology_esp(config, workdir=tmp_path, years=[1955, 1956])
    np.testing.assert_array_equal(esp.storage.member, [1955, 1956])
    assert len(esp.hydrograph.time) == minimal_emulator.duration + 1


def test_hindcast_climatology_esp(minimal_emulator, tmp_path, yangtze):
    config = minimal_emulator.model_copy(deep=True)
    hc = hindcast_climatology_esp(
        config,
        warm_up_duration=5,
        years=[1955, 1956, 1957],
        hindcast_years=[1990, 1991],
        workdir=tmp_path,
    )
    assert hc.sizes == {
        "lead": config.duration + 1,
        "init": 2,
        "member": 3,
        "nbasins": 1,
    }

    # Construct climpred HindcastEnsemble
    salmon_file = yangtze.fetch(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    qobs = xr.open_dataset(salmon_file).qobs

    qobs.name = "q_sim"
    hce = HindcastEnsemble(hc).add_observations(qobs)

    rank_histo_verif = hce.verify(
        metric="rank_histogram",
        comparison="m2o",
        dim=["member", "init"],
        alignment="same_inits",
    )
    assert rank_histo_verif.q_sim.shape[0] == config.duration + 1


def test_forecasting_GEPS(numeric_config, yangtze):
    """Test to perform a forecast using auto-queried ECCC data aggregated on THREDDS."""
    name, wup = numeric_config
    if name != "GR4JCN":
        pytest.skip("Test only for GR4JCN model.")

    geps = yangtze.fetch("eccc_forecasts/geps_watershed.nc")

    # Prepare a RAVEN model run using historical data, GR4JCN in this case.
    # This is a dummy run to get initial states. In a real forecast situation,
    # this run would end on the day before the forecast, but process is the same.
    wup.start_date = dt.datetime(2000, 1, 1)
    wup.duration = 30

    e = Emulator(wup)
    e.run()

    # Extract the final states that will be used as the next initial states
    conf = e.resume()

    # Set run parameters
    conf.start_date = "2020-12-10"
    conf.end_date = None
    conf.duration = 9

    # Collect test forecast data for location and climate model (20 members)
    nm = 20
    data_kwds = {
        "ALL": {
            "TimeShift": -0.25,
            "Latitude": conf.hrus[0].latitude,
            "Longitude": conf.hrus[0].longitude,
        },
        "PRECIP": {"Deaccumulate": True},
    }

    # Data types to extract from netCDF
    data_type = ["TEMP_AVE", "PRECIP"]

    # Run the 20 members
    out = ensemble_prediction(
        conf,
        geps,
        ens_dim="member",
        data_kwds=data_kwds,
        data_type=data_type,
    )

    # The model now has the forecast data generated and has 10 days of forecasts.
    assert len(out.hydrograph.time) == 10

    # Also see if GEPS has 20 members produced.
    assert len(out.hydrograph.q_sim.member) == nm

    # Check all members are different (checking snow because data in winter)
    assert len(set(out.storage.Snow.isel(time=-1).values)) == nm
