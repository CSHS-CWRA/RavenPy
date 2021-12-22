import datetime as dt

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr
import xskillscore as xss

from ravenpy.models import GR4JCN
from ravenpy.utilities.data_assimilation import (
    assimilation_initialization,
    perturb_full_series,
    perturbation,
    sequential_assimilation,
)
from ravenpy.utilities.testdata import get_local_testdata


def test_perturbation():
    ts = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    ds = xr.open_dataset(ts)

    tmax = ds.tmax.isel(time=slice(0, 10))
    p_tmax = perturbation(tmax, "norm", 0.01, members=50)
    np.testing.assert_allclose(p_tmax.mean("members"), tmax, rtol=0.1)

    rain = ds.rain.isel(time=slice(30, 60))
    p_rain = perturbation(rain, "gamma", 0.01, members=50)
    np.testing.assert_allclose(p_rain.mean("members"), rain, rtol=0.1)

    assert p_tmax.attrs == ds.tmax.attrs
    assert p_rain.attrs == ds.rain.attrs


class TestAssimilationGR4JCN:
    def test_simple(self):

        # get timeseries
        ts = get_local_testdata(
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
        )

        # set number of members. Using 7 here to make it easier to find and debug.
        n_members = 7

        # Perturbation parameters for the assimilation, keyed by standard_name
        std = {
            "rainfall": 0.30,
            "prsn": 0.30,
            "tasmin": 2.0,
            "tasmax": 2.0,
            "water_volume_transport_in_river_channel": 0.10,
        }

        # Perturbation distribution
        dists = {
            "pr": "gamma",
            "rainfall": "gamma",
            "prsn": "gamma",
            "water_volume_transport_in_river_channel": "rnorm",
        }

        qkey = "water_volume_transport_in_river_channel"
        if qkey not in std:
            raise ValueError("Assimilation requires perturbing the flow variable.")

        # Assimilation variables (from HRUStateVariable)
        assim_var = ("soil0", "soil1")

        # Assimilation period (days between each assimilation step)
        assim_step_days = 3

        # GR4JCN model instance
        model = GR4JCN()

        # set the start and end dates for the first assimilation period, warm-up
        start_date = dt.datetime(1996, 9, 1)
        end_date = dt.datetime(1996, 9, 30)

        # Catchment properties to populate model
        area = 4250.6
        elevation = 843.0
        latitude = 54.4848
        longitude = -123.3659
        params = (0.1353389, -0.005067198, 576.8007, 6.986121, 1.102917, 0.9224778)

        # Do the first assimilation pass to get hru_states and basin_states.
        # Can be skipped if there is already this data from a previous run.
        model, xa, hru_states, basin_states = assimilation_initialization(
            model,
            ts,
            start_date=start_date,
            end_date=start_date + dt.timedelta(days=assim_step_days - 1),
            area=area,
            elevation=elevation,
            latitude=latitude,
            longitude=longitude,
            params=params,
            assim_var=assim_var,
            n_members=n_members,
        )

        # Perturb the inputs for the rest of the assimilation
        perturbed = perturb_full_series(
            model,
            std=std,
            start_date=start_date,
            end_date=end_date,
            dists=dists,
            n_members=n_members,
        )

        # Get observed streamflow for computing results later
        q_obs = xr.open_dataset(ts)["qobs"].sel(time=slice(start_date, end_date))

        # Create netcdf for the model.
        p_fn = model.workdir / "perturbed_forcing.nc"
        perturbed = xr.Dataset(perturbed)
        perturbed.to_netcdf(p_fn, mode="w")

        # Run the sequential assimilation for the entire period.
        q_assim, hru_states, basin_states = sequential_assimilation(
            model,
            hru_states,
            basin_states,
            p_fn,
            q_obs,
            assim_var,
            start_date=start_date + dt.timedelta(days=assim_step_days),
            end_date=end_date,
            n_members=n_members,
            assim_step_days=assim_step_days,
        )

        # ==== Reference run ====
        model.config.rvi.run_name = "ref"
        model.config.rvi.start_date = start_date
        model.config.rvi.end_date = end_date

        model.config.rvc.hru_states = {}
        model.config.rvc.basin_states = {}
        model.config.rvc.soil0 = None
        model.config.rvc.soil1 = 15

        model([ts])

        # We can now plot everything!
        plt.plot(q_assim.T, "r", label="Assimilated")  # plot the assimilated flows
        plt.plot(q_obs.T, "b", label="Observed")  # plot the observed flows
        plt.plot(
            model.q_sim, "g", label="Simulated"
        )  # plot the open_loop (simulation with no assimilation)
        # plt.legend()
        # plt.show()

        # print('RMSE - Assimilated: ' + str(xss.rmse(q_assim.mean(dim='state').T,q_obs[0:q_assim.shape[1]].T).data))
        # print('RMSE - Open-Loop: ' + str(xss.rmse(model.q_sim[0:q_assim.shape[1],0],q_obs[0:q_assim.shape[1]].T).data))

        assert q_assim.shape[0] == n_members
