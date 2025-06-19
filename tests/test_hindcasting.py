import datetime as dt
import pathlib
import sys

import pytest
import xarray as xr

from ravenpy import Emulator, EnsembleReader
from ravenpy.config import commands as rc
from ravenpy.config.emulators import GR4JCN
from ravenpy.utilities.forecasting import (
    hindcast_climatology_esp,
    to_climpred_hindcast_ensemble,
    warm_up,
)

"""
Test to perform a hindcast using Caspar data on THREDDS.
Currently only runs GEPS, eventually will run GEPS, GDPS, REPS and RDPS.
To do so will need to add the actual data from Caspar, currently being downloaded
but this is a good proof of concept.
"""


class TestHindcasting:
    def test_hindcasting_GEPS(self, salmon_hru, tmp_path, yangtze):
        ts20 = yangtze.fetch("caspar_eccc_hindcasts/geps_watershed.nc")

        hru = salmon_hru["land"]
        data_kwds = {
            "ALL": {"Latitude": hru["latitude"], "Longitude": hru["longitude"]},
            "PRECIP": {
                "Deaccumulate": True,
                "TimeShift": -0.25,
                "LinearTransform": {"scale": 1000},
            },
        }

        pre_conf = GR4JCN(
            StartDate=dt.datetime(2018, 6, 10),
            EndDate=dt.datetime(2018, 6, 14),
            HRUs=[
                hru,
            ],
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
            GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
            Gauge=[
                rc.Gauge.from_nc(
                    ts20, data_type=["PRECIP", "TEMP_AVE"], data_kwds=data_kwds
                ),
            ],
        )

        init_conf = warm_up(pre_conf, duration=8, workdir=tmp_path / "wup")

        # Create gauges for each member and run the model
        nm = 3
        ens = []
        for i in range(nm):
            g = rc.Gauge.from_nc(
                ts20,
                station_idx=i + 1,
                data_type=["PRECIP", "TEMP_AVE"],
                data_kwds=data_kwds,
            )
            assert g.data[0].read_from_netcdf.station_idx == i + 1
            for d in g.data:
                if d.data_type == "PRECIP":
                    assert d.read_from_netcdf.linear_transform.scale == 1000

            mem_conf = init_conf.model_copy(
                update=dict(
                    gauge=[g],
                ),
            )

            # This segfaults
            e = Emulator(
                config=mem_conf,
                workdir=tmp_path / f"m{i:02}",
            )
            ens.append(e.run())

        out = EnsembleReader(runs=ens)

        # The model now has the forecast data generated, and it has 5 days of forecasts.
        assert len(out.hydrograph.member) == 3
        assert len(out.hydrograph.time) == 5

    # Skip if using Python3.10
    @pytest.mark.skipif(
        (3, 11) > sys.version_info >= (3, 10),
        reason="climpred is unstable in Python 3.10",
    )
    def test_climpred_hindcast_verif(self, salmon_hru, tmp_path, yangtze):
        ts = pathlib.Path(
            yangtze.fetch(
                "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
            )
        )
        # Make a local copy to evade double-ownership of file - first file
        ts_tmp1 = tmp_path / "salmon_river_near_prince_george-tmp1.nc"
        ts_tmp1.write_bytes(ts.read_bytes())
        # Make a local copy to evade double-ownership of file - second file
        ts_tmp2 = tmp_path / "salmon_river_near_prince_george-tmp2.nc"
        ts_tmp2.write_bytes(ts.read_bytes())

        # This is the forecast start date, on which the forecasts will be launched.
        start_date = dt.datetime(1980, 6, 1)

        # Provide the length of the forecast, in days:
        forecast_duration = 100

        # Define HRU to build the hydrological model
        hru = salmon_hru["land"]

        # Set alternative names for netCDF variables
        alt_names = {
            "TEMP_MIN": "tmin",
            "TEMP_MAX": "tmax",
            "RAINFALL": "rain",
            "SNOWFALL": "snow",
        }

        # Data types to extract from netCDF
        data_type = ["TEMP_MAX", "TEMP_MIN", "RAINFALL", "SNOWFALL"]
        data_kwds = {
            "ALL": {
                "elevation": hru[
                    "elevation"
                ],  # No need for lat/lon as they are included in the netcdf file already
            }
        }

        # Model configuration
        model_config = GR4JCN(
            params=[0.529, -3.396, 407.29, 1.072, 16.9, 0.947],
            Gauge=[
                rc.Gauge.from_nc(
                    ts_tmp1,
                    data_type=data_type,
                    alt_names=alt_names,
                    data_kwds=data_kwds,
                )
            ],
            HRUs=[hru],
            StartDate=start_date,
            Duration=forecast_duration,
            GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
        )

        hindcasts = hindcast_climatology_esp(
            config=model_config,  # Note that the forecast duration is already set-up in the model_config above.
            warm_up_duration=365,  # number of days for the warm-up
            years=[1985, 1986, 1987, 1988, 1989, 1990],
            hindcast_years=[2001, 2002, 2003, 2004, 2005, 2006, 2007],
        )

        q_obs = xr.open_dataset(ts_tmp2)

        # However, our simulated streamflow is named "q_sim" and climpred requires the observation to be named the same thing
        # so let's rename it. While we're at it, we need to make sure that the identifier is the same. In our observation
        # dataset, it is called "nstations" but in our simulated streamflow it's called "nbasins". Here we standardize.
        q_obs = q_obs.rename({"qobs": "q_sim", "nstations": "nbasins"})

        # Make the hindcasting object we can use to compute statistics and metrics
        hindcast_object = to_climpred_hindcast_ensemble(hindcasts, q_obs)

        # This function is used to convert to binary to see if yes/no forecast is larger than observations
        def pos(x):
            return x > 0  # Check for binary outcome

        # Rank histogram verification metric
        rank_histo_verif = hindcast_object.verify(
            metric="rank_histogram",
            comparison="m2o",
            dim=["member", "init"],
            alignment="same_inits",
        )
        assert "q_sim" in rank_histo_verif
        assert rank_histo_verif.q_sim.shape[0] - 1 == forecast_duration

        # CRPS verification metric
        crps_verif = hindcast_object.verify(
            metric="crps",
            comparison="m2o",
            dim=["member", "init"],
            alignment="same_inits",
        )
        assert "q_sim" in crps_verif
        assert crps_verif.q_sim.shape[0] - 1 == forecast_duration

        # TODO: Fix this part
        """
        reliability_verif = hindcast_object.verify(
            metric="reliability",
            comparison="m2o",
            dim=["member", "init"],
            alignment="same_inits",
            logical=pos,
        )

        assert "flow" in reliability_verif
        assert reliability_verif.flow.shape[0] == forecast_duration
        """
