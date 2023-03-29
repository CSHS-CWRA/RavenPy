#!/usr/bin/env python3
"""
Created on Fri Jul 17 09:11:58 2020

@author: ets
"""
import collections
import datetime as dt
import logging
from pathlib import Path
from typing import List

import pandas as pd
import xarray as xr
from climpred import HindcastEnsemble

from ravenpy import Emulator, EnsembleReader
from ravenpy.new_config.rvs import Config

LOGGER = logging.getLogger("PYWPS")


def climatology_esp(
    config, path, years: List[int] = None, overwrite=False
) -> EnsembleReader:
    """
    Ensemble Streamflow Prediction based on historical variability.

    Run the model using forcing for different years.
    No model warm-up is performed by this function, make sure the initial states are
    consistent with the start date.

    Parameters
    ----------
    config : Config
      Model configuration.
    path: str, Path
      Path to rv files and model outputs.
    years : List[int]
      Years from which forcing time series will be drawn. If None, run for all years where
      forcing data is available.

    Returns
    -------
    EnsembleReader
      Class facilitating the analysis of multiple Raven outputs.
    """
    path = Path(path)

    # Time from meteo forcing files.
    time = _get_time(config)

    if years is None:
        # Create list of years
        years = list(range(time[0].dt.year.values, time[-1].dt.year.values + 1))

    # Define duration and end_date
    # Make sure configuration uses duration instead of end_date
    if config.duration is None:
        if config.end_date is None or config.start_date is None:
            raise ValueError("StartDate and Duration must be configured.")
        duration = config.end_date - config.start_date
        end = config.end_date
        config.duration = duration.days
        config.end_date = None

    else:
        if config.start_date is None:
            raise ValueError("StartDate must be configured.")
        duration = config.duration
        end = config.start_date + dt.timedelta(days=duration)

    # Define start date
    start = config.start_date
    if start.month == 2 and start.day == 29:
        start = start.replace(day=28)

    # Run the model for each year
    ensemble = []
    for year in years:
        # Prepare model instance
        config.start_date = s = start.replace(year=year)
        if len(time.sel(time=slice(str(s), str(end.replace(year=year))))) < duration:
            continue

        out = Emulator(config=config, workdir=path / f"Y{year}").run(
            overwrite=overwrite
        )

        # Format output so that `time` is set to the forecast time, not the forcing time.
        if "hydrograph" in out.files:
            out.files["hydrograph"] = _shift_esp_time(
                out.files["hydrograph"], start.year
            )
        if "storage" in out.files:
            out.files["storage"] = _shift_esp_time(out.files["storage"], start.year)
        ensemble.append(out)

    return EnsembleReader(config.run_name, runs=ensemble)


def to_climpred_hindcast_ensemble(hindcast, observations):
    """
    Ceates a hindcasting object that can be used by the `climpred`
     toolbox for hindcast verification.

    Parameters
    ----------
    hindcast : xarray.Dataset
      The hindcasted streamflow data for a given period.
    observations : xarray.Dataset
      The streamflow observations that are used to verify the hindcasts.

    Returns
    -------
    climpred.HindcastEnsemble object
      The hindcast ensemble formatted to be used in climpred.

    """

    # Make sure the variable names are the same.
    var_names = set(hindcast.data_vars().keys()).intersection(
        observations.data_vars().keys()
    )
    if len(var_names) == 0:
        raise ValueError(
            "`hindcast` and `observations` should have a variable with the same name."
        )

    # Todo: Add verification of sizes, catch and return message.

    # Make the hindcastEnsemble object for the hindcast data
    hindcast_obj = HindcastEnsemble(hindcast)

    # Add the observations to that hindcastEnsemble object for verification.
    hindcast_obj = hindcast_obj.add_observations(observations)

    return hindcast_obj


def warm_up(config, duration: int, path: Path, overwrite=False):
    """Run the model on a time series preceding the start date.

    Parameters
    ----------
    config : Config
      Model configuration.
    duration : int
      Number of days the warm-up simulation should last *before* the start date.
    path: Path
      Work directory.
    overwrite: bool
      If True, overwrite existing files.

    Returns
    -------
    Config
      Model configuration with initial state set by running the model prior to the start date.
    """
    wup = config.copy()
    wup.start_date = config.start_date - dt.timedelta(days=duration)
    wup.duration = duration
    wup.end_date = None

    out = Emulator(wup, workdir=path).run(overwrite=overwrite)
    return config.set_solution(out.files["solution"])


def hindcast_climatology_esp(
    config,
    path,
    warm_up_duration,
    years=None,
    hindcast_years=None,
    overwrite: bool = False,
) -> xr.DataArray:
    """Hindcast of Ensemble Prediction Streamflow.

    This function runs an emulator initialized for each year in `hindcast_years`,
    using the forcing time series for each year in `years`. This allows an assessment
    of the performance of the ESP. The total number of simulations is given by
    `len(years) * len(hindcasts_years)`.

    Parameters
    ----------
    config: Config
      Model configuration. Initial states will be overwritten.
    path: Path
      Work directory.
    warm_up_duration : int
      Number of days to run the model prior to the starting date to initialize the state variables.
    years : List[int]
      Years from which forcing time series will be drawn. If None, run for all years where
      forcing data is available.
    hindcast_years:  List[int]
      Years for which the model will be initialized and the `climatology_esp` function run.
      Defaults to all years when forcing data is available.
    overwrite: bool
      If True, overwrite existing files.

    Returns
    -------
    xarray.DataArray
      The array containing the (init, member, lead) dimensions ready for using in climpred. (qsim)

    Notes
    -----
    The dataset output dimensions are
     - `init`: hindcast issue date,
     - `member`: ESP members of the hindcasting experiment,
     - `lead`: number of lead days of the forecast.
    """
    path = Path(path)

    # `hindcast_years` defaults to all available years
    time = _get_time(config)
    if hindcast_years is None:
        # Create list of years
        hindcast_years = list(
            range(time[0].dt.year.values, time[-1].dt.year.values + 1)
        )

    # Define start date
    start = config.start_date
    if start.month == 2 and start.day == 29:
        start = start.replace(day=28)

    # Run climatology_esp for each year in `hindcast_years`, with initial conditions set
    # by a warm-up simulation.
    q_sims = {}
    for year in hindcast_years:
        # Compute initial state for each year
        config.start_date = e = start.replace(year=year)
        s = e - dt.timedelta(days=warm_up_duration)
        if len(time.sel(time=slice(str(s), str(e)))) < warm_up_duration:
            continue

        # Run warm-up simulation to set initial state
        conf = warm_up(
            config,
            duration=warm_up_duration,
            path=path / "wup" / f"Y{year}",
            overwrite=overwrite,
        )

        # Run climatology ESP
        esp = climatology_esp(
            conf, years=years, path=path / "esp" / f"Y{year}", overwrite=overwrite
        )

        # Format results for climpred
        q_sim = esp.hydrograph.q_sim.expand_dims(
            init=[
                config.start_date,
            ]
        )
        q_sim = q_sim.rename({"time": "lead"})
        q_sim.lead.attrs["units"] = "days"
        q_sim = q_sim.assign_coords(lead=list(range(1, len(q_sim.lead) + 1)))

        q_sims[year] = q_sim

    # Format results for climpred
    q_sims = xr.concat(q_sims.values(), dim="init")
    q_sims.lead.attrs["units"] = "days"
    # q_sims = q_sims.assign_coords(init=q_sims.init)
    return q_sims


def _shift_esp_time(nc, year, dim="member"):
    """Modify netCDF time to facilitate ensemble analysis out of ESP forecasts.

    Modify time such that it starts on the given year, and
    add member dimension with original year as value.
    """

    ds = xr.open_dataset(nc, use_cftime=True)

    # Create new time coordinate
    start = ds.time.data[0]
    freq = xr.infer_freq(ds.time)
    ds["time"] = xr.cftime_range(
        start.replace(year=year), periods=len(ds.time), freq=freq
    )

    # New coordinate dimension to store the original year
    out = ds.expand_dims(
        {
            dim: [
                start.year,
            ]
        }
    )

    # Write to disk
    fn = nc.with_stem(nc.stem + "_shifted")
    out.to_netcdf(fn, mode="w")
    return fn


def _get_time(config):
    """Return time DataArray."""
    if config.gauge is not None:
        time = config.gauge[0].ds.time
    elif config.station_forcing is not None:
        time = config.station_forcing[0].da.time
    elif config.gridded_forcing is not None:
        time = config.gridded_forcing[0].da.time
    else:
        raise NotImplementedError
    return time