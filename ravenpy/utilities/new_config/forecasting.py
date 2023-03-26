#!/usr/bin/env python3
"""
Created on Fri Jul 17 09:11:58 2020

@author: ets
"""
import collections
import datetime as dt
import logging
import re
import warnings
from pathlib import Path
from typing import List, Sequence, Tuple

import cftime
import numpy as np
import pandas as pd
import xarray as xr
from climpred import HindcastEnsemble

from ravenpy import Emulator, EnsembleReader
from ravenpy.new_config.rvs import Config

from .. import gis_import_error_message

try:
    import rioxarray
except (ImportError, ModuleNotFoundError) as e:
    msg = gis_import_error_message.format(Path(__file__).stem)
    raise ImportError(msg) from e

from ravenpy.models.emulators import get_model

LOGGER = logging.getLogger("PYWPS")


def climatology_esp(config, path, years: List[int] = None) -> List[Emulator]:
    """
    Ensemble Streamflow Prediction based on historical variability.

    Run the model using forcing for different years.
    The initial states should be consistent with the start date.


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
    """
    path = Path(path)

    # Time from meteo forcing files.
    time = config.gauges.ds.time

    if years is None:
        # Create list of years
        years = list(range(time[0].dt.year.values, time[-1].dt.year.values + 1))

    # Define duration and end_date
    # Make sure configuration uses duration instead of end_date
    if config.duration is None:
        duration = config.end_date - config.start_date
        end = config.end_date
        config.duration = duration.days
        config.end_date = None

    else:
        duration = config.duration
        end = config.start_date + dt.timedelta(days=duration)

    # Define start date
    start = config.start_date
    if start.month == 2 and start.day == 29:
        start = start.replace(day=28)

    # Select only years with forcing data
    valid_years = []
    for year in years:
        s = start.replace(year=year)
        e = end.replace(year=year)
        if len(time.sel(slice(s, e))) >= duration:
            # This assumes time series is daily or higher resolution - not fool proof
            valid_years.append(year)

    ensemble = []
    for year in valid_years:
        # Prepare model instance
        config.start_date = config.start_date.replace(year=year)
        e = Emulator(config=config, path=path / f"Y{year}")
        e.run()
        if "hydrograph" in e.output_files:
            _shift_esp_time(e.output_files["hydrograph"])
        if "storage" in e.output_files:
            _shift_esp_time(e.output_files["storage"])
        ensemble.append(e)

    return EnsembleReader(config.run_name, paths=[e.path for e in ensemble])


def make_climpred_hindcast_object(hindcast, observations):
    """
    This function takes a hindcasting dataset of streamflow as well as associated
    observations and creates a hindcasting object that can be used by the
    climpred toolbox for hindcast verification.

    Parameters
    ----------
    hindcast : xarray.Dataset
      The hindcasting streamflow data for a given period
    observations : xarray.Dataset
      The streamflow observations that are used to verify the hindcasts

    Returns
    -------
    hindcast_obj : climpred.HindcastEnsemble object
      The hindcast ensemble formatted to be used in climpred.

    """

    # Todo: Add verification that the variable names are the same
    # Todo: Add verification of sizes, catch and return message.

    # Make the hindcastEnsemble object for the hindcast data
    hindcast_obj = HindcastEnsemble(hindcast)
    # Add the observations to that hindcastEnsemble object for verification.
    hindcast_obj = hindcast_obj.add_observations(observations)

    return hindcast_obj


# TODO: convert to new config
def make_ESP_hindcast_dataset(
    model_name, forecast_date, included_years, forecast_duration, **kwargs
) -> Tuple[xr.Dataset, xr.Dataset]:
    """Make a hindcast using ESP dataset.

    This function takes the required information to run a RAVEN model run in ESP
    mode, and runs RAVEN for a series of queried initialization dates. For each
    requested date, it will perform an ESP forecast using all climate data at
    its disposal. Returns the hindcast and observations datasets required by climpred.

    Parameters
    ----------
    model_name: string
      Name of the RAVEN model setup (GR4JCN, MOHYSE, HMETS or HBVEC).
    forecast_date: datetime.datetime
      Calendar date that is used as the base to perform hindcasts.
      The year component is updated for each initialization date.
    included_years: List[int]
      the years that we want to perform a hindcasting for, on the calendar date in "foreacast_date"
    forecast_duration: int
      Duration of the forecast, in days. Refers to the longest lead-time required.
    kwargs: dict
      All other parameters needed to run RAVEN.

    Returns
    -------
    xarray.Dataset
      The dataset containing the (init, member, lead) dimensions ready for using in climpred. (qsim)
    xarray.Dataset
      The dataset containing all observations over the verification period.
      The only dimension is 'time', as required by climpred. No other processing
      is required as climpred will find corresponding verification dates on its own. (qobs)

    Notes
    -----
    init = hindcast issue date
    member = ESP members of the hindcasting experiment
    lead = number of lead days of the forecast.
    """

    # Use the ESP functions to generate the ESP hindcasts for the first period (first initialization)
    qsims = climatology_esp(
        model_name,
        forecast_date.replace(year=included_years[0]),
        forecast_duration,
        **kwargs,
    )

    # climpred needs the 'time' variable name to be 'lead' for the hindcasts.
    # Also add the coordinates for the lead and member dimensions.
    qsims = qsims.rename({"time": "lead"})
    qsims = qsims.assign_coords(lead=list(range(1, forecast_duration + 1)))
    qsims = qsims.assign_coords(member=list(range(1, qsims.data.shape[0] + 1)))

    # Repeat the process for all hindcast years required. Could be parallelized by a pro!
    for i in included_years[1:]:
        qsims_tmp = climatology_esp(
            model_name, forecast_date.replace(year=i), forecast_duration, **kwargs
        )
        qsims_tmp = qsims_tmp.rename({"time": "lead"})
        qsims_tmp = qsims_tmp.assign_coords(lead=list(range(1, forecast_duration + 1)))
        qsims_tmp = qsims_tmp.assign_coords(
            member=list(range(1, qsims_tmp.data.shape[0] + 1))
        )

        # Concatenate in the 'init' dimension.
        qsims = xr.concat([qsims, qsims_tmp], dim="init")

    qsims["lead"] = list(range(1, forecast_duration + 1))

    # Make the list of hindcast dates for populating the 'init' dimension coordinates.
    date_list = [
        forecast_date.replace(year=included_years[x])
        for x in range(len(included_years))
    ]

    # Other processing required by climpred.
    qsims = qsims.assign_coords(init=pd.to_datetime(date_list))
    qsims = qsims.to_dataset(name="flow")

    # Here units are days and always will be!
    qsims["lead"].attrs["units"] = "days"

    # Extract and prepare the observations for verification.
    qobs = xr.open_dataset(kwargs["ts"]).qobs

    # The hindcast and observation variable names MUST be the same. Here we
    # use "flow" but it could be abything else, really.
    qobs = qobs.to_dataset(name="flow")

    return qsims, qobs


def _shift_esp_time(nc, year, dim="member"):
    """Modify netCDF time to facilitate ensemble analysis out of ESP forecasts.

    Modify time such that it starts on the given year, and
    add member dimension with original year as value.
    """

    ds = xr.open_dataset(nc, mode="a", use_cftime=True)

    # Create new time coordinate
    start = ds.time.data[0]
    freq = xr.infer_freq(ds.time)
    ds["time"] = xr.date_range(
        start.replace(year=year), periods=len(ds.time), freq=freq
    )

    # New coordinate dimension to store the original year
    out = ds.expand_dims({dim: start.year})

    # Write to disk
    out.to_netcdf(nc)


def _shift_esp_time_v1(nc, year, dim="member"):
    """Modify netCDF time to facilite ensemble analysis out of ESP forecasts.

    Modify time such that it starts on the given year, and
    add member dimension with original year as value.
    """

    ds = xr.open_dataset(nc, mode="a")
    calendar = ds.time.encoding["calendar"]
    orig_year = int(ds.time[0].dt.year)

    # Original start and time delta
    # Simpler but only valid for daily:
    # new_member["time"] = pd.date_range(forecast_date_main, periods=forecast_duration)
    time = ds.time.values.astype("datetime64[s]")
    start = time[0]
    delta = time - time[0]

    # New start and time series
    tt = cftime._parse_date(str(start))
    new_start = cftime.datetime(*tt, calendar=calendar).replace(year=year)
    new_time = cftime.num2date(
        delta, units=f"seconds since {new_start}", calendar=calendar
    )

    # New coordinate dimension to store the original year
    ds["time"] = new_time
    out = ds.expand_dims({dim: orig_year})

    # Write to disk
    out.to_netcdf(nc)
