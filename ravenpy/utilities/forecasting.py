#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 09:11:58 2020

@author: ets
"""

import datetime as dt
import logging
import warnings
from typing import Dict, List, Tuple

import climpred
import numpy as np
import pandas as pd
import xarray as xr
from climpred import HindcastEnsemble

from ravenpy.models.emulators import get_model

LOGGER = logging.getLogger("PYWPS")


# TODO: Complete docstrings

# This function gets model states after running the model (i.e. states at the end of the run).
def get_raven_states(model: str, workdir=None, **kwargs) -> Dict:
    """Get the RAVEN states file (.rvc file) after a model run.

    Parameters
    ----------
    model : {'HMETS', 'GR4JCN', 'MOHYSE', 'HBVEC'}
      Model name.
    kwargs : dict
      Model configuration parameters, including the forcing files (ts).

    Returns
    -------
    dict
      Raven model forcing file

    """
    # Run the model and get the rvc file for future hot start.
    m = get_model(model)(workdir=workdir)
    m(overwrite=True, **kwargs)
    rvc = m.outputs["solution"]

    return rvc


# Do the actual forecasting step
def perform_forecasting_step(rvc, model, workdir=None, **kwds):
    """
    Function that might be useful eventually to do a forecast from a model setup.
    """
    # kwds includes 'ts', the forecast timeseries data
    # Setup the model
    m = get_model(model)(workdir=workdir)

    # Force the initial conditions
    m.resume(rvc)

    # Set the parameters, start dates, etc. required to run the model and run
    m(overwrite=True, **kwds)

    return m.q_sim


def perform_climatology_esp(
    model_name: str,
    forecast_date: dt.datetime,
    forecast_duration: int,
    workdir=None,
    **kwargs,
) -> xr.DataArray:
    """
    This function takes the model setup and name as well as forecast data and duration and returns
    an ESP forecast netcdf. The data comes from the climatology data and thus there is a mechanism
    to get the correct data from the time series and exclude the current year.

    Parameters
    ----------
    model_name : {'HMETS', 'MOHYSE', 'GR4JCN', 'HBVEC'}
      Model name to instantiate Raven model.
    forecast_date : datetime.datetime
      Date of the forecast issue.
    forecast_duration : int
      Number of days of forecast, forward looking.
    kwargs : dict
      Raven model configuration parameters.

    Returns
    -------
    xarray.DataArray
      Array of streamflow values from the ESP method along with list of member years

    """
    # Get the timeseries
    tsnc = xr.open_dataset(kwargs["ts"])

    # Prepare model instance
    m = get_model(model_name)(workdir=workdir)

    # Now find the periods of time for warm-up and forecast and add to the model keywords as the defaults are failing
    # (nanoseconds datetimes do not like the year 0001...)
    start_date = pd.to_datetime(tsnc["time"][0].values)
    start_date = start_date.to_pydatetime()

    kwargs["start_date"] = start_date

    # Forecasting from Feb 29th is not ideal, we will replace with Feb 28th.
    # Should not change much in a climatological forecast.
    if forecast_date.month == 2 and forecast_date.day == 29:
        forecast_date.replace(day=28)

    # Check to make sure forecast date is not in the first year as we need model warm-up.
    # We cannot use timedelta because if the dataset happens to start on a leap
    # year, then the timedelta=365 days will not be robust. (and we cannot use timedelta(years=1)...)
    dateLimit = start_date.replace(year=start_date.year + 1)
    if dateLimit > forecast_date:
        msg = (
            "Forecast date is within the warm-up period. Select another forecast date."
        )
        warnings.warn(msg)

    # initialize the array of forecast variables
    qsims = []

    # list of unique years in the dataset:
    avail_years = list(np.unique(tsnc["time.year"].data))

    # Take a copy of the forecast initial date before overwriting in the forecast step.
    forecast_date_main = forecast_date

    # Remove the year that we are forecasting. Or else it's cheating!
    avail_years.remove(forecast_date.year)

    # Update the forecast end-date, which will be the day prior to the forecast date.
    # So forecasts warm-up will be from day 1 in the dataset to the forecast date.
    kwargs["end_date"] = forecast_date - dt.timedelta(days=1)

    # Get RVC file if it exists, else compute it.
    if "rvc" in kwargs and len(kwargs["rvc"]) > 0:
        rvc = kwargs.pop("rvc")
    else:
        # Run model to get rvc file after warm-up using base meteo
        rvc = get_raven_states(model_name, workdir=workdir, **kwargs)

    # We need to check which years are long enough (ex: wrapping years, 365-day forecast starting in
    # September 2015 will need data up to August 2016 at least)
    for years in avail_years:
        if forecast_date.replace(year=years) + dt.timedelta(
            days=forecast_duration - 1
        ) > pd.to_datetime(tsnc["time"][-1].values):
            avail_years.remove(years)
            msg = (
                f"Year {years} has been removed because it is the last year in the dataset and does not cover the "
                f"forecast duration."
            )
            warnings.warn(msg)

    # We will iterate this for all forecast years
    for years in avail_years:

        # Replace the forecast period start and end dates with the climatological ESP dates for the
        # current member (year)
        forecast_date = forecast_date.replace(year=years)
        kwargs["start_date"] = forecast_date
        kwargs["end_date"] = forecast_date + dt.timedelta(days=forecast_duration - 1)

        # Setup the initial states from the warm-up and run the model.
        # Note that info on start/end dates and timeseries are in the kwds.
        m.resume(rvc)
        m(run_name=f"run_{years}", **kwargs)

        # Add member to the ensemble and retag the dates to the real forecast dates
        # (or else we will get dates from the climate dataset that cover all years)
        new_member = m.q_sim.copy(deep=True)
        new_member["time"] = pd.date_range(
            forecast_date_main, periods=forecast_duration
        )
        qsims.append(new_member)

    # Concatenate the members through a new dimension for the members and remove unused dims.
    qsims = xr.concat(qsims, dim="member")
    qsims = qsims.squeeze()

    # Add the number of the forecast year as member ID
    qsims["member"] = (["member"], avail_years)

    return qsims


def make_climpred_hindcast_object(
    hindcast: xr.Dataset, observations: xr.Dataset
) -> climpred.HindcastEnsemble:
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
    climpred.HindcastEnsemble
      The hindcast ensemble formatted to be used in climpred.

    """

    # Todo: Add verification that the variable names are the same
    # Todo: Add verification of sizes, catch and return message.

    # Make the hindcastEnsemble object for the hindcast data
    hindcast_obj = HindcastEnsemble(hindcast)
    # Add the observations to that hindcastEnsemble object for verification.
    hindcast_obj = hindcast_obj.add_observations(observations)

    return hindcast_obj


def make_ESP_hindcast_dataset(
    model_name: str,
    forecast_date: dt.datetime,
    included_years: List[int],
    forecast_duration: int,
    **kwargs,
) -> Tuple[xr.Dataset, xr.Dataset]:
    """Make a hindcast using ESP dataset.

    This function takes the required information to run a RAVEN model run in ESP
    mode, and runs RAVEN for a series of queried initialization dates. For each
    requested date, it will perform an ESP forecast using all climate data at
    its disposal. Returns the hindcast and observations datasets required by climpred.

    Parameters
    ----------
    model_name: str
      Name of the RAVEN model setup (GR4JCN, MOHYSE, HMETS or HBVEC).
    forecast_date: datetime.datetime
      Calendar date that is used as the base to perform hindcasts.
      The year component is updated for each initialization date.
    included_years: List[int]
      the years that we want to perform a hindcasting for, on the calendar date in "foreacast_date"
    forecast_duration: int
      Duration of the forecast, in days. Refers to the longest lead-time required.
    **kwargs
      All other parameters needed to run RAVEN.

    Returns
    -------
    qsims: xarray.Dataset
      The dataset containing the (init, member, lead) dimensions ready for using in climpred. (qsim)
    qobs: xarray.Dataset
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
    qsims = perform_climatology_esp(
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

        qsims_tmp = perform_climatology_esp(
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
