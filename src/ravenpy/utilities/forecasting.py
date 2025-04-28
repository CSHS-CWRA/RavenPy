#!/usr/bin/env python3
# Created on Fri Jul 17 09:11:58 2020
# @author: ets
import datetime as dt
import logging
import os
import tempfile
import warnings
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union
from urllib.parse import urlparse

import climpred
import xarray as xr
from climpred import HindcastEnsemble

from ravenpy import Emulator, EnsembleReader
from ravenpy.config import commands as rc
from ravenpy.config.rvs import Config

LOGGER = logging.getLogger("PYWPS")

# Can be set at runtime with `$ env RAVENPY_THREDDS_URL=https://xx.yy.zz/thredds/ ...`.
THREDDS_URL = os.environ.get(
    "RAVENPY_THREDDS_URL", "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/"
)
if not THREDDS_URL.endswith("/"):
    THREDDS_URL = f"{THREDDS_URL}/"


def climatology_esp(
    config,
    years: Optional[list[int]] = None,
    workdir: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
) -> EnsembleReader:
    """Ensemble Streamflow Prediction based on historical variability.

    Run the model using forcing for different years.
    No model warm-up is performed by this function, make sure the initial states are
    consistent with the start date.

    Parameters
    ----------
    config : ravenpy.config.rvs.Config
        Model configuration.
    years : List[int], optional
        Years from which forcing time series will be drawn.
        If None, run for all years where forcing data is available.
    workdir : str or Path
        The path to rv files and model outputs. If None, create a temporary directory.
    overwrite : bool
        Whether to overwrite existing values or not. Default: False.

    Returns
    -------
    EnsembleReader
        Class facilitating the analysis of multiple Raven outputs.
    """
    workdir = Path(workdir or tempfile.mkdtemp())

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

    # Get the start date that is immutable
    startdate_fixed = deepcopy(start)

    # Run the model for each year
    ensemble = []
    for year in years:
        # Prepare model instance
        config.start_date = s = startdate_fixed.replace(year=year)
        if len(time.sel(time=slice(str(s), str(end.replace(year=year))))) < duration:
            continue

        out = Emulator(config=config, workdir=workdir / f"Y{year}").run(
            overwrite=overwrite
        )

        # Format output so that `time` is set to the forecast time, not the forcing time.
        if "hydrograph" in out.files:
            out.files["hydrograph"] = _shift_esp_time(
                out.files["hydrograph"], startdate_fixed.year
            )
        if "storage" in out.files:
            out.files["storage"] = _shift_esp_time(
                out.files["storage"], startdate_fixed.year
            )
        ensemble.append(out)

    return EnsembleReader(run_name=config.run_name, runs=ensemble)


def to_climpred_hindcast_ensemble(
    hindcast: xr.Dataset, observations: xr.Dataset
) -> climpred.HindcastEnsemble:
    """
    Create a hindcasting object that can be used by the `climpred` toolbox for hindcast verification.

    Parameters
    ----------
    hindcast : xarray.Dataset
        The hindcasted streamflow data for a given period.
    observations : xarray.Dataset
        The streamflow observations that are used to verify the hindcasts.

    Returns
    -------
    climpred.HindcastEnsemble
        The hindcast ensemble formatted to be used in climpred.
    """
    # Make sure the variable names are the same.
    var_names = set(hindcast.data_vars.keys()).intersection(
        observations.data_vars.keys()
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


def warm_up(
    config,
    duration: int,
    workdir: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
) -> Config:
    """Run the model on a time series preceding the start date.

    Parameters
    ----------
    config : ravenpy.config.rvs.Config
        Model configuration.
    duration : int
        Number of days the warm-up simulation should last *before* the start date.
    workdir : Path
        Work directory.
    overwrite : bool
        If True, overwrite existing files.

    Returns
    -------
    ravenpy.config.rvs.Config
        Model configuration with initial state set by running the model prior to the start date.
    """
    wup = config.model_copy(deep=True)
    wup.start_date = config.start_date - dt.timedelta(days=duration)
    wup.duration = duration
    wup.end_date = None

    out = Emulator(wup, workdir=workdir).run(overwrite=overwrite)
    return config.set_solution(out.files["solution"])


def hindcast_climatology_esp(
    config: Config,
    warm_up_duration: int,
    years: Optional[list[int]] = None,
    hindcast_years: Optional[list[int]] = None,
    workdir: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
) -> xr.Dataset:
    """Hindcast of Ensemble Prediction Streamflow.

    This function runs an emulator initialized for each year in `hindcast_years`,
    using the forcing time series for each year in `years`. This allows an assessment
    of the performance of the ESP. The total number of simulations is given by
    `len(years) * len(hindcasts_years)`.

    Parameters
    ----------
    config : ravenpy.config.rvs.Config
        Model configuration. Initial states will be overwritten.
    warm_up_duration : int
        Number of days to run the model prior to the starting date to initialize the state variables.
    years : List[int]
        Years from which forcing time series will be drawn.
        If None, run for all years where forcing data is available.
    hindcast_years : List[int]
        Years for which the model will be initialized and the `climatology_esp` function run.
        Defaults to all years when forcing data is available.
    workdir : Path
        Work directory. If None, creates a temporary directory.
    overwrite : bool
        If True, overwrite existing files.

    Returns
    -------
    xarray.DataArray
        The array containing the (init, member, lead) dimensions ready for using in climpred (qsim).

    Notes
    -----
    The dataset output dimensions are:
     - `init`: hindcast issue date,
     - `member`: ESP members of the hindcasting experiment,
     - `lead`: number of lead days of the forecast.
    """
    workdir = Path(workdir or tempfile.mkdtemp())

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
            workdir=workdir / "wup" / f"Y{year}",
            overwrite=overwrite,
        )

        # Run climatology ESP
        esp = climatology_esp(
            conf, years=years, workdir=workdir / "esp" / f"Y{year}", overwrite=overwrite
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
    q_sims = q_sims.to_dataset()

    return q_sims


def _shift_esp_time(nc, year, dim="member"):
    """Modify netCDF time to facilitate ensemble analysis out of ESP forecasts.

    Modify time such that it starts on the given year, and add member dimension with original year as value.
    """
    ds = xr.open_dataset(nc)

    # Create new time coordinate
    start = ds.indexes["time"][0]
    freq = xr.infer_freq(ds.time)
    ds["time"] = xr.date_range(
        start.replace(year=year),
        periods=len(ds.time),
        freq=freq,
        calendar=ds.time.dt.calendar,
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
    fn = nc.parent.joinpath(f"{nc.stem}_shifted{nc.suffix}")
    out.to_netcdf(fn, mode="w")
    return fn


def _get_time(config: Config) -> xr.DataArray:
    """Return time DataArray."""
    if config.gauge is not None:
        time = config.gauge[0].ds.time
    elif config.station_forcing is not None:
        time = config.station_forcing[0].da.time
    elif config.gridded_forcing is not None:
        time = config.gridded_forcing[0].da.time
    else:
        raise NotImplementedError()
    return time


def ensemble_prediction(
    config,
    forecast: Union[str, Path],
    ens_dim: str = "member",
    workdir=None,
    overwrite=True,
    **kwds,
) -> EnsembleReader:
    r"""Ensemble Streamflow Prediction based on historical weather forecasts (CASPAR or other).

    Run the model using forcing for different years.
    No model warm-up is performed by this function, make sure the initial states are consistent with the start date.

    Parameters
    ----------
    config : ravenpy.config.rvs.Config
        Model configuration.
    forecast : str or Path
        Forecast subsetted to the catchment location (.nc).
    ens_dim : str
        Name of dimension to iterate over.
    workdir : str or Path
        The path to rv files and model outputs. If None, create temporary directory.
    overwrite : bool
        Overwrite files when writing to disk.
    \*\*kwds : dict
        Keywords for the `Gauge.from_nc` function.

    Returns
    -------
    EnsembleReader
        Class facilitating the analysis of multiple Raven outputs.
    """
    workdir = Path(workdir or tempfile.mkdtemp())

    # Run the model for each year
    ensemble = []

    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    forecast_ds = xr.open_dataset(forecast, decode_times=time_coder)

    for member in range(0, len(forecast_ds[ens_dim])):
        # Prepare model instance

        out = Emulator(
            config=config.model_copy(
                update={
                    "gauge": [
                        rc.Gauge.from_nc(
                            forecast,
                            station_idx=member + 1,
                            **kwds,
                        ),
                    ]
                },
            ),
            workdir=workdir / f"Y{member}",
        ).run(overwrite=overwrite)

        # Append to the ensemble.
        ensemble.append(out)

    return EnsembleReader(run_name=config.run_name, runs=ensemble, dim=ens_dim)


# Alias
hindcast_from_meteo_forecast = ensemble_prediction


def compute_forecast_flood_risk(
    forecast: xr.Dataset, flood_level: float, thredds: str = THREDDS_URL
) -> xr.Dataset:
    """Return the empirical exceedance probability for each forecast day based on a flood level threshold.

    Parameters
    ----------
    forecast : xr.Dataset
        Ensemble or deterministic streamflow forecast.
    flood_level : float
        Flood level threshold. Will be used to determine if forecasts exceed
        this specified flood threshold. Should be in the same units as the forecasted streamflow.
    thredds : str
        The thredds server url. Default: "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/".

    Returns
    -------
    xr.Dataset
        Time series of probabilities of flood level exceedance.
    """
    if thredds[-1] != "/":
        warnings.warn(
            "The thredds url should end with a slash. Appending it to the url."
        )
        thredds = f"{thredds}/"

    # ---- Calculations ---- #
    # Ensemble: for each day, calculate the percentage of members that are above the threshold
    if "member" in forecast.coords:
        # Get number of members originally
        number_members = len(forecast.member)

        # now compute the ratio of cases that are above the threshold
        pct = (
            forecast.where(forecast > flood_level)
            .notnull()
            .sum(dim="member", keep_attrs=True)
            / number_members
        )

    # it's deterministic:
    else:
        pct = (
            forecast.where(forecast > flood_level).notnull() / 1.0
        )  # This is needed to return values instead of floats

    domain = urlparse(thredds).netloc

    out = pct.to_dataset(name="exceedance_probability")
    out.attrs["source"] = f"PAVICS-Hydro flood risk forecasting tool, {domain}"
    out.attrs["history"] = (
        f"File created on {dt.datetime.now(dt.UTC).strftime('%Y-%m-%d %H:%M:%S')} "
        f"UTC on the PAVICS-Hydro service available at {domain}."
    )
    out.attrs["title"] = (
        "Identification of ensemble members that exceed a certain flow threshold."
    )

    return out
