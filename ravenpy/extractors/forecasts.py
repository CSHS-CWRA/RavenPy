import datetime as dt
import logging
import re
from pathlib import Path
from typing import Any, List, Tuple, Union

import pandas as pd
import xarray as xr
from pandas import DatetimeIndex, Series, Timestamp
from xarray import Dataset

from . import gis_import_error_message

try:
    import fiona
except (ImportError, ModuleNotFoundError) as e:
    msg = gis_import_error_message.format(Path(__file__).stem)
    raise ImportError(msg) from e

LOGGER = logging.getLogger("PYWPS")


def get_hindcast_day(region_coll: fiona.Collection, date, climate_model="GEPS"):
    """Generate a forecast dataset that can be used to run raven.

    Data comes from the CASPAR archive and must be aggregated such that each file
    contains forecast data for a single day, but for all forecast timesteps and
    all members.

    The code takes the region shapefile, the forecast date required, and the
    climate_model to use, here GEPS by default, but eventually could be GEPS, GDPS, REPS or RDPS.
    """

    # Get the file locations and filenames as a function of the climate model and date
    [ds, times] = get_CASPAR_dataset(climate_model, date)

    return get_subsetted_forecast(region_coll, ds, times, True)


def get_CASPAR_dataset(
    climate_model: str, date: dt.datetime
) -> Tuple[
    xr.Dataset, List[Union[Union[DatetimeIndex, Series, Timestamp, Timestamp], Any]]
]:
    """Return CASPAR dataset."""

    if climate_model == "GEPS":
        d = dt.datetime.strftime(date, "%Y%m%d")
        file_url = f"https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/disk2/caspar/daily/GEPS_{d}.nc"
        ds = xr.open_dataset(file_url)
        # Here we also extract the times at 6-hour intervals as Raven must have
        # constant timesteps and GEPS goes to 6 hours
        start = pd.to_datetime(ds.time[0].values)
        times = [start + dt.timedelta(hours=n) for n in range(0, 384, 6)]
    else:
        # Eventually: GDPS, RDPS and REPS
        raise NotImplementedError("Only the GEPS model is currently supported")

    # Checking that these exist.
    for f in ["pr", "tas"]:
        if f not in ds:
            raise AttributeError(f"'{f}' not present in dataset")

    return ds, times


def get_ECCC_dataset(
    climate_model: str,
) -> Tuple[
    Dataset, List[Union[Union[DatetimeIndex, Series, Timestamp, Timestamp], Any]]
]:
    """Return latest GEPS forecast dataset."""
    if climate_model == "GEPS":
        # Eventually the file will find a permanent home, until then let's use the test folder.
        file_url = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/forecasts/eccc_geps/GEPS_latest.ncml"

        ds = xr.open_dataset(file_url)
        # Here we also extract the times at 6-hour intervals as Raven must have
        # constant timesteps and GEPS goes to 6 hours
        start = pd.to_datetime(ds.time[0].values)
        times = [start + dt.timedelta(hours=n) for n in range(0, 384, 6)]
    else:
        # Eventually: GDPS, RDPS and REPS
        raise NotImplementedError("Only the GEPS model is currently supported")

    # Checking that these exist. IF the files are still processing, possible that one or both are not available!
    for f in ["pr", "tas"]:
        if f not in ds:
            raise AttributeError(f"'{f}' not present in dataset")

    return ds, times


def get_recent_ECCC_forecast(
    region_coll: fiona.Collection, climate_model: str = "GEPS"
) -> xr.Dataset:
    """Generate a forecast dataset that can be used to run raven.

    Data comes from the ECCC datamart and collected daily. It is aggregated
    such that each file contains forecast data for a single day, but for all
    forecast timesteps and all members.

    The code takes the region shapefile and the climate_model to use, here GEPS
    by default, but eventually could be GEPS, GDPS, REPS or RDPS.

    Parameters
    ----------
    region_coll : fiona.Collection
        The region vectors.
    climate_model : str
        Type of climate model, for now only "GEPS" is supported.

    Returns
    -------
    xr.Dataset
        The forecast dataset.
    """

    [ds, times] = get_ECCC_dataset(climate_model)

    # Make the variable name compatible with the hindcasting tools.
    ds = ds.rename({"member": "members"})

    return get_subsetted_forecast(region_coll, ds, times, False)


def get_subsetted_forecast(
    region_coll: fiona.Collection,
    ds: xr.Dataset,
    times: Union[dt.datetime, xr.DataArray],
    is_caspar: bool,
) -> xr.Dataset:
    """
    This function takes a dataset, a region and the time sampling array and returns
    the subsetted values for the given region and times

    Parameters
    ----------
    region_coll : fiona.Collection
        The region vectors.
    ds : xr.Dataset
        The dataset containing the raw, worldwide forecast data
    times : dt.datetime or xr.DataArray
        The array of times required to do the forecast.
    is_caspar : bool
        True if the data comes from Caspar, false otherwise.
        Used to define lat/lon on rotated grid.

    Returns
    -------
    xr.Dataset
        The forecast dataset.

    """
    # Extract the bounding box to subset the entire forecast grid to something
    # more manageable
    lon_min = region_coll.bounds[0]
    lon_max = region_coll.bounds[2]
    lat_min = region_coll.bounds[1]
    lat_max = region_coll.bounds[3]

    # Add a very simple lon wraparound if data suggests for it
    if ((ds.lon.min() >= 0) and (ds.lon.max() <= 360)) and (lon_max < 0):
        lon_min += 360
        lon_max += 360

    # Subset the data to the desired location (bounding box) and times
    ds = ds.where(
        (ds.lon <= lon_max)
        & (ds.lon >= lon_min)
        & (ds.lat <= lat_max)
        & (ds.lat >= lat_min),
        drop=True,
    ).sel(time=times)

    # Rioxarray requires CRS definitions for variables
    # Get CRS, e.g. 4326
    crs = int(re.match(r"epsg:(\d+)", region_coll.crs["init"]).group(1))

    # Here the name of the variable could differ based on the Caspar file processing
    tas = ds.tas.rio.write_crs(crs)
    pr = ds.pr.rio.write_crs(crs)
    ds = xr.merge([tas, pr])

    # Now apply the mask of the basin contour and average the values to get a single time series
    if is_caspar:
        ds.rio.set_spatial_dims("rlon", "rlat")
        ds["rlon"] = ds["rlon"] - 360
        # clip the netcdf and average across space.
        shdf = [next(iter(region_coll))["geometry"]]
        forecast = ds.rio.clip(shdf, crs=crs)
        forecast = forecast.mean(dim={"rlat", "rlon"}, keep_attrs=True)

    else:
        ds.rio.set_spatial_dims("lon", "lat")
        ds["lon"] = ds["lon"] - 360
        # clip the netcdf and average across space.
        shdf = [next(iter(region_coll))["geometry"]]
        forecast = ds.rio.clip(shdf, crs=crs)
        forecast = forecast.mean(dim={"lat", "lon"}, keep_attrs=True)

    return forecast
