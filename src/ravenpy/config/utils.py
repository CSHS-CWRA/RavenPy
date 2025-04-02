import os
import typing
from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union

import cf_xarray  # noqa: F401
import numpy as np
import xarray as xr

from .conventions import CF_RAVEN, MonthlyAverages
from .defaults import RAVEN_NO_DATA_VALUE


def nc_specs(
    fn: Union[str, os.PathLike[str]],
    data_type: str,
    station_idx: Optional[int] = None,
    alt_names: Optional[Union[str, Sequence[str]]] = None,
    mon_ave: bool = False,
    engine: str = "h5netcdf",
):
    """Extract specifications from netCDF file.

    Parameters
    ----------
    fn : str, Path
        NetCDF file path or DAP link.
    data_type : str
        Raven data type.
    station_idx : int, optional
        Index along station dimension. Starts at 1.
    alt_names : str or list of str, optional
        Alternative variable names for data type if not the CF standard default.
    mon_ave : bool
        If True, compute the monthly average.
    engine : str
        The engine used to open the dataset. Default is 'h5netcdf'.

    Returns
    -------
    dict

    Notes
    -----
    Attributes: var_name_nc, dim_names_nc, units, linear_transform, latitude_var_name_nc, longitude_var_name_nc,
    elevation_var_name_nc
    latitude, longitude, elevation, name
    """
    from ravenpy.utilities.coords import infer_scale_and_offset

    if isinstance(fn, str) and str(fn)[:4] == "http":
        pass
    elif Path(fn).exists():
        fn = os.path.realpath(fn, strict=True)
    else:
        raise ValueError("NetCDF file not found.")

    if alt_names is None:
        alt_names = ()
    elif isinstance(alt_names, str):
        alt_names = (alt_names,)

    attrs = {
        "file_name_nc": fn,
        "data_type": data_type,
        "forcing_type": data_type,
        "station_idx": station_idx,
    }

    with xr.open_dataset(fn, engine=engine) as ds:
        var_names = CF_RAVEN.get(data_type, ()) + tuple(alt_names)
        if len(var_names) == 0:
            raise ValueError(
                f"Provide alternative variable name with `alt_names` for {data_type}"
            )

        for v in var_names:
            if v in ds.data_vars:
                nc_var = ds[v]
                attrs["var_name_nc"] = v
                attrs["dim_names_nc"] = nc_var.dims
                attrs["_time_dim_name_nc"] = ds.cf["time"].name
                attrs["_dim_size_nc"] = dict(zip(nc_var.dims, nc_var.shape))
                attrs["units"] = nc_var.attrs.get("units")
                if attrs["units"] is not None:
                    s, o = infer_scale_and_offset(nc_var, data_type)
                    attrs["linear_transform"] = dict(scale=s, offset=o)
                if mon_ave:
                    ma = MonthlyAverages.get(data_type)
                    if ma:
                        attrs[ma] = nc_var.groupby("time.month").mean().values.tolist()

                break
        else:
            raise ValueError(f"No variable found for {data_type}.\n {ds.data_vars}")

        if station_idx is not None:
            # Convert to NumPy 0-based indexing
            i = station_idx - 1

            try:
                attrs["latitude_var_name_nc"] = ds.cf["latitude"].name
                attrs["longitude_var_name_nc"] = ds.cf["longitude"].name

                attrs["latitude"] = ds.cf["latitude"][i]
                attrs["longitude"] = ds.cf["longitude"][i]

            except KeyError:  # noqa: S110
                pass

            try:
                nc_elev = ds.cf["vertical"].name
            except KeyError:
                nc_elev = "elevation" if "elevation" in ds else None

            finally:
                if nc_elev is not None:
                    attrs["elevation_var_name_nc"] = nc_elev
                    try:
                        attrs["elevation"] = ds[nc_elev][i]
                    except IndexError:
                        # Elevation is a scalar
                        attrs["elevation"] = ds[nc_elev].item(0)

            if "station_id" in ds:
                if ds["station_id"].shape and len(ds["station_id"]) > i:
                    attrs["name"] = ds["station_id"].values[i]

    return attrs


def filter_for(kls, attrs, **kwds):
    """Return attributes that are fields of dataclass.

    Notes
    -----
    If attrs includes an attribute name and its Raven alias, e.g. `linear_transform` and `LinearTransform`, the latter will have priority.
    """
    from pydantic._internal._model_construction import ModelMetaclass

    from .commands import Command

    attrs.update(kwds)

    if hasattr(kls, "__dataclass_fields__"):
        return {
            k: v for (k, v) in attrs.items() if k in kls.__dataclass_fields__.keys()
        }
    else:
        out = {}
        for key, field in kls.model_fields.items():
            if field.alias in attrs:
                out[field.alias] = attrs[field.alias]
            elif key in attrs:
                out[key] = attrs[key]
            elif isinstance(field.annotation, ModelMetaclass) and issubclass(
                field.annotation, Command
            ):
                out[key] = filter_for(field.annotation, attrs)
        return out
        # return {k: v for (k, v) in attrs.items() if k in kls.__fields__.keys()}


def get_average_annual_runoff(
    nc_file_path: Union[str, os.PathLike[str]],
    area_in_m2: float,
    time_dim: str = "time",
    obs_var: str = "qobs",
    na_value: Union[int, float] = RAVEN_NO_DATA_VALUE,
):
    """Compute the average annual runoff from observed data."""
    with xr.open_dataset(nc_file_path) as ds:
        q_obs = ds.where(ds[obs_var] != na_value)[obs_var]
        q_obs *= 86400.0  # convert m**3/s to m**3/d
        axis = q_obs.dims.index(time_dim)
        # avg daily runoff [m3/d] for each year in record
        q_year = np.nanmean(q_obs.groupby("time.year").mean("time"), axis=axis)
        q_year = q_year / area_in_m2 * 365 * 1000.0  # [mm/yr] for each year in record
        q_year = np.mean(q_year)  # [mm/yr] mean over all years in record

    return q_year


def get_annotations(a):
    """Return all annotations inside [] or Union[...]."""
    for arg in typing.get_args(a):
        if typing.get_origin(arg) == Union:
            yield from get_annotations(arg)
        else:
            yield arg
