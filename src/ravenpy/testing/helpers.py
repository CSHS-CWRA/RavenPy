from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

__all__ = [
    "count_pixels",
    "make_bnds",
    "synthetic_gr4j_inputs",
]


def count_pixels(stats: dict, numeric_categories=False) -> int:
    category_counts = 0
    for key, val in stats.items():
        if numeric_categories:
            try:
                int(key)
            except ValueError:  # noqa: S112
                continue
        if key in ["count", "min", "max", "mean", "median", "sum", "nodata"]:
            continue
        category_counts += val
    return category_counts


def make_bnds(params, delta):
    """
    Return low and high parameter bounds by subtracting and adding delta*params to params.

    Parameters
    ----------
    params : sequence
        Parameters.
    delta : float [0,1]
        Relative delta to subtract and add to parameters.

    Returns
    -------
    (tuple, tuple)
        Low and high bounds for parameters.
    """
    arr = np.asarray(params)
    d = np.abs(arr * delta)
    return tuple(arr - d), tuple(arr + d)


def synthetic_gr4j_inputs(path):
    time = pd.date_range(start="2000-07-01", end="2002-07-01", freq="D")

    pr = 3 * np.ones(len(time))
    pr = xr.DataArray(pr, coords={"time": time}, dims="time", name="pr")
    pr.to_netcdf(Path(path).joinpath("pr.nc"))

    tas = 280 + 20 * np.cos(np.arange(len(time)) * 2 * np.pi / 365.0)
    tas = xr.DataArray(tas, coords={"time": time}, dims="time", name="tas")
    tas.to_netcdf(Path(path).joinpath("tas.nc"))

    evap = 3 + 3 * np.cos(-30 + np.arange(len(time)) * 2 * np.pi / 365.0)
    evap = xr.DataArray(evap, coords={"time": time}, dims="time", name="evap")
    evap.to_netcdf(Path(path).joinpath("evap.nc"))
