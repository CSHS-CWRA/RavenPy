"""
Library to perform graphs for the streamflow time series analysis.

The following graphs can be plotted:
    - hydrograph
    - mean_annual_hydrograph
    - spaghetti_annual_hydrograph
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Union

import matplotlib.pyplot
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import pyplot as plt
from scipy import stats
from xclim.core.units import units2pint

from ravenpy.utilities.mk_test import mk_test_calc

# TODO: Review docstrings, ensure numpydoc convention compliance


def hydrograph(file_list: Sequence[Union[str, Path]]):
    """Create a graphic of the hydrograph for each model simulation.

    Parameters
    ----------
    file_list : Sequence of str or Path
        Raven output files containing simulated streamflows.
    """
    ds = [xr.open_dataset(file) for file in file_list]

    # Get time data for the plot
    dates = pd.DatetimeIndex(ds[0].time.values)
    first_date = dates.min().strftime("%Y/%m/%d")
    last_date = dates.max().strftime("%Y/%m/%d")

    basin_name = ds[0].basin_name.values[0]  # selected basin name

    fig, ax = plt.subplots()  # initialize figure

    # Plot the observed streamflows if available
    if hasattr(ds[0], "q_obs"):
        q_obs = ds[0].q_obs
        plt.plot(dates, q_obs, linewidth=2, label="obs")

    # Plot the simulated streamflows for each hydrological model
    q_sim = [h.q_sim for h in ds]
    for sim in q_sim:
        plt.plot(dates, sim, linewidth=2, label="sim: " + basin_name)

    # plt.xlim([first_date, last_date])
    plt.ylim(bottom=0, top=None)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$Streamflow [m^3s^{{-1}}]$")
    ax.set_title(
        f"Hydrograph between {first_date} and {last_date}\nSelected basin: {basin_name} basin."
    )
    ax.legend()
    ax.grid()

    plt.xticks(rotation=90)
    plt.tight_layout()

    return fig


def mean_annual_hydrograph(file_list: Sequence[Union[str, Path]]):
    """
    Create a graphic of the mean hydrological cycle for each model simulation.

    Parameters
    ----------
    file_list : Sequence of str or Path
        Raven output files containing simulated streamflows.
    """
    # Time series for the plot
    ds = [xr.open_dataset(file) for file in file_list]

    # Get time data for the plot
    dates = pd.DatetimeIndex(ds[0].time.values)
    first_date = dates.min().strftime("%Y/%m/%d")
    last_date = dates.max().strftime("%Y/%m/%d")

    basin_name = ds[0].basin_name.values[0]  # selected basin name

    fig, ax = plt.subplots()  # initialize figure

    # Plot the observed streamflows if available
    mah_obs = None
    if hasattr(ds[0], "q_obs"):
        q_obs = ds[0].q_obs
        mah_obs = q_obs.groupby("time.dayofyear").mean()

        plt.plot(mah_obs.dayofyear, mah_obs, linewidth=2, label="obs")

    # Plot the simulated streamflows for each hydrological model
    q_sim = [h.q_sim for h in ds]
    mahs_sim = [q.groupby("time.dayofyear").mean() for q in q_sim]

    for mah in mahs_sim:
        plt.plot(mah.dayofyear, mah, linewidth=2, label="sim: " + basin_name)

    plt.xticks(
        np.linspace(0, 365, 13)[:-1],
        (
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec",
        ),
    )

    if mah_obs is not None:
        plt.xlim(0, mah_obs.shape[0])
    plt.ylim(bottom=0, top=None)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$Streamflow [m^3s^{{-1}}]$")
    ax.set_title(
        f"Hydrograph between {first_date} and {last_date}\nSelected basin: {basin_name} basin."
    )
    ax.legend()
    ax.grid()

    plt.tight_layout()

    return fig


def spaghetti_annual_hydrograph(file: Union[str, Path]):
    """Create a spaghetti plot of the mean hydrological cycle for one model simulations.

    The mean simulation is also displayed.

    Parameters
    ----------
    file : str or Path
        Raven output files containing simulated streamflows of one model.
    """
    # Time series for the plot
    ds = xr.open_dataset(file)

    # Get time data for the plot
    dates = pd.DatetimeIndex(ds.time.values)
    first_date = dates.min().strftime("%Y/%m/%d")
    last_date = dates.max().strftime("%Y/%m/%d")

    basin_name = ds.basin_name.values[0]  # selected basin name

    fig, ax = plt.subplots()  # initialize figure

    # Plot the observed streamflows if available
    if hasattr(ds, "q_obs"):
        q_obs = ds.q_obs
        mah_obs = q_obs.groupby("time.year")
        mah_obs_mean = q_obs.groupby("time.dayofyear").mean()

        for year in mah_obs.groups.keys():
            plt.plot(
                np.arange(1, q_obs.values[mah_obs.groups[year]].shape[0] + 1, 1),
                q_obs.values[mah_obs.groups[year]],
                linewidth=1,
                color="C0",
            )

        plt.plot(
            mah_obs_mean.dayofyear, mah_obs_mean, linewidth=2, color="C0", label="obs"
        )

    # Plot the simulated streamflows for each hydrological model
    q_sim = ds.q_sim
    mah_sim = q_sim.groupby("time.year")
    mah_sim_mean = q_sim.groupby("time.dayofyear").mean()

    for year in mah_sim.groups.keys():
        plt.plot(
            np.arange(1, q_sim.values[mah_sim.groups[year]].shape[0] + 1, 1),
            q_sim.values[mah_sim.groups[year]],
            linewidth=1,
            color="C1",
        )

    plt.plot(
        mah_sim_mean.dayofyear,
        mah_sim_mean,
        linewidth=2,
        color="C1",
        label="sim: " + "<model_name>",
    )

    plt.xticks(
        np.linspace(0, 365, 13)[:-1],
        (
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec",
        ),
    )

    plt.xlim(0, mah_sim_mean.shape[0])
    plt.ylim(bottom=0, top=None)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$Streamflow [m^3s^{{-1}}]$")
    ax.set_title(
        f"Spaghetti annual hydrograph between {first_date} and {last_date}"
        f"\n Selected basin: {basin_name} basin."
    )
    ax.legend()
    ax.grid()

    plt.tight_layout()

    return fig


def ts_graphs(file, trend: bool = True, alpha: float = 0.05):
    """Create a figure with the statistics so one can see a trend in the data.

    Graphs for time series statistics.

    Parameters
    ----------
    file : str or Path
        xarray-compatible file containing streamflow statistics for one run.
    trend : bool
        If True, the slope will be plotted.
    alpha : float
        Significance level for the Mann-Kendall test.
    """
    # Time series for the plot
    ds = xr.open_dataset(file)
    titlename = ds["ts_stats"].description
    values = ds["ts_stats"].values[:]

    values[0] = 2
    values[1] = 1
    values[2] = 3
    values[3] = 3
    values[4] = 3
    values[5] = 5
    values[6] = 6
    values[7] = 7
    values[8] = 10
    values[9] = 4

    # Get time data for the plot
    dates = pd.DatetimeIndex(ds.time.values)
    # first_date = dates.min()#.strftime('%Y/%m/%d')
    # last_date = dates.max()#.strftime('%Y/%m/%d')
    res = None

    if trend:
        res = stats.theilslopes(values, dates)
        trd, h, p, z = mk_test_calc(values, alpha=alpha)
        titlename = (
            titlename
            + ", Mann-Kendall h="
            + str(h)
            + ", p-value="
            + str(np.round(p, 4))
        )

    fig, ax = plt.subplots()
    ax.plot(dates, values, label="time-series index")

    # plt.xlim([first_date, last_date])
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$Streamflow [m^3s^{{-1}}]$")

    ax.set_title(titlename)

    ax.grid()
    plt.tight_layout()

    if trend:
        # TODO: This does not work yet, trying to compute the y-value of a datetime x-axis value * a slope...
        ax.plot(
            dates,
            res[1] + res[0] * dates,
            linestyle="--",
            linewidth=2,
            label="Sen's Slope",
        )

    ax.legend()

    return fig


def ts_fit_graph(ts: xr.DataArray, params: xr.DataArray) -> matplotlib.pyplot.Figure:
    """Create graphic showing a histogram of the data and the distribution fitted to it.

    The graphic contains one panel per watershed.

    Parameters
    ----------
    ts : xr.DataArray
        Stream flow time series with dimensions (time, nbasins).
    params : xr.DataArray
        Fitted distribution parameters returned by `xclim.land.fit` indicator.

    Returns
    -------
    matplotlib.pyplot.Figure
        Figure showing a histogram and the parameterized pdf.
    """
    from xclim.indices.stats import get_dist

    n = ts.nbasins.size
    dist = params.attrs["scipy_dist"]

    fig, axes = plt.subplots(n, figsize=(10, 6), squeeze=False)
    if params.isnull().any():
        raise ValueError("Null values found in `params`.")

    for i in range(n):
        ax = axes.flat[i]
        ax2 = plt.twinx(ax)
        p = params.isel(nbasins=i)
        t = ts.isel(nbasins=i).dropna(dim="time")

        # Plot histogram of time series as density then as a normal count.
        density, bins, patches = ax.hist(
            t,
            alpha=0.5,
            density=True,
            bins="auto",
            label="__nolabel__",
        )
        ax2.hist(
            t,
            facecolor="none",
            bins=bins,
            edgecolor="grey",
            linewidth=1,
        )

        # Plot pdf of distribution
        dc = get_dist(dist)(*params.isel(nbasins=i))
        mn = dc.ppf(0.01)
        mx = dc.ppf(0.99)
        q = np.linspace(mn, mx, 200)
        pdf = dc.pdf(q)

        ps = ", ".join([f"{x:.1f}" for x in p.values])
        ax.plot(q, pdf, "-", label=f"{params.attrs['scipy_dist']}({ps})")

        # Labels
        ax.set_xlabel(f"{ts.long_name} (${units2pint(ts.units):~P}$)")
        ax.set_ylabel("Probability density")
        ax2.set_ylabel("Histogram count")

        ax.legend(frameon=False)

    plt.tight_layout()
    return fig


def forecast(
    file: Union[str, Path], fcst_var: str = "q_sim"
) -> matplotlib.pyplot.Figure:
    """Create a graphic of the hydrograph for each forecast member.

    Parameters
    ----------
    file : str or Path
        Raven output file containing simulated streamflows.
    fcst_var : str
        Name of the streamflow variable.

    Returns
    -------
    matplotlib.pyplot.Figure
    """
    ds = xr.open_dataset(file)

    # Get time data for the plot
    dates = pd.DatetimeIndex(ds.time.values)
    start = dates.min()
    end = dates.max()

    fig, ax = plt.subplots()  # initialize figure

    # Plot the simulated streamflows for each hydrological model
    ds[fcst_var].plot.line("b", x="time", add_legend=False)

    # plt.xlim([first_date, last_date])
    ax.set_xlabel("Time")
    ax.set_ylabel("Streamflow (m³/s)")
    ax.set_title(f"Forecasted hydrograph between {start:%Y/%m/%d} and {end:%Y/%m/%d}.")
    ax.grid()

    plt.xticks(rotation=90)
    plt.tight_layout()

    return fig


def hindcast(
    file: Union[str, Path], fcst_var: str, qobs: Union[str, Path], qobs_var: str
) -> matplotlib.pyplot.Figure:
    """Create a graphic of the hydrograph for each hindcast member.

    Parameters
    ----------
    file : str or Path
        Raven output file containing simulated streamflows.
    fcst_var : str
        Name of the streamflow variable.
    qobs : str or Path
        Streamflow observation file, with times matching the hindcast.
    qobs_var : str
        Name of the streamflow observation variable.

    Returns
    -------
    matplotlib.pyplot.Figure
    """
    ds = xr.open_dataset(file)
    ds2 = xr.open_dataset(qobs)

    # Get time data for the plot
    dates = pd.DatetimeIndex(ds.time.values)
    start = dates.min()
    end = dates.max()

    fig, ax = plt.subplots()  # initialize figure

    # Plot the simulated streamflows for each hydrological model
    hh = ds[fcst_var].plot.line("b", x="time", label="Hindcasts")
    ho = ds2[qobs_var].plot.line("r", label="Observations")

    # plt.xlim([first_date, last_date])
    ax.set_xlabel("Time")
    ax.set_ylabel("Streamflow (m³/s)")
    ax.set_title(f"Hindcasted hydrograph between {start:%Y/%m/%d} and {end:%Y/%m/%d}.")

    # Add legend
    handles = hh[:1] + ho
    ax.legend(handles, [h.get_label() for h in handles])
    ax.grid()

    plt.xticks(rotation=90)
    plt.tight_layout()

    return fig
