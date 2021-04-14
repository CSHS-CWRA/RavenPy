"""
This module contains functions creating web-friendly interactive graphics using holoviews.

The graphic outputs are meant to be displayed in a notebook.
In a console, use `hvplot.show(fig)` to render the figures.
"""
from pathlib import Path

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

try:
    import holoviews as hv
    import hvplot.xarray

    # from holoviews.streams import Buffer
    # from bokeh.models import Range1d, LinearAxis
except (ImportError, ModuleNotFoundError) as e:
    raise ImportError("Install holoviews and hvplot.") from e

hv.extension("bokeh")


def hydrographs(ds):
    """Return a graphic showing the discharge simulations and observations."""

    basin_name = ds.basin_name.values[0]  # selected basin name
    g = ds.q_sim.hvplot.line(
        x="time",
        line_width=1,
        label="Simulations",
        ylabel="Streamflow (m³/s)",
        xlabel="Time",
        title=basin_name,
    )

    # Plot the observed streamflows if available
    if hasattr(ds, "q_obs"):
        g *= ds.q_obs.hvplot.line(
            x="time", line_width=1.5, color="k", label="Observations"
        )

    return g


def mean_annual_hydrograph(ds):
    """Return a graphic showing the discharge simulations and observations."""

    basin_name = ds.basin_name.values[0]  # selected basin name
    mq_sim = ds.q_sim.groupby("time.dayofyear").mean()

    g = mq_sim.hvplot.line(
        x="dayofyear",
        line_width=1.5,
        label="Mean simulation",
        ylabel="Mean streamflow (m³/s)",
        xlabel="Day of year",
        title=basin_name,
    )

    # Plot the observed streamflows if available
    if hasattr(ds, "q_obs"):
        mq_obs = ds.q_obs.groupby("time.dayofyear").mean()
        g *= mq_obs.hvplot.line(
            x="dayofyear", line_width=2, color="k", label="Mean observations"
        )

    return g


def spaghetti_annual_hydrograph(ds):
    """
    Create a spaghetti plot of the mean hydrological cycle for one model
    simulations.
    """
    basin_name = ds.basin_name.values[0]  # selected basin name
    mq_sim = ds.q_sim.groupby("time.dayofyear").mean()
    g1 = mq_sim.hvplot.line(
        x="dayofyear",
        line_width=1.5,
        label="Mean simulation",
        ylabel="Mean streamflow (m³/s)",
        xlabel="Day of year",
        title=basin_name,
    )

    g1 *= ds.q_sim.hvplot.line(
        x="time.dayofyear",
        by="time.year",
        line_width=1,
        label="Simulations",
        legend=False,
    )

    # Plot the observed streamflows if available
    if hasattr(ds, "q_obs"):
        mq_obs = ds.q_obs.groupby("time.dayofyear").mean()
        g2 = mq_obs.hvplot.line(
            x="dayofyear", line_width=2, color="k", label="Mean observations"
        )

        g2 *= ds.q_obs.hvplot.line(
            x="time.dayofyear",
            by="time.year",
            line_width=1,
            color="k",
            label="Observations",
            legend=False,
        )

    return g1 + g2


def ts_fit_graph(ts, params):
    """Create graphic showing an histogram of the data and the distribution fitted to it.

    The graphic contains one panel per watershed.

    Parameters
    ----------
    ts : xr.DataArray
      Stream flow time series with dimensions (time, nbasins).
    params : xr.DataArray
      Fitted distribution parameters returned by `xclim.land.fit` indicator.

    Returns
    -------
    fig
      Figure showing a histogram and the parameterized pdf.
    """
    # Note: The hover tool could be customized to show the histogram count in addition to the frequency.
    from xclim.indices.stats import get_dist

    n = ts.nbasins.size
    if n > 1:
        raise NotImplementedError

    ts = ts.isel(nbasins=0)
    params = params.isel(nbasins=0)

    # Using matplotlib's default binning strategy
    hist, bins, mh = plt.hist(ts, bins="auto", density=True)

    # Histogram graphic object
    h = hv.Histogram((hist, bins), kdims=ts.name, label="Histogram")

    # PDF domain
    mn = np.min(bins)
    mx = np.max(bins)
    q = np.linspace(mn, mx, 200)

    # Compute PDF
    dist = params.attrs["scipy_dist"]
    dc = get_dist(dist)(*params)  # Works because dparams is the first dimension.
    pdf = xr.DataArray(
        data=dc.pdf(q),
        dims=(ts.name,) + params.dims[1:],
        coords={ts.name: q},
        name="pdf",
    )

    # PDF line label
    ps = ", ".join(
        [
            f"{key}={x:.1f}".format(x)
            for (key, x) in zip(params.dparams.data, params.values)
        ]
    )
    label = "{}({})".format(dist, ps)

    # PDF graphic object
    p = pdf.hvplot.line(label=label, xlabel=ts.attrs["long_name"], color="orange")

    # Layout
    return (h * p).opts(hv.opts.Histogram(tools=["hover"]))
