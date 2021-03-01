"""Same graphics as in graphs, but using holoviews to generate web interactive graphics.

This is meant to be used in a notebook. In a console, use

hvplot.show(fig)

"""

import holoviews as hv
import hvplot.xarray
from holoviews.streams import Buffer
from bokeh.models import Range1d, LinearAxis
hv.extension('bokeh')


def hydrographs(ds):
    """Return a graphic showing the discharge simulations and observations."""

    basin_name = ds.basin_name.values[0]  # selected basin name
    g = ds.q_sim.hvplot.line(x="time", line_width=1,
                             label="Simulations",
                             ylabel="Streamflow (m³/s)",
                             xlabel="Time",
                             title=basin_name)

    # Plot the observed streamflows if available
    if hasattr(ds, "q_obs"):
        g *= ds.q_obs.hvplot.line(x="time", line_width=1.5, color="k", label="Observations")

    return g


def mean_annual_hydrograph(ds):
    """Return a graphic showing the discharge simulations and observations."""

    basin_name = ds.basin_name.values[0]  # selected basin name
    mq_sim = ds.q_sim.groupby("time.dayofyear").mean()

    g = mq_sim.hvplot.line(x="dayofyear",
                           line_width=1.5,
                           label="Mean simulation",
                           ylabel="Mean streamflow (m³/s)",
                           xlabel="Day of year",
                           title=basin_name)

    # Plot the observed streamflows if available
    if hasattr(ds, "q_obs"):
        mq_obs = ds.q_obs.groupby("time.dayofyear").mean()
        g *= mq_obs.hvplot.line(x="dayofyear",
                                line_width=2,
                                color="k",
                                label="Mean observations")

    return g


def spaghetti_annual_hydrograph(ds):
    """
    Create a spaghetti plot of the mean hydrological cycle for one model
    simulations.
    """
    basin_name = ds.basin_name.values[0]  # selected basin name
    mq_sim = ds.q_sim.groupby("time.dayofyear").mean()
    g1 = mq_sim.hvplot.line(x="dayofyear",
                           line_width=1.5,
                           label="Mean simulation",
                           ylabel="Mean streamflow (m³/s)",
                           xlabel="Day of year",
                           title=basin_name)

    g1 *= ds.q_sim.hvplot.line(x="time.dayofyear", by="time.year",
                             line_width=1,
                             label="Simulations",
                             legend=False,
                             )


    # Plot the observed streamflows if available
    if hasattr(ds, "q_obs"):
        mq_obs = ds.q_obs.groupby("time.dayofyear").mean()
        g2 = mq_obs.hvplot.line(x="dayofyear",
                                line_width=2,
                                color="k",
                                label="Mean observations")


        g2 *= ds.q_obs.hvplot.line(x="time.dayofyear", by="time.year",
                                line_width=1,
                                color="k",
                                label="Observations",
                                legend=False)


    return g1 + g2



# h_el = t.hvplot.hist(normed=False)
# c_el = hv.Curve(zip(q, pdf)).opts(framewise=True, hooks=[plot_secondary], color='red')
