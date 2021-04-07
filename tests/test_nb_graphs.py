import xarray as xr

from ravenpy.utilities import nb_graphs as nbg
from ravenpy.utilities.testdata import get_local_testdata

fn = get_local_testdata(
    "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc"
)
ds = xr.open_dataset(fn)


def test_hydrograph():
    nbg.hydrographs(ds)


def test_mean_annual_hydrograph():
    nbg.mean_annual_hydrograph(ds)


def test_spaghetti_annual_hydrograph():
    nbg.spaghetti_annual_hydrograph(ds)


def test_ts_fit_graph():
    from xclim.indicators.land import fit, stats

    ts = stats(ds.q_sim, op="max", freq="M")
    params = fit(ts, dist="gamma")
    nbg.ts_fit_graph(ts, params)
