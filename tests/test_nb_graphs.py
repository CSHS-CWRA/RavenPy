import pytest
import xarray as xr
from xclim import set_options


class TestNBGraphs:
    hydrographs = "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc"

    nbg = pytest.importorskip("ravenpy.utilities.nb_graphs")

    def test_hydrograph(self, yangtze):
        with xr.open_dataset(yangtze.fetch(self.hydrographs)) as ds:
            self.nbg.hydrographs(ds)

    def test_mean_annual_hydrograph(self, yangtze):
        with xr.open_dataset(yangtze.fetch(self.hydrographs)) as ds:
            self.nbg.mean_annual_hydrograph(ds)

    def test_spaghetti_annual_hydrograph(self, yangtze):
        with xr.open_dataset(yangtze.fetch(self.hydrographs)) as ds:
            self.nbg.spaghetti_annual_hydrograph(ds)

    def test_ts_fit_graph(self, yangtze):
        from xclim.indicators.generic import fit, stats

        with xr.open_dataset(yangtze.fetch(self.hydrographs)) as ds:
            ts = stats(ds.q_sim.load(), op="max", freq="ME")
        with set_options(check_missing="skip"):
            params = fit(ts, dist="gamma")
        self.nbg.ts_fit_graph(ts, params)
