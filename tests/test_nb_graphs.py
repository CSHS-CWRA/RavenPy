import pytest
import xarray as xr

from ravenpy.utilities.testdata import get_local_testdata


class TestNBGraphs:
    fn = get_local_testdata(
        "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc"
    )
    ds = xr.open_dataset(fn)
    nbg = pytest.importorskip("ravenpy.utilities.nb_graphs")

    def test_hydrograph(self):
        self.nbg.hydrographs(self.ds)

    def test_mean_annual_hydrograph(self):
        self.nbg.mean_annual_hydrograph(self.ds)

    def test_spaghetti_annual_hydrograph(self):
        self.nbg.spaghetti_annual_hydrograph(self.ds)

    def test_ts_fit_graph(self):
        from xclim.indicators.land import fit, stats

        ts = stats(self.ds.q_sim, op="max", freq="M")
        params = fit(ts, dist="gamma")
        self.nbg.ts_fit_graph(ts, params)
