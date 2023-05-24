import pytest
import xarray as xr


class TestNBGraphs:
    hydrographs = "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc"

    nbg = pytest.importorskip("ravenpy.utilities.nb_graphs")

    def test_hydrograph(self, get_local_testdata):
        self.nbg.hydrographs(xr.open_dataset(get_local_testdata(self.hydrographs)))

    def test_mean_annual_hydrograph(self, get_local_testdata):
        self.nbg.mean_annual_hydrograph(
            xr.open_dataset(get_local_testdata(self.hydrographs))
        )

    def test_spaghetti_annual_hydrograph(self, get_local_testdata):
        self.nbg.spaghetti_annual_hydrograph(
            xr.open_dataset(get_local_testdata(self.hydrographs))
        )

    def test_ts_fit_graph(self, get_local_testdata):
        from xclim.indicators.generic import fit, stats

        ds = xr.open_dataset(get_local_testdata(self.hydrographs))

        ts = stats(ds.q_sim, op="max", freq="M")
        params = fit(ts, dist="gamma")
        self.nbg.ts_fit_graph(ts, params)
