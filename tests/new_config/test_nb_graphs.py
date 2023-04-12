import pytest

from ravenpy.utilities.testdata import open_dataset


class TestNBGraphs:
    hydrographs = "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc"

    nbg = pytest.importorskip("ravenpy.utilities.nb_graphs")

    def test_hydrograph(self, threadsafe_data_dir):
        self.nbg.hydrographs(
            open_dataset(self.hydrographs, cache_dir=threadsafe_data_dir)
        )

    def test_mean_annual_hydrograph(self, threadsafe_data_dir):
        self.nbg.mean_annual_hydrograph(
            open_dataset(self.hydrographs, cache_dir=threadsafe_data_dir)
        )

    def test_spaghetti_annual_hydrograph(self, threadsafe_data_dir):
        self.nbg.spaghetti_annual_hydrograph(
            open_dataset(self.hydrographs, cache_dir=threadsafe_data_dir)
        )

    def test_ts_fit_graph(self, threadsafe_data_dir):
        from xclim.indicators.land import fit, stats

        ds = open_dataset(self.hydrographs, cache_dir=threadsafe_data_dir)

        ts = stats(ds.q_sim, op="max", freq="M")
        params = fit(ts, dist="gamma")
        self.nbg.ts_fit_graph(ts, params)
