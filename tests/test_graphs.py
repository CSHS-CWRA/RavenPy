import numpy as np
import xarray as xr
from xclim.indicators.land import fit, stats

from ravenpy.utilities import graphs
from ravenpy.utilities.testdata import open_dataset


class TestGraph:
    def test_ts_fit_graph(self, threadsafe_data_dir):
        ds = open_dataset(
            "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
            cache_dir=threadsafe_data_dir,
        )

        ts = stats(ds.q_sim, op="max", freq="M")
        p = fit(ts)
        np.testing.assert_array_equal(p.isnull(), False)

        fig = graphs.ts_fit_graph(ts, p)
        return fig
