import numpy as np
import xarray as xr
from xclim.indicators.generic import fit, stats

from ravenpy.utilities import graphs


class TestGraph:
    def test_ts_fit_graph(self, get_local_testdata):
        f = get_local_testdata(
            "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
        )
        with xr.open_dataset(f) as ds:
            ts = stats(ds.q_sim, op="max", freq="M")
            p = fit(ts)
            np.testing.assert_array_equal(p.isnull(), False)

            fig = graphs.ts_fit_graph(ts, p)
            return fig
