from shutil import copyfile

import numpy as np
import xarray as xr
from xclim import set_options
from xclim.indicators.generic import fit, stats

from ravenpy.utilities import graphs


class TestGraph:
    def test_ts_fit_graph(self, get_local_testdata, tmp_path):
        raven_hydrograph = get_local_testdata(
            "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
        )
        file = tmp_path / "raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc"

        copyfile(raven_hydrograph, file)

        with xr.open_dataset(file) as ds:
            ts = stats(ds.q_sim, op="max", freq="ME")
            with set_options(check_missing="skip"):
                p = fit(ts)
            np.testing.assert_array_equal(p.isnull(), False)

            fig = graphs.ts_fit_graph(ts, p)
            return fig
