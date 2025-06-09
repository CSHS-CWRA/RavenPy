from shutil import copyfile

import numpy as np
from xarray import open_dataset
from xclim import set_options
from xclim.indicators.generic import fit, stats

from ravenpy.utilities import graphs


class TestGraph:
    def test_ts_fit_graph(self, tmp_path, yangtze):
        raven_hydrograph = yangtze.fetch(
            "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
        )
        file = tmp_path / "raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc"

        copyfile(raven_hydrograph, file)

        with open_dataset(file) as ds:
            ts = stats(ds.q_sim, op="max", freq="ME")
            with set_options(check_missing="skip"):
                p = fit(ts)
            np.testing.assert_array_equal(p.isnull(), False)

            fig = graphs.ts_fit_graph(ts, p)
            assert fig
