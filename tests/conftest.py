from pathlib import Path

import pytest
import xarray as xr
from xclim.indicators.land import fit, stats

from ravenpy.utilities.testdata import get_file


@pytest.fixture(autouse=True, scope="session")
def threadsafe_data_dir(tmp_path_factory) -> Path:
    return Path(tmp_path_factory.getbasetemp().joinpath("data"))


@pytest.fixture
def q_sim_1(tmp_path, threadsafe_data_dir):
    """A file storing a Raven streamflow simulation over one basin."""
    return get_file(
        "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
        cache_dir=threadsafe_data_dir,
    )


@pytest.fixture
def ts_stats(q_sim_1, tmp_path):
    q = xr.open_dataset(q_sim_1).q_sim
    ts = stats(q, op="max")
    fn = tmp_path / "ts_stats.nc"
    ts.to_netcdf_(fn)
    return fn


@pytest.fixture
def params(ts_stats, tmp_path):
    ds = xr.open_dataset(ts_stats)
    name = list(ds.data_vars.keys()).pop()
    q = ds[name]
    p = fit(q, dist="gumbel_r")
    fn = tmp_path / "fit.nc"
    p.to_netcdf(fn)
    return fn
