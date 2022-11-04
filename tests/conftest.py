from pathlib import Path
from typing import Optional

import pytest
import xarray as xr
from filelock import FileLock
from xclim.indicators.land import fit, stats

from ravenpy.utilities.testdata import (
    _default_cache_dir,
    get_file,
    get_local_testdata,
    query_folder,
)


def populate_testing_data(
    temp_folder: Optional[Path] = None, _local_cache: Path = _default_cache_dir
):
    models = [
        "gr4j-cemaneige",
        "hbv-ec",
        "hmets",
        "mohyse",
    ]

    data_entries = dict()
    query_entries = dict()
    query_entries["ostrich-{model}-rv"] = {"ostrich-{model}": (r".rv\w", r"\.t\w{2}")}
    query_entries["raven-{model}-rv"] = {"raven-{model}": (r"\.rv\w",)}
    query_entries["raven-{model}-salmon_river"] = {
        "raven-{model}": (r"Salmon-River-Near-Prince-George_\w+.rvt",)
    }
    for model in models:
        for entry, query in query_entries.items():
            queried = []
            for folder, patterns in query.items():
                queried.extend(
                    [
                        query_folder(folder.format(model=model), pattern)
                        for pattern in patterns
                    ]
                )
            data_entries[entry.format(model=model)] = queried

    data_entries_simple = [
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc",
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
        "raven-gr4j-cemaneige/raven-gr4j-salmon.rvt",
        "gr4j_cemaneige/solution.rvc",
        "ostrich-gr4j-cemaneige/raven-gr4j-salmon.rvp.tpl",
        "ostrich-gr4j-cemaneige/OstRandomNumbers.txt",
        "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp85_nex-gddp_2070-2071_subset.nc",
        "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp45_nex-gddp_1971-1972_subset.nc",
        "nrcan/NRCAN_1971-1972_subset.nc",
        "watershed_vector/LSJ_LL.zip",
    ]

    data = dict()
    for filepattern in data_entries_simple:
        if temp_folder is None:
            data[filepattern] = get_file(filepattern, cache_dir=_local_cache)
        else:
            data[filepattern] = get_local_testdata(
                filepattern, temp_folder=temp_folder, _local_cache=_local_cache
            )

    for name, filepattern in data_entries.items():
        if temp_folder is None:
            data[name] = [
                get_file(f, cache_dir=_default_cache_dir) for f in filepattern
            ]
        else:
            data[name] = [
                get_local_testdata(
                    f, temp_folder=temp_folder, _local_cache=_local_cache
                )
                for f in filepattern
            ]

    return data


@pytest.fixture(scope="session")
def threadsafe_data_dir(tmp_path_factory) -> Path:
    """Constructor for worker-session temporary data folders"""
    return Path(tmp_path_factory.getbasetemp().joinpath("data"))


@pytest.fixture(scope="session")
def session_data(threadsafe_data_dir, worker_id):

    test_data_being_written = FileLock(_default_cache_dir.joinpath(".lock"))

    with test_data_being_written as fl:
        if worker_id in ["master", "gw0"]:
            return populate_testing_data()
        fl.acquire()
        return populate_testing_data(temp_folder=threadsafe_data_dir)


@pytest.fixture(scope="session")
def q_sim_1(threadsafe_data_dir):
    """A file storing a Raven streamflow simulation over one basin."""
    return get_file(
        "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
        cache_dir=threadsafe_data_dir,
    )


@pytest.fixture(scope="session")
def ts_stats(q_sim_1, tmp_path):
    q = xr.open_dataset(q_sim_1).q_sim
    ts = stats(q, op="max")
    fn = tmp_path / "ts_stats.nc"
    ts.to_netcdf_(fn)
    return fn


@pytest.fixture(scope="session")
def params(ts_stats, tmp_path):
    ds = xr.open_dataset(ts_stats)
    name = list(ds.data_vars.keys()).pop()
    q = ds[name]
    p = fit(q, dist="gumbel_r")
    fn = tmp_path / "fit.nc"
    p.to_netcdf(fn)
    return fn
