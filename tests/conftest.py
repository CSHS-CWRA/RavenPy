from pathlib import Path
from typing import Optional

import pytest
import xarray as xr
from filelock import FileLock
from xclim.indicators.land import fit, stats

from ravenpy.utilities.testdata import _default_cache_dir, get_file, get_local_testdata


def populate_testing_data(
    temp_folder: Optional[Path] = None,
    branch: str = "master",
    _local_cache: Path = _default_cache_dir,
):
    if _local_cache.joinpath(".data_written").exists():
        # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
        return

    models = [
        "gr4j-cemaneige",
        "hbv-ec",
        "hmets",
        "mohyse",
    ]

    data_entries = list()
    entries = [
        "ostrich-{model}/raven-{model}-salmon.rvc",
        "ostrich-{model}/raven-{model}-salmon.rvc.tpl",
        "ostrich-{model}/raven-{model}-salmon.rvh",
        "ostrich-{model}/raven-{model}-salmon.rvh.tpl",
        "ostrich-{model}/raven-{model}-salmon.rvi",
        "ostrich-{model}/raven-{model}-salmon.rvi.tpl",
        "ostrich-{model}/raven-{model}-salmon.rvp",
        "ostrich-{model}/raven-{model}-salmon.rvp.tpl",
        "ostrich-{model}/raven-{model}-salmon.rvt",
        "ostrich-{model}/raven-{model}-salmon.rvt.tpl",
        "raven-{model}/raven-{model}-salmon.rvc",
        "raven-{model}/raven-{model}-salmon.rvh",
        "raven-{model}/raven-{model}-salmon.rvi",
        "raven-{model}/raven-{model}-salmon.rvp",
        "raven-{model}/raven-{model}-salmon.rvt",
    ]
    for model in models:
        for entry in entries:
            data_entries.append(entry.format(model=model))

    data_entries.extend(
        [
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
            "raven-hmets/Salmon-River-Near-Prince-George_meteo_daily.rvt",
            "raven-hmets/Salmon-River-Near-Prince-George_Qobs_daily.rvt"
            "raven-mohyse/Salmon-River-Near-Prince-George_Qobs_daily.rvt",
            "raven-mohyse/Salmon-River-Near-Prince-George_meteo_daily.rvt",
            "raven-hbv-ec/Salmon-River-Near-Prince-George_meteo_daily.rvt",
            "raven-hbv-ec/Salmon-River-Near-Prince-George_Qobs_daily.rvt",
            "gr4j_cemaneige/solution.rvc",
            "ostrich-gr4j-cemaneige/OstRandomNumbers.txt",
            "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical%2Brcp85_nex-gddp_2070-2071_subset.nc",
            "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical%2Brcp45_nex-gddp_1971-1972_subset.nc",
            "nrcan/NRCAN_1971-1972_subset.nc",
            "watershed_vector/LSJ_LL.zip",
        ]
    )

    data = dict()
    for filepattern in data_entries:
        if temp_folder is None:
            try:
                data[filepattern] = get_file(
                    filepattern, branch=branch, cache_dir=_local_cache
                )
            except FileNotFoundError:
                continue
        elif temp_folder:
            try:
                data[filepattern] = get_local_testdata(
                    filepattern,
                    temp_folder=temp_folder,
                    branch=branch,
                    _local_cache=_local_cache,
                )
            except FileNotFoundError:
                continue

    if temp_folder is None:
        # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
        _local_cache.joinpath(".data_written").touch()

    return


@pytest.fixture(scope="session")
def threadsafe_data_dir(tmp_path_factory) -> Path:
    """Constructor for worker-session temporary data folders."""
    return Path(tmp_path_factory.getbasetemp().joinpath("data"))


@pytest.fixture(scope="session", autouse=True)
def gather_session_data(threadsafe_data_dir, worker_id):
    """Gather testing data on pytest run.

    When running pytest with multiple workers, one worker will copy data remotely to _default_cache_dir while
    other workers wait using lockfile. Once the lock is released, all workers will copy data to their local
    threadsafe_data_dir."""
    if worker_id == "master":
        return populate_testing_data(branch="add_full_model_name")
    else:
        _default_cache_dir.mkdir(exist_ok=True)
        test_data_being_written = FileLock(_default_cache_dir.joinpath(".lock"))
        with test_data_being_written.acquire() as fl:
            if fl.is_locked:
                populate_testing_data(branch="add_full_model_name")
                fl.release()
            fl.acquire()
            populate_testing_data(
                temp_folder=threadsafe_data_dir, branch="add_full_model_name"
            )


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing file once we are finished.

    This flag prevents remote data from being downloaded multiple times in the same pytest run.
    """

    def remove_data_written_flag():
        flag = _default_cache_dir.joinpath(".data_written")
        if flag.exists():
            flag.unlink()

    request.addfinalizer(remove_data_written_flag)


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
