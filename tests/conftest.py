import datetime as dt
import os
import shutil
from pathlib import Path
from typing import Optional, Union

import pytest
import xarray as xr
from filelock import FileLock
from xclim.indicators.generic import fit, stats

from ravenpy.utilities.testdata import _default_cache_dir
from ravenpy.utilities.testdata import get_file as _get_file
from ravenpy.utilities.testdata import get_local_testdata as _get_local_testdata

from .common import _convert_2d, _convert_3d

# Additional pytest fixtures for emulators
from .emulators import (  # noqa: F401
    config_rv,
    gr4jcn_config,
    minimal_emulator,
    numeric_config,
    symbolic_config,
)

RAVEN_TESTING_DATA_BRANCH = os.getenv("RAVEN_TESTING_DATA_BRANCH", "master")
SKIP_TEST_DATA = os.getenv("RAVENPY_SKIP_TEST_DATA")
DEFAULT_CACHE = Path(_default_cache_dir)


def populate_testing_data(
    temp_folder: Optional[Path] = None,
    branch: str = RAVEN_TESTING_DATA_BRANCH,
    _local_cache: Path = DEFAULT_CACHE,
) -> None:
    if _local_cache.joinpath(".data_written").exists():
        # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
        return

    models = [
        "gr4j-cemaneige",
        "hbvec",
        "hmets",
        "mohyse",
    ]

    data_entries = list()
    entries = [
        "raven-{model}/Salmon-River-Near-Prince-George_Qobs_daily.rvt",
        "raven-{model}/Salmon-River-Near-Prince-George_meteo_daily.rvt",
        "raven-{model}/raven-{model0}-salmon.rvc",
        "raven-{model}/raven-{model0}-salmon.rvh",
        "raven-{model}/raven-{model0}-salmon.rvi",
        "raven-{model}/raven-{model0}-salmon.rvp",
        "raven-{model}/raven-{model0}-salmon.rvt",
    ]
    for model in models:
        for entry in entries:
            data_entries.append(entry.format(model=model, model0=model.split("-")[0]))

    data_entries.extend(
        [
            "caspar_eccc_hindcasts/geps_watershed.nc",
            "eccc_forecasts/geps_watershed.nc",
            "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp45_nex-gddp_1971-1972_subset.nc",
            "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp85_nex-gddp_2070-2071_subset.nc",
            "famine/famine_input.nc",
            "gr4j_cemaneige/solution.rvc",
            "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
            "nasa/Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff",
            "nrcan/NRCAN_1971-1972_subset.nc",
            "nrcan/NRCAN_2006-2007_subset.nc",
            "polygons/mars.geojson",
            "polygons/mars.zip",
            "polygons/Saskatoon.geojson",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc",
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
            "raven-routing-sample/OTT_sub.zip",
            "raven-routing-sample/VIC_streaminputs.nc",
            "raven-routing-sample/VIC_streaminputs_weights.rvt",
            "raven-routing-sample/VIC_temperatures.nc",
            "raven-routing-sample/VIC_test_nodata.nc",
            "raven-routing-sample/VIC_test_nodata_weights.rvt",
            "raven-routing-sample/WSC02LE024.nc",
            "raven-routing-sample/era5-test-dataset-crop.nc",
            "raven-routing-sample/finalcat_hru_info.zip",
            "raven-routing-sample/lievre_hrus_v21.zip",
            "watershed_vector/LSJ_LL.zip",
            "basinmaker/drainage_region_0175_v2-1/finalcat_info_v2-1.zip",
            "matapedia/Qobs_Matapedia_01BD009.nc",
            "matapedia/Matapedia_meteo_data_2D.nc",
            "matapedia/Matapedia_meteo_data_stations.nc",
        ]
    )

    data = dict()
    for filepattern in data_entries:
        if temp_folder is None:
            try:
                data[filepattern] = _get_file(
                    filepattern, branch=branch, cache_dir=_local_cache
                )
            except FileNotFoundError:
                continue
        elif temp_folder:
            try:
                data[filepattern] = _get_local_testdata(
                    filepattern,
                    temp_folder=temp_folder,
                    branch=branch,
                    _local_cache=_local_cache,
                )
            except FileNotFoundError:
                continue

    return


@pytest.fixture(scope="session")
def threadsafe_data_dir(tmp_path_factory) -> Path:
    """Constructor for worker-session temporary data folders."""
    return Path(tmp_path_factory.getbasetemp().joinpath("data"))


@pytest.fixture(scope="session")
def get_file(threadsafe_data_dir):
    def _get_session_scoped_file(file: Union[str, Path]):
        return _get_file(
            file, cache_dir=threadsafe_data_dir, branch=RAVEN_TESTING_DATA_BRANCH
        )

    return _get_session_scoped_file


@pytest.fixture(scope="session")
def get_local_testdata(threadsafe_data_dir):
    def _get_session_scoped_local_testdata(file: Union[str, Path]):
        return _get_local_testdata(
            file,
            temp_folder=threadsafe_data_dir,
            branch=RAVEN_TESTING_DATA_BRANCH,
            _local_cache=DEFAULT_CACHE,
        )

    return _get_session_scoped_local_testdata


@pytest.fixture(scope="session", autouse=True)
def gather_session_data(threadsafe_data_dir, worker_id):
    """Gather testing data on pytest run.

    When running pytest with multiple workers, one worker will copy data remotely to DEFAULT_CACHE while
    other workers wait using lockfile. Once the lock is released, all workers will copy data to their local
    threadsafe_data_dir."""
    if worker_id == "master":
        if not SKIP_TEST_DATA:
            populate_testing_data(branch=RAVEN_TESTING_DATA_BRANCH)
    else:
        if not SKIP_TEST_DATA:
            DEFAULT_CACHE.mkdir(exist_ok=True)
            test_data_being_written = FileLock(DEFAULT_CACHE.joinpath(".lock"))
            with test_data_being_written as fl:
                # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
                populate_testing_data(branch=RAVEN_TESTING_DATA_BRANCH)
                DEFAULT_CACHE.joinpath(".data_written").touch()
            fl.acquire()
        shutil.copytree(DEFAULT_CACHE, threadsafe_data_dir)


@pytest.fixture(scope="session")
def q_sim_1(threadsafe_data_dir, get_local_testdata):
    """A file storing a Raven streamflow simulation over one basin."""
    return get_local_testdata(
        "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
    )


@pytest.fixture(scope="session")
def ts_stats(q_sim_1, tmp_path):
    with xr.open_dataset(q_sim_1) as ds:
        q = ds.q_sim
        ts = stats(q, op="max")
        fn = tmp_path / "ts_stats.nc"
        ts.to_netcdf_(fn)
    return fn


@pytest.fixture(scope="session")
def fit_params(ts_stats, threadsafe_data_dir):
    fn = threadsafe_data_dir / "fit.nc"

    if not fn.exists():
        with xr.open_dataset(ts_stats) as ds:
            name = list(ds.data_vars.keys()).pop()
            q = ds[name]
            p = fit(q, dist="gumbel_r")
            p.to_netcdf(fn)

    yield fn


@pytest.fixture(scope="session")
def bad_netcdf(get_local_testdata, threadsafe_data_dir):
    fn = threadsafe_data_dir / "bad_netcdf.nc"

    salmon_file = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )

    if not fn.exists():
        # Test with scalar elevation. Should normally have a station dimension, but not always the case.
        with xr.open_dataset(salmon_file) as ds:
            ds["station_name"] = ds["station_name"].astype("str")
            ds["elevation"] = 1.0
            ds.to_netcdf(fn)

    return fn


@pytest.fixture(scope="session")
def salmon_hru():
    out = {
        "land": dict(
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            hru_type="land",
        )
    }
    yield out


# Used in test_emulators.py
@pytest.fixture(scope="session")
def input2d(salmon, threadsafe_data_dir):
    """Convert 1D input to 2D output by copying all the time series along a new region dimension."""
    fn_out = threadsafe_data_dir / "input2d.nc"

    if not fn_out.exists():
        _convert_2d(salmon).to_netcdf(fn_out)

    return fn_out


@pytest.fixture(scope="session")
def input3d(salmon, threadsafe_data_dir):
    """Convert 1D input to 2D output by copying all the time series along a new region dimension."""
    fn_out = threadsafe_data_dir / "input3d.nc"

    if not fn_out.exists():
        ds = _convert_3d(salmon)
        ds = ds.drop_vars("qobs")
        ds.to_netcdf(fn_out)
        ds.close()

    return fn_out


@pytest.fixture
def dummy_config():
    """Return an almost empty config class and the parameter dataclass."""
    from pydantic import Field
    from pydantic.dataclasses import dataclass

    from ravenpy.config import options as o
    from ravenpy.config.base import Params, Sym, SymConfig, Variable
    from ravenpy.config.rvs import Config

    @dataclass(config=SymConfig)
    class P(Params):
        X1: Sym = Variable("X1")

    class TestConfig(Config):
        params: P = P()
        calendar: o.Calendar = Field("JULIAN", alias="Calendar")
        air_snow_coeff: Optional[Sym] = Field(1 - P.X1, alias="AirSnowCoeff")

    return TestConfig, P


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing file once we are finished.

    This flag prevents remote data from being downloaded multiple times in the same pytest run.
    """

    def remove_data_written_flag():
        flag = DEFAULT_CACHE.joinpath(".data_written")
        if flag.exists():
            flag.unlink()

    request.addfinalizer(remove_data_written_flag)


if __name__ == "__main__":
    populate_testing_data(branch=RAVEN_TESTING_DATA_BRANCH)
