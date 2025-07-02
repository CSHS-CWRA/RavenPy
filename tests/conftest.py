from __future__ import annotations

import logging
import os
from pathlib import Path

import pytest
import xarray as xr

# Additional pytest fixtures for emulators
from emulators import (  # noqa: F401
    config_rv,
    gr4jcn_config,
    minimal_emulator,
    numeric_config,
    symbolic_config,
)
from xclim.indicators.generic import fit, stats

from ravenpy.testing.utils import (
    TESTDATA_BRANCH,
    TESTDATA_CACHE_DIR,
    TESTDATA_REPO_URL,
    default_testdata_cache,
    gather_testing_data,
)
from ravenpy.testing.utils import open_dataset as _open_dataset
from ravenpy.testing.utils import (
    testing_setup_warnings,
)
from ravenpy.testing.utils import yangtze as _yangtze


def convert_2d(fn):
    """Take the 1D Salmon time series and convert it to a 2D time series.

    Example
    -------
    >>> fn = "./testdata/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    >>> fn2 = "./testdata/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc"
    >>> _convert_2d(fn).to_netcdf(fn2, "w")
    """
    features = {
        "name": "Salmon",
        "area": 4250.6,
        "elevation": 843.0,
        "latitude": 54.4848,
        "longitude": -123.3659,
    }
    ds = xr.open_dataset(fn, decode_times=False).rename({"nstations": "region"})

    out = xr.Dataset(
        coords={
            "lon": ds.lon.expand_dims("lon").squeeze("region"),
            "lat": ds.lat.expand_dims("lat").squeeze("region"),
            "time": ds.time,
        }
    )

    for v in ds.data_vars:
        if v not in ["lon", "lat"]:
            out[v] = ds[v].expand_dims("region", axis=1)

    # Add geometry feature variables
    for key, val in features.items():
        out[key] = xr.DataArray(name=key, data=[val], dims="region")

    return out


def convert_3d(fn):
    """Take the 1D Salmon time series and convert it to a 3D time series.

    Example
    -------
    >>> fn = "./testdata/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    >>> fn3 = "./testdata/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc"
    >>> _convert_3d(fn).to_netcdf(fn3, "w")
    """
    elevation = [[843.0]]
    ds = xr.open_dataset(fn, decode_times=False)

    out = xr.Dataset(
        coords={
            "lon": ds.lon.expand_dims("lon").squeeze("nstations"),
            "lat": ds.lat.expand_dims("lat").squeeze("nstations"),
            "time": ds.time,
        }
    )

    for v in ds.data_vars:
        if v not in ["lon", "lat", "time"]:
            out[v] = ds[v]
            out[v] = out[v].expand_dims(
                ["lon", "lat"]
            )  # Needs to be in other step to keep attributes

    out["elevation"] = xr.DataArray(
        data=elevation,
        dims=["lon", "lat"],
        attrs={"units": "m", "standard_name": "altitude"},
    )

    return out


@pytest.fixture(scope="session")
def threadsafe_data_dir(tmp_path_factory) -> Path:
    return Path(tmp_path_factory.getbasetemp().joinpath("data"))


@pytest.fixture(scope="session")
def yangtze(threadsafe_data_dir, worker_id):
    return _yangtze(
        repo=TESTDATA_REPO_URL,
        branch=TESTDATA_BRANCH,
        cache_dir=(
            TESTDATA_CACHE_DIR if worker_id == "master" else threadsafe_data_dir
        ),
    )


@pytest.fixture(scope="session")
def open_dataset(threadsafe_data_dir, worker_id):
    def _open_session_scoped_file(file: str | os.PathLike, **xr_kwargs):
        _yangtze_kwargs = {
            "branch": TESTDATA_BRANCH,
            "repo": TESTDATA_REPO_URL,
            "cache_dir": (
                TESTDATA_CACHE_DIR if worker_id == "master" else threadsafe_data_dir
            ),
        }
        xr_kwargs.setdefault("cache", True)
        xr_kwargs.setdefault("engine", "h5netcdf")
        return _open_dataset(
            file,
            _yangtze_kwargs=_yangtze_kwargs,
            **xr_kwargs,
        )

    return _open_session_scoped_file


@pytest.fixture(scope="session")
def q_sim_1(yangtze):
    """A file storing a Raven streamflow simulation over one basin."""
    return yangtze.fetch(
        "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
    )


@pytest.fixture(scope="session")
def ts_stats(q_sim_1, threadsafe_data_dir):
    fn = threadsafe_data_dir / "ts_stats.nc"

    if not fn.exists():
        with xr.open_dataset(q_sim_1) as ds:
            q = ds.q_sim
            ts = stats(q, op="max")
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
    return fn


@pytest.fixture
def bad_netcdf(yangtze, tmp_path):
    fn = tmp_path / "bad_netcdf.nc"

    salmon_file = yangtze.fetch(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )

    # Test with scalar elevation. Should normally have a station dimension, but not always the case.
    with xr.open_dataset(salmon_file) as ds:
        ds["station_name"] = ds["station_name"].astype("str")
        ds["elevation"] = 1.0
        ds.to_netcdf(fn)
    return fn


@pytest.fixture(scope="session")
def salmon_hru():
    out = {
        "land": {
            "area": 4250.6,
            "elevation": 843.0,
            "latitude": 54.4848,
            "longitude": -123.3659,
            "hru_type": "land",
        }
    }
    return out


# Used in test_emulators.py
@pytest.fixture(scope="session")
def input2d(salmon, threadsafe_data_dir):
    """Convert 1D input to 2D output by copying all the time series along a new region dimension."""
    fn_out = threadsafe_data_dir / "input2d.nc"

    if not fn_out.exists():
        convert_2d(salmon).to_netcdf(fn_out)

    return fn_out


@pytest.fixture(scope="session")
def input3d(salmon, threadsafe_data_dir):
    """Convert 1D input to 2D output by copying all the time series along a new region dimension."""
    fn_out = threadsafe_data_dir / "input3d.nc"

    if not fn_out.exists():
        ds = convert_3d(salmon)
        ds = ds.drop_vars("qobs")
        ds.to_netcdf(fn_out)
        ds.close()

    return fn_out


@pytest.fixture
def dummy_config():
    """Return an almost empty config class and the parameter dataclass."""
    from pydantic import Field
    from pydantic.dataclasses import dataclass, rebuild_dataclass

    from ravenpy.config import options as o
    from ravenpy.config.base import Params, Sym, SymConfig, Variable
    from ravenpy.config.rvs import Config

    @dataclass(config=SymConfig)
    class P(Params):
        X1: Sym = Variable("X1")

    rebuild_dataclass(P)

    class TestConfig(Config):
        params: P = P()
        calendar: o.Calendar = Field("JULIAN", alias="Calendar")
        air_snow_coeff: Sym = Field(1 - P.X1, alias="AirSnowCoeff")

    return TestConfig, P


@pytest.fixture(autouse=True, scope="session")
def gather_session_data(request, yangtze, worker_id):
    """
    Gather testing data on pytest run.

    When running pytest with multiple workers, one worker will copy data remotely to the default cache dir while
    other workers wait using a lockfile. Once the lock is released, all workers will then copy data to their local
    threadsafe_data_dir. As this fixture is scoped to the session, it will only run once per pytest run.
    """
    testing_setup_warnings()
    gather_testing_data(worker_cache_dir=yangtze.path, worker_id=worker_id)

    def remove_data_written_flag():
        """Clean up the cache folder once we are finished."""
        flag = default_testdata_cache.joinpath(".data_written")
        if flag.exists():
            try:
                flag.unlink()
            except FileNotFoundError:
                logging.info(
                    "Teardown race condition occurred: .data_written flag already removed. Lucky!"
                )
                pass

    request.addfinalizer(remove_data_written_flag)
