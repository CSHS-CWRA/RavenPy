import datetime as dt
import os
import shutil
from pathlib import Path
from typing import Optional, Union

import pytest
import xarray as xr
from filelock import FileLock
from xclim.indicators.generic import fit, stats

from ravenpy.config import commands as rc
from ravenpy.config.emulators import (
    GR4JCN,
    HBVEC,
    HMETS,
    HYPR,
    SACSMA,
    Blended,
    CanadianShield,
    Mohyse,
)
from ravenpy.utilities.testdata import _default_cache_dir
from ravenpy.utilities.testdata import get_file as _get_file
from ravenpy.utilities.testdata import get_local_testdata as _get_local_testdata

from .common import _convert_2d, _convert_3d

TESTDATA_BRANCH = os.getenv("RAVENPY_TESTDATA_BRANCH", "master")
SKIP_TEST_DATA = os.getenv("RAVENPY_SKIP_TEST_DATA")


def populate_testing_data(
    temp_folder: Optional[Path] = None,
    branch: str = TESTDATA_BRANCH,
    _local_cache: Path = _default_cache_dir,
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
        return _get_file(file, cache_dir=threadsafe_data_dir, branch=TESTDATA_BRANCH)

    return _get_session_scoped_file


@pytest.fixture(scope="session")
def get_local_testdata(threadsafe_data_dir):
    def _get_session_scoped_local_testdata(file: Union[str, Path]):
        return _get_local_testdata(
            file,
            temp_folder=threadsafe_data_dir,
            branch=TESTDATA_BRANCH,
            _local_cache=_default_cache_dir,
        )

    return _get_session_scoped_local_testdata


@pytest.fixture(scope="session", autouse=True)
def gather_session_data(threadsafe_data_dir, worker_id):
    """Gather testing data on pytest run.

    When running pytest with multiple workers, one worker will copy data remotely to _default_cache_dir while
    other workers wait using lockfile. Once the lock is released, all workers will copy data to their local
    threadsafe_data_dir."""
    if worker_id == "master":
        if not SKIP_TEST_DATA:
            populate_testing_data(branch=TESTDATA_BRANCH)
    else:
        if not SKIP_TEST_DATA:
            _default_cache_dir.mkdir(exist_ok=True)
            test_data_being_written = FileLock(_default_cache_dir.joinpath(".lock"))
            with test_data_being_written as fl:
                # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
                populate_testing_data(branch=TESTDATA_BRANCH)
                _default_cache_dir.joinpath(".data_written").touch()
            fl.acquire()
        shutil.copytree(_default_cache_dir, threadsafe_data_dir)


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
    return _get_file(
        "hydro_simulations/raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc",
        cache_dir=threadsafe_data_dir,
        branch=TESTDATA_BRANCH,
    )


@pytest.fixture(scope="session")
def ts_stats(q_sim_1, tmp_path):
    q = xr.open_dataset(q_sim_1).q_sim
    ts = stats(q, op="max")
    fn = tmp_path / "ts_stats.nc"
    ts.to_netcdf_(fn)
    return fn


@pytest.fixture(scope="session")
def fit_params(ts_stats, tmp_path):
    fn = tmp_path / "fit.nc"

    if not fn.exists():
        ds = xr.open_dataset(ts_stats)
        name = list(ds.data_vars.keys()).pop()
        q = ds[name]
        p = fit(q, dist="gumbel_r")
        p.to_netcdf(fn)

    return fn


@pytest.fixture(scope="session")
def salmon_hru():
    out = {}
    out["land"] = dict(
        area=4250.6,
        elevation=843.0,
        latitude=54.4848,
        longitude=-123.3659,
        hru_type="land",
    )
    return out


@pytest.fixture(scope="session")
def salmon_meteo(get_local_testdata):
    return get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )


@pytest.fixture(scope="session")
def salmon(threadsafe_data_dir):
    """A file storing a Raven streamflow simulation over one basin."""
    return _get_file(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
        cache_dir=threadsafe_data_dir,
        branch=TESTDATA_BRANCH,
    )


# Used in test_emulators.py
@pytest.fixture(scope="session")
def input2d(threadsafe_data_dir, salmon):
    """Convert 1D input to 2D output by copying all the time series along a new region dimension."""

    fn_out = threadsafe_data_dir / "input2d.nc"
    if not fn_out.exists():
        _convert_2d(salmon).to_netcdf(fn_out)

    return fn_out


# Used in test_emulators.py
@pytest.fixture(scope="session")
def input3d(threadsafe_data_dir, salmon):
    """Convert 1D input to 2D output by copying all the time series along a new region dimension."""

    fn_out = threadsafe_data_dir / "input3d.nc"
    if not fn_out.exists():
        ds = _convert_3d(salmon)
        ds = ds.drop_vars("qobs")
        ds.to_netcdf(fn_out)

    return fn_out


# Alternative names for variables in meteo forcing file
alt_names = {
    "RAINFALL": "rain",
    "TEMP_MIN": "tmin",
    "TEMP_MAX": "tmax",
    "PET": "pet",
    "HYDROGRAPH": "qobs",
    "SNOWFALL": "snow",
}

# Names of emulators to be included in tests
names = [
    "GR4JCN",
    "HMETS",
    "Mohyse",
    "HBVEC",
    "CanadianShield",
    "HYPR",
    "SACSMA",
    "Blended",
]

# Mapping of emulator names to configuration classes
configs = {
    "GR4JCN": GR4JCN,
    "HMETS": HMETS,
    "Mohyse": Mohyse,
    "HBVEC": HBVEC,
    "CanadianShield": CanadianShield,
    "HYPR": HYPR,
    "SACSMA": SACSMA,
    "Blended": Blended,
}

# Model parameters that match NSE criteria computed by Juliane
params = {
    "GR4JCN": [0.529, -3.396, 407.29, 1.072, 16.9, 0.947],
    "HMETS": [
        9.5019,
        0.2774,
        6.3942,
        0.6884,
        1.2875,
        5.4134,
        2.3641,
        0.0973,
        0.0464,
        0.1998,
        0.0222,
        -1.0919,
        2.6851,
        0.3740,
        1.0000,
        0.4739,
        0.0114,
        0.0243,
        0.0069,
        310.7211,
        916.1947,
    ],
    "Mohyse": (
        1.0,
        0.0468,
        4.2952,
        2.658,
        0.4038,
        0.0621,
        0.0273,
        0.0453,
        0.9039,
        5.6167,
    ),
    "HBVEC": (
        0.05984519,
        4.072232,
        2.001574,
        0.03473693,
        0.09985144,
        0.506052,
        3.438486,
        38.32455,
        0.4606565,
        0.06303738,
        2.277781,
        4.873686,
        0.5718813,
        0.04505643,
        0.877607,
        18.94145,
        2.036937,
        0.4452843,
        0.6771759,
        1.141608,
        1.024278,
    ),
    "CanadianShield": (
        4.72304300e-01,
        8.16392200e-01,
        9.86197600e-02,
        3.92699900e-03,
        4.69073600e-02,
        4.95528400e-01,
        6.803492000e00,
        4.33050200e-03,
        1.01425900e-05,
        1.823470000e00,
        5.12215400e-01,
        9.017555000e00,
        3.077103000e01,
        5.094095000e01,
        1.69422700e-01,
        8.23412200e-02,
        2.34595300e-01,
        7.30904000e-02,
        1.284052000e00,
        3.653415000e00,
        2.306515000e01,
        2.402183000e00,
        2.522095000e00,
        5.80344900e-01,
        1.614157000e00,
        6.031781000e00,
        3.11129800e-01,
        6.71695100e-02,
        5.83759500e-05,
        9.824723000e00,
        9.00747600e-01,
        8.04057300e-01,
        1.179003000e00,
        7.98001300e-01,
    ),
    "HYPR": (
        -1.856410e-01,
        2.92301100e00,
        3.1194200e-02,
        4.3982810e-01,
        4.6509760e-01,
        1.1770040e-01,
        1.31236800e01,
        4.0417950e-01,
        1.21225800e00,
        5.91273900e01,
        1.6612030e-01,
        4.10501500e00,
        8.2296110e-01,
        4.15635200e01,
        5.85111700e00,
        6.9090140e-01,
        9.2459950e-01,
        1.64358800e00,
        1.59920500e00,
        2.51938100e00,
        1.14820100e00,
    ),
    "SACSMA": (
        0.0100000,
        0.0500000,
        0.3000000,
        0.0500000,
        0.0500000,
        0.1300000,
        0.0250000,
        0.0600000,
        0.0600000,
        1.0000000,
        40.000000,
        0.0000000,
        0.0000000,
        0.1000000,
        0.0000000,
        0.0100000,
        1.5000000,
        0.4827523,
        4.0998200,
        1.0000000,
        1.0000000,
    ),
    "Blended": (
        2.930702e-02,
        2.211166e00,
        2.166229e00,
        0.0002254976,
        2.173976e01,
        1.565091e00,
        6.211146e00,
        9.313578e-01,
        3.486263e-02,
        0.251835,
        0.0002279250,
        1.214339e00,
        4.736668e-02,
        0.2070342,
        7.806324e-02,
        -1.336429e00,
        2.189741e-01,
        3.845617e00,
        2.950022e-01,
        4.827523e-01,
        4.099820e00,
        1.283144e01,
        5.937894e-01,
        1.651588e00,
        1.705806,
        3.719308e-01,
        7.121015e-02,
        1.906440e-02,
        4.080660e-01,
        9.415693e-01,
        -1.856108e00,
        2.356995e00,
        1.0e00,
        1.0e00,
        7.510967e-03,
        5.321608e-01,
        2.891977e-02,
        9.605330e-01,
        6.128669e-01,
        9.558293e-01,
        1.008196e-01,
        9.275730e-02,
        7.469583e-01,
    ),
}


@pytest.fixture(scope="session")
def gr4jcn_config(salmon_meteo, salmon_hru) -> (GR4JCN, params):
    """Return symbolic config and params for basic gr4jcn."""

    yield GR4JCN(
        Gauge=[
            rc.Gauge.from_nc(
                salmon_meteo,
                data_type=["RAINFALL", "TEMP_MIN", "TEMP_MAX", "SNOWFALL"],
                alt_names=alt_names,
                data_kwds={"ALL": {"elevation": salmon_hru["land"]["elevation"]}},
            ),
        ],
        ObservationData=[rc.ObservationData.from_nc(salmon_meteo, alt_names="qobs")],
        HRUs=[salmon_hru["land"]],
        StartDate=dt.datetime(2000, 1, 1),
        Duration=15,
        EvaluationMetrics=("NASH_SUTCLIFFE",),
    ), params["GR4JCN"]


@pytest.fixture(scope="session", params=names)
def symbolic_config(salmon_meteo, salmon_hru, request):
    """Emulator configuration instantiated with symbolic parameters."""
    name = request.param
    cls = configs[name]
    data_type = ["RAINFALL", "TEMP_MIN", "TEMP_MAX", "SNOWFALL"]

    # Extra attributes for gauges
    gextras = {"ALL": {"elevation": salmon_hru["land"]["elevation"]}}

    if name in ["HBVEC", "HYPR"]:
        with xr.open_dataset(salmon_meteo) as ds:
            gextras["ALL"]["monthly_ave_temperature"] = (
                ((ds.tmin + ds.tmax) / 2).groupby("time.month").mean().values.tolist()
            )
            gextras["ALL"]["monthly_ave_evaporation"] = (
                ds.pet.groupby("time.month").mean().values.tolist()
            )

    # Extra attributes for emulator
    extras = {}

    if name in ["CanadianShield", "HYPR", "SACSMA", "Blended"]:
        salmon_hru["slope"] = 0.01234

    hrus = [salmon_hru["land"]]

    if name in ["GR4JCN"]:
        extras.update(
            dict(
                GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
            )
        )
    if name in ["HMETS"]:
        data_type.extend(["PET"])

    if name in ["CanadianShield"]:
        hrus = [salmon_hru["land"], salmon_hru["land"]]
        extras["SuppressOutput"] = True

    yield name, cls(
        Gauge=[
            rc.Gauge.from_nc(
                salmon_meteo,
                data_type=data_type,
                alt_names=alt_names,
                data_kwds=gextras,
            ),
        ],
        ObservationData=[rc.ObservationData.from_nc(salmon_meteo, alt_names="qobs")],
        HRUs=hrus,
        StartDate=dt.datetime(2000, 1, 1),
        EndDate=dt.datetime(2002, 1, 1),
        RunName="test",
        EvaluationMetrics=("NASH_SUTCLIFFE",),
        **extras,
    )


@pytest.fixture(scope="session")
def numeric_config(symbolic_config):
    """Emulator configuration instantiated with numeric parameters."""
    name, conf = symbolic_config
    yield name, conf.set_params(params[name])


@pytest.fixture(scope="session")
def minimal_emulator(salmon_meteo, salmon_hru):
    """Return the config for a single emulator."""
    cls = configs["HMETS"]
    data_type = ["RAINFALL", "TEMP_MIN", "TEMP_MAX", "SNOWFALL", "PET"]

    return cls(
        params=params["HMETS"],
        Gauge=[
            rc.Gauge.from_nc(
                salmon_meteo,
                data_type=data_type,
                alt_names=alt_names,
            ),
        ],
        ObservationData=[rc.ObservationData.from_nc(salmon_meteo, alt_names="qobs")],
        HRUs=[salmon_hru["land"]],
        StartDate=dt.datetime(2000, 1, 1),
        Duration=15,
    )


@pytest.fixture(scope="session")
def config_rv(tmp_path_factory, numeric_config):
    """Directory with emulator configuration files written to disk."""
    name, conf = numeric_config
    out = tmp_path_factory.mktemp(name) / "config"
    conf.write_rv(out)
    yield name, out


@pytest.fixture
def dummy_config():
    """Return a almost empty config class and the parameter dataclass."""
    from pydantic import Field
    from pydantic.dataclasses import dataclass

    from ravenpy.config import options as o
    from ravenpy.config.base import Sym, SymConfig, Variable
    from ravenpy.config.rvs import Config

    @dataclass(config=SymConfig)
    class P:
        X1: Sym = Variable("X1")

    class TestConfig(Config):
        params: P = P()
        calendar: o.Calendar = Field("JULIAN", alias="Calendar")
        air_snow_coeff: Sym = Field(1 - P.X1, alias="AirSnowCoeff")

    return TestConfig, P


#
if __name__ == "__main__":
    populate_testing_data(branch=TESTDATA_BRANCH)
