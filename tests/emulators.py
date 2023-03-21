import datetime as dt

import pytest
import xarray as xr

from ravenpy.new_config import commands as rc
from ravenpy.new_config.emulators import GR4JCN, HBVEC, HMETS, Mohyse

alt_names = {
    "RAINFALL": "rain",
    "TEMP_MIN": "tmin",
    "TEMP_MAX": "tmax",
    "PET": "pet",
    "HYDROGRAPH": "qobs",
    "SNOWFALL": "snow",
}


names = ["GR4JCN", "HMETS", "Mohyse", "HBVEC"][:]
configs = {"GR4JCN": GR4JCN, "HMETS": HMETS, "Mohyse": Mohyse, "HBVEC": HBVEC}
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
}


@pytest.fixture
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


@pytest.fixture(scope="session", params=names)
def symbolic_config(salmon_meteo, salmon_hru, request):
    """Config with symbolic parameters"""
    name = request.param
    cls = configs[name]
    data_type = ["RAINFALL", "TEMP_MIN", "TEMP_MAX", "SNOWFALL"]

    # Extra attributes for gauges
    gextras = {1: {"elevation": salmon_hru["land"]["elevation"]}}

    if name == "HBVEC":
        ds = xr.open_dataset(salmon_meteo)
        gextras[1]["monthly_ave_temperature"] = (
            ((ds.tmin + ds.tmax) / 2).groupby("time.month").mean().values.tolist()
        )
        gextras[1]["monthly_ave_evaporation"] = (
            ds.pet.groupby("time.month").mean().values.tolist()
        )

    # Extra attributes for emulator
    extras = {}

    if name in ["GR4JCN"]:
        extras.update(
            dict(
                GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
            )
        )
    if name in ["HMETS"]:
        data_type.extend(["SNOWFALL", "PET"])

    yield cls(
        Gauge=rc.Gauge.from_nc(
            salmon_meteo,
            data_type=data_type,
            alt_names=alt_names,
            extra=gextras,
        ),
        ObservationData=rc.ObservationData.from_nc(salmon_meteo, alt_names="qobs"),
        HRUs=[salmon_hru["land"]],
        StartDate=dt.datetime(2000, 1, 1),
        EndDate=dt.datetime(2002, 1, 1),
        RunName="test",
        EvaluationMetrics=("NASH_SUTCLIFFE",),
        **extras,
    )


@pytest.fixture(scope="session", params=names)
def numeric_config(symbolic_config, request):
    name = request.param
    yield symbolic_config.set_params(params[name])


@pytest.fixture(scope="session")
def config_rv(tmp_path_factory, numeric_config):
    name = numeric_config.__class__.__name__
    out = tmp_path_factory.mktemp(name) / "config"
    numeric_config.build(out)
    yield out
