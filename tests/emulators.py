import datetime as dt

import pytest
import xarray as xr

from ravenpy.new_config import commands as rc
from ravenpy.new_config.emulators import (
    GR4JCN,
    HBVEC,
    HMETS,
    HYPR,
    SACSMA,
    CanadianShield,
    Mohyse,
)

alt_names = {
    "RAINFALL": "rain",
    "TEMP_MIN": "tmin",
    "TEMP_MAX": "tmax",
    "PET": "pet",
    "HYDROGRAPH": "qobs",
    "SNOWFALL": "snow",
}


names = ["GR4JCN", "HMETS", "Mohyse", "HBVEC", "CanadianShield", "HYPR", "SACSMA"][-1:]
configs = {
    "GR4JCN": GR4JCN,
    "HMETS": HMETS,
    "Mohyse": Mohyse,
    "HBVEC": HBVEC,
    "CanadianShield": CanadianShield,
    "HYPR": HYPR,
    "SACSMA": SACSMA,
}
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
}


@pytest.fixture(scope="session", params=names)
def symbolic_config(salmon_meteo, salmon_hru, request):
    """Config with symbolic parameters"""
    name = request.param
    cls = configs[name]
    data_type = ["RAINFALL", "TEMP_MIN", "TEMP_MAX", "SNOWFALL"]

    # Extra attributes for gauges
    gextras = {1: {"elevation": salmon_hru["land"]["elevation"]}}

    if name in ["HBVEC", "HYPR"]:
        ds = xr.open_dataset(salmon_meteo)
        gextras[1]["monthly_ave_temperature"] = (
            ((ds.tmin + ds.tmax) / 2).groupby("time.month").mean().values.tolist()
        )
        gextras[1]["monthly_ave_evaporation"] = (
            ds.pet.groupby("time.month").mean().values.tolist()
        )

    # Extra attributes for emulator
    extras = {}
    hrus = [salmon_hru["land"]]
    if name in ["GR4JCN"]:
        extras.update(
            dict(
                GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
            )
        )
    if name in ["HMETS"]:
        data_type.extend(["SNOWFALL", "PET"])

    if name in ["CanadianShield"]:
        hrus = [salmon_hru["land"], salmon_hru["land"]]
        extras["SuppressOutput"] = True

    yield cls(
        Gauge=rc.Gauge.from_nc(
            salmon_meteo,
            data_type=data_type,
            alt_names=alt_names,
            extra=gextras,
        ),
        ObservationData=rc.ObservationData.from_nc(salmon_meteo, alt_names="qobs"),
        HRUs=hrus,
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