import datetime as dt
import re
from textwrap import dedent
from typing import Any, Sequence, Type

from pydantic import Field
from pydantic.dataclasses import dataclass

from ravenpy.config import options as o
from ravenpy.new_config import commands as rc
from ravenpy.new_config.base import RV, Command, Sym, SymConfig, Variable, encoder
from ravenpy.new_config.emulators.gr4jcn import GR4JCN
from ravenpy.new_config.rvs import Config

alt_names = {
    "PRECIP": "rain",
    "TEMP_MIN": "tmin",
    "TEMP_MAX": "tmax",
    "PET": "pet",
    "HYDROGRAPH": "qobs",
}

salmon_land_hru_1 = dict(
    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659, hru_type="land"
)


@dataclass(config=SymConfig)
class P:
    X1: Sym = Variable("X1")


def test_evaluation_metrics():
    class Test(RV):
        em: Sequence[o.EvaluationMetrics] = Field(
            ["RMSE", "NASH_SUTCLIFFE"], alias="EvaluationMetrics"
        )

    assert Test().to_rv() == ":EvaluationMetrics    RMSE NASH_SUTCLIFFE"


def test_emulator():
    class TestConfig(Config):
        params: P
        calendar: o.Calendar = Field("JULIAN", alias="Calendar")
        air_snow_coeff: Sym = Field(1 - P.X1, alias="AirSnowCoeff")

    t = TestConfig(params=[0.5], Calendar="NOLEAP")
    assert t.air_snow_coeff == 0.5


def test_gr4j(tmpdir, get_local_testdata):
    f = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )

    m = GR4JCN(
        params=[0.529, -3.396, 407.29, 1.072, 16.9, 0.947],
        Gauge=rc.Gauge.from_nc(
            f,
            data_type=["PRECIP", "TEMP_MIN", "TEMP_MAX"],
            alt_names=alt_names,
            extra={1: {"elevation": salmon_land_hru_1["elevation"]}},
        ),
        ObservationData=rc.ObservationData.from_nc(f, alt_names="qobs"),
        HRUs=[salmon_land_hru_1],
        StartDate=dt.datetime(2000, 1, 1),
        EndDate=dt.datetime(2002, 1, 1),
        RunName="test",
        CustomOutput=rc.CustomOutput("YEARLY", "AVERAGE", "PRECIP", "ENTIRE_WATERSHED"),
        GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
    )
    m.write("/tmp/test_config")


def test_hru():
    hru = rc.HRUs.parse_obj([salmon_land_hru_1])
    hru.to_rv()


def test_soil_classes():
    from ravenpy.new_config.commands import SoilClasses

    c = SoilClasses(names=["TOPSOIL", "FAST_RES", "SLOW_RES"])
    assert dedent(c.to_rv()) == dedent(
        """
    :SoilClasses
      TOPSOIL
      FAST_RES
      SLOW_RES
    :EndSoilClasses
    """
    )


def test_vegetation_classes():
    from ravenpy.new_config.commands import VegetationClasses

    class Test(RV):
        vegetation_classes: VegetationClasses = Field(
            VegetationClasses.parse_obj([{"name": "VEG_ALL"}, {"name": "VEG_WATER"}]),
            alias="VegetationClasses",
        )

    s = dedent(Test().to_rv())
    assert (
        s
        == dedent(
            """
    :VegetationClasses
      :Attributes     ,        MAX_HT,       MAX_LAI, MAX_LEAF_COND
      :Units          ,             m,          none,      mm_per_s
      VEG_ALL         ,           0.0,           0.0,           0.0
      VEG_WATER       ,           0.0,           0.0,           0.0
    :EndVegetationClasses"""
        ).strip()
    )


def test_soil_profiles(x=42.0):
    from ravenpy.new_config.commands import SoilProfile

    c = SoilProfile.parse_obj(
        {
            "name": "DEFAULT_P",
            "soil_classes": ["TOPSOIL", "FAST_RES", "SLOW_RES"],
            "thicknesses": [x, 100.0, 100.0],
        }
    )

    pat = r"DEFAULT_P\s*,\s*3,\s*TOPSOIL,\s*(.+),\s*FAST_RES,\s*100.0,\s*SLOW_RES,\s*100.0"

    m = re.findall(pat, c.to_rv())[0]
    assert m == "42.0"


def test_hydrologic_processes():
    from ravenpy.new_config.commands import Process
    from ravenpy.new_config.processes import Precipitation, SnowTempEvolve

    hp = [
        Precipitation(algo="PRECIP_RAVEN", source=["ATMOS_PRECIP", "MULTIPLE"]),
        SnowTempEvolve(algo="SNOTEMP_NEWTONS", source=["SNOW_TEMP"]),
    ]

    class Test(RV):
        hydrologic_processes: Sequence[Process] = Field(hp, alias="HydrologicProcesses")

    out = dedent(Test().to_rv())
    assert (
        out
        == dedent(
            """
    :HydrologicProcesses
      :Precipitation ATMOS_PRECIP MULTIPLE
      :SnowTempEvolve SNOW_TEMP
    :EndHydrologicProcesses
    """
        ).strip()
    )


def test_read_from_netcdf(get_local_testdata):
    from ravenpy.new_config.commands import ReadFromNetCDF

    f = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    c = ReadFromNetCDF.from_nc(f, "PRECIP", station_idx=1, alt_names=("rain",))
    s = dedent(c.to_rv())

    pat = re.compile(
        r"""
    :ReadFromNetCDF\s+
      :FileNameNC\s+.+/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc\s+
      :VarNameNC\s+rain\s+
      :DimNamesNC\s+time\s+
      :StationIdx\s+1\s+
      :LatitudeVarNameNC\s+lat\s+
      :LongitudeVarNameNC\s+lon\s+
    :EndReadFromNetCDF
    """,
        re.VERBOSE + re.MULTILINE,
    )
    assert pat.search(s) is not None


def test_station_forcing(get_local_testdata):
    from ravenpy.new_config.commands import StationForcing

    f = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc"
    )
    c = StationForcing.from_nc(f, "PRECIP", station_idx=1, alt_names="rain")
    s = dedent(c.to_rv())


def test_gauge(get_local_testdata):
    from ravenpy.new_config.commands import Gauge

    f = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    g = Gauge.from_nc(f, alt_names={"PRECIP": "rain"})

    s = dedent(g[0].to_rv())
    assert "Data" in s


def a_test_defaults():
    import pint

    ureg = pint.UnitRegistry()

    for name, u in defaults.units.items():
        ureg(u)
