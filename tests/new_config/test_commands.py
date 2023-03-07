import re
from textwrap import dedent
from typing import Sequence

from pydantic import Field

from ravenpy.new_config import commands as rc
from ravenpy.new_config import options as o
from ravenpy.new_config.base import RV

salmon_land_hru_1 = dict(
    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659, hru_type="land"
)


def test_evaluation_metrics():
    class Test(RV):
        em: Sequence[o.EvaluationMetrics] = Field(
            ["RMSE", "NASH_SUTCLIFFE"], alias="EvaluationMetrics"
        )

    assert Test().to_rv() == ":EvaluationMetrics    RMSE NASH_SUTCLIFFE"


def test_hru():
    hru = rc.HRUs.parse_obj([salmon_land_hru_1])
    hru.to_rv()


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
        Precipitation(algo="PRECIP_RAVEN", source="ATMOS_PRECIP", to="MULTIPLE"),
        SnowTempEvolve(algo="SNOTEMP_NEWTONS", source="SNOW_TEMP"),
    ]

    class Test(RV):
        hydrologic_processes: Sequence[Process] = Field(hp, alias="HydrologicProcesses")

    out = dedent(Test().to_rv())
    assert (
        out
        == dedent(
            """
    :HydrologicProcesses
      :Precipitation        PRECIP_RAVEN        ATMOS_PRECIP        MULTIPLE
      :SnowTempEvolve       SNOTEMP_NEWTONS     SNOW_TEMP
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
    dedent(c.to_rv())


def test_gauge(get_local_testdata):
    from ravenpy.new_config.commands import Gauge

    f = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    g = Gauge.from_nc(f, alt_names={"PRECIP": "rain"})

    s = dedent(g[0].to_rv())
    assert "Data" in s
