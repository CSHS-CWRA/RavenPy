import re
from textwrap import dedent
from typing import Any, Sequence, Type

from pydantic import Field
from pydantic.dataclasses import dataclass

from ravenpy.config import options as o
from ravenpy.new_config.base import Command, Sym, SymConfig, Variable, encoder
from ravenpy.new_config.rvs import Config


@dataclass(config=SymConfig)
class P:
    X1: Sym = Variable("X1")


def test_emulator():
    class TestConfig(Config):
        params: P
        calendar: o.Calendar = Field("JULIAN", alias="Calendar")
        air_snow_coeff: Sym = Field(1 - P.X1, alias="AirSnowCoeff")

    t = TestConfig(params=[0.5], Calendar="NOLEAP")
    assert t.air_snow_coeff == 0.5


def test_gr4j(tmpdir):
    from ravenpy.new_config.emulators.gr4jcn import GR4JCN

    m = GR4JCN(params=[1, 2, 3, 4, 5, 6])
    m.write("/tmp/test_config")
    print(tmpdir)


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

    class Test(Command):
        vegetation_classes: VegetationClasses = Field(
            VegetationClasses.parse_obj([{"name": "VEG_ALL"}, {"name": "VEG_WATER"}]),
            alias="VegetationClasses",
        )

    t = Test()
    t.to_rv()


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

    class TestConfig(Config):
        params: P
        hydrologic_processes: Sequence[Process] = Field(hp, alias="HydrologicProcess")

    t = TestConfig(params=[0.5])

    out = t.to_rv("rvi")
    # assert out == "\n:HydrologicProcesses\n  :Precipitation ATMOS_PRECIP MULTIPLE  \n  :SnowTempEvolve
    # SNOW_TEMP\n:EndHydrologicProcesses\n"


def a_test_defaults():
    import pint

    ureg = pint.UnitRegistry()

    for name, u in defaults.units.items():
        ureg(u)
