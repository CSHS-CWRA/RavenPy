import datetime as dt

from pydantic import Field
from pydantic.dataclasses import dataclass

from ravenpy.new_config import commands as rc
from ravenpy.new_config import options as o
from ravenpy.new_config.base import Sym, SymConfig, Variable, encoder
from ravenpy.new_config.emulators.gr4jcn import GR4JCN
from ravenpy.new_config.rvs import Config

alt_names = {
    "PRECIP": "rain",
    "TEMP_MIN": "tmin",
    "TEMP_MAX": "tmax",
    "PET": "pet",
    "HYDROGRAPH": "qobs",
}


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


def test_gr4j(gr4jcn_config_rv):
    assert (gr4jcn_config_rv / "test.rvi").exists()
