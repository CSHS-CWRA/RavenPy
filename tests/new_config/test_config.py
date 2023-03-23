import datetime as dt

import pytest
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

    # Test with symbolic params
    st = TestConfig(params=P(), Calendar="NOLEAP")

    with pytest.raises(ValueError):
        assert st.rvi

    # Set params
    t1 = st.set_params([0.5])
    assert t1.air_snow_coeff == 0.5

    with pytest.raises(ValueError):
        t1.set_params([0.6])

    # Instantiate with numerical params
    nt = TestConfig(params=[0.5], Calendar="NOLEAP")
    assert nt.air_snow_coeff == 0.5
