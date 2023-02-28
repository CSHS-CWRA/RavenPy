from ravenpy.new_config.rvs import Config
from ravenpy.config import options as o
from pydantic import Field
from pydantic.dataclasses import dataclass

def test_emulator():
    from ravenpy.new_config.base import Sym, SymConfig, Variable

    @dataclass(config=SymConfig)
    class P:
        X1: Sym = Variable("X1")

    class TestConfig:
        params: P
        calendar: o.Calendar = Field("JULIAN", alias="Calendar")
        air_snow_coeff = Field(1 - P.X1, alias="AirSnowCoeff")

    t = TestConfig(params=[.5], calendar="NOLEAP")




def test_gr4j():
    from ravenpy.new_config.emulators import GR4JCN
    GR4JCN(params=[1,2,3,4,5,6])


def test_gr4jcn(tmpdir):
    from ravenpy.models.emulators import GR4JCN
    m = GR4JCN(params=[1,2,3,4,5,6])
    m.write(tmpdir)


def test_defaults():
    import pint

    ureg = pint.UnitRegistry()

    for name, u in defaults.units.items():
        ureg(u)
