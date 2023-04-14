import datetime as dt
from typing import Union

import cftime
import pytest
from pydantic import Field, ValidationError
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

from ravenpy.config import commands as rc
from ravenpy.config.rvs import RVI, Config


def test_rvi_datetime():
    exp = cftime.datetime(1990, 1, 1, calendar="PROLEPTIC_GREGORIAN")

    rvi = RVI(start_date=dt.datetime(1990, 1, 1))
    assert rvi.start_date == exp

    rvi = RVI(start_date="1990-01-01")
    assert rvi.start_date == exp

    rvi = RVI(start_date="1990-01-01T00:00:00")
    assert rvi.start_date == exp

    rvi = RVI(end_date=cftime.datetime(1990, 1, 1))
    assert rvi.end_date == exp

    with pytest.raises(ValidationError):
        rvi = RVI(start_date=(dt.datetime(1990, 1, 1),))


def test_duplicate():
    conf = Config(start_date="1990-01-01")

    # Updating values with an alias and an attribute name
    out = conf.duplicate(Duration=10, debug_mode=True)

    assert isinstance(out.start_date, cftime.datetime)
    assert out.duration == 10
    assert out.debug_mode


def test_set_params():
    @dataclass(config=dict(arbitrary_types_allowed=True))
    class P:
        X01: Union[Variable, float] = Variable("X01")

    class MySymbolicEmulator(Config):
        params: P = P()
        rain_snow_transition: rc.RainSnowTransition = Field(
            default=rc.RainSnowTransition(temp=P.X01, delta=2),
            alias="RainSnowTransition",
        )

    exp = MySymbolicEmulator(params=[0.5])
    s = MySymbolicEmulator()
    assert s.set_params([0.5]) == exp
    s.params = [0.5]
    assert s.rain_snow_transition == exp.rain_snow_transition
    assert s.rvp == exp.rvp


def test_solution(get_local_testdata):
    sol = get_local_testdata("gr4j_cemaneige/solution.rvc")
    conf = Config().set_solution(sol)
    assert len(conf.hru_state_variable_table) == 1
    assert conf.hru_state_variable_table[0].data["ATMOSPHERE"] == 821.98274
    assert conf.hru_state_variable_table[0].data["ATMOS_PRECIP"] == -1233.16

    assert len(conf.basin_state_variables) == 1
    assert conf.basin_state_variables[0].channel_storage == 0
    assert conf.basin_state_variables[0].qout == (1, 13.21660, 13.29232)

    assert ":BasinIndex 1 watershed" in conf.rvc


def test_rvh_from_extractor(get_local_testdata):
    from ravenpy.extractors import BasinMakerExtractor, open_shapefile

    shp = get_local_testdata(
        "basinmaker/drainage_region_0175_v2-1/finalcat_info_v2-1.zip"
    )
    bm = BasinMakerExtractor(open_shapefile(shp))

    # Smoke test
    Config(**bm.extract(hru_from_sb=True))


def test_config(dummy_config):
    cls, P = dummy_config
    conf = cls(Calendar="NOLEAP")

    assert conf.__config__.allow_mutation

    with pytest.raises(ValueError):
        assert conf.rvi

    # Set params
    num = conf.set_params([0.5])
    assert num.air_snow_coeff == 0.5

    with pytest.raises(ValueError):
        num.set_params([0.6])

    # Instantiate with numerical params
    nt = cls(params=[0.5], Calendar="NOLEAP")
    assert nt.air_snow_coeff == 0.5
