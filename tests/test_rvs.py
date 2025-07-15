import datetime as dt
import pathlib
from typing import Union

import cftime
import pytest
from pydantic import ConfigDict, Field, ValidationError
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

from ravenpy.config import commands as rc
from ravenpy.config import options
from ravenpy.config.rvs import RV, RVI, Config, optfield


def test_optfield():
    class Test(RV):
        a: bool = optfield(alias="a")

    t = Test()
    assert not t.__class__.model_fields["a"].is_required()


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

    rvi = RVI(start_date=dt.datetime(1990, 1, 1), calendar="NOLEAP")
    assert rvi.start_date == cftime.datetime(1990, 1, 1, calendar="NOLEAP")

    with pytest.raises(ValidationError):
        rvi = RVI(start_date=(dt.datetime(1990, 1, 1),))


def test_duplicate_simple():
    conf = Config(start_date="1990-01-01")

    # Updating values with an alias and an attribute name
    out = conf.duplicate(Duration=10, debug_mode=True)

    assert isinstance(out.start_date, cftime.datetime)
    assert out.duration == 10
    assert out.debug_mode


def test_duplicate_emulator(gr4jcn_config):
    conf, params = gr4jcn_config
    conf.duplicate()

    conf.duplicate(params=params)


def test_set_params():
    @dataclass(config=ConfigDict(arbitrary_types_allowed=True))
    class P:
        X01: Union[Variable, float] = Variable("X01")

    class MySymbolicEmulator(Config):
        params: P = P()
        rain_snow_transition: rc.RainSnowTransition = Field(
            default=rc.RainSnowTransition(temp=P.X01, delta=2),
            alias="RainSnowTransition",
        )

    # Assignment through set_params -> new instance
    s = MySymbolicEmulator()
    num = s.set_params([0.5])
    assert num.rain_snow_transition.temp == 0.5

    s = MySymbolicEmulator()
    num = s.set_params({"X01": 0.5})
    assert num.rain_snow_transition.temp == 0.5


def test_solution(yangtze):
    sol = pathlib.Path(yangtze.fetch("gr4j_cemaneige/solution.rvc"))
    conf = Config().set_solution(sol)
    assert len(conf.hru_state_variable_table) == 1
    assert conf.hru_state_variable_table[0].data["ATMOSPHERE"] == 821.98274
    assert conf.hru_state_variable_table[0].data["ATMOS_PRECIP"] == -1233.16

    assert len(conf.basin_state_variables) == 1
    assert conf.basin_state_variables[0].channel_storage == 0
    assert conf.basin_state_variables[0].qout == (1, 13.21660, 13.29232)

    assert ":BasinIndex 1 watershed" in conf.rvc


def test_rvh_from_extractor(yangtze):
    from ravenpy.extractors import BasinMakerExtractor, open_shapefile

    shp = yangtze.fetch("basinmaker/drainage_region_0175_v2-1/finalcat_info_v2-1.zip")
    bm = BasinMakerExtractor(open_shapefile(shp))

    # Smoke test
    Config(**bm.extract(hru_from_sb=True))


def test_config(dummy_config):
    cls, P = dummy_config
    conf = cls(Calendar="NOLEAP")

    # Set params
    num = conf.set_params([0.5])
    assert num.air_snow_coeff == 0.5

    with pytest.raises(ValueError):
        num.set_params([0.6])

    # Instantiate with numerical params
    assert conf.model_config["populate_by_name"]
    # nt = cls(params=[0.5], Calendar="NOLEAP")
    # assert nt.air_snow_coeff == 0.5


def test_custom_subclass(dummy_config, tmp_path):
    """Test that users can subclass RV and Config."""
    cls, P = dummy_config

    # Custom RVI
    class myRVI(RVI):
        run_name: str = Field("myRunName", alias="RunName")

    # Custom config with custom RVI
    class MyConfig(myRVI, cls):
        params: P = P()
        enkf_mode: options.EnKFMode = optfield(alias="EnKFMode")

    # Make sure rv files can be written
    conf = MyConfig(EnKFMode="ENKF_SPINUP").set_params([0.5])
    conf.write_rv(workdir=tmp_path)
    assert conf.run_name == "myRunName"
    assert "myRunName" in conf._rv("RVI")
    assert "EnKFMode" not in conf._rv("RVI")


def test_hru_filter():
    """Test that unrecognized HRU types are filtered out."""
    from ravenpy.config.emulators.gr4jcn import LakeHRU
    from ravenpy.config.emulators.hbvec import HRUs, LandHRU

    with pytest.warns(UserWarning):
        hrus = HRUs([LandHRU(), LandHRU(), LakeHRU()])

    # The GR4J lake HRU is not part of the HBVEC config.
    assert len(hrus.root) == 2
