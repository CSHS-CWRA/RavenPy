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

salmon_land_hru_1 = dict(
    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659, hru_type="land"
)


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
