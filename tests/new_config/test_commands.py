import datetime as dt
from textwrap import dedent
from typing import Sequence
from pydantic import Field
from ravenpy.new_config.base import RV, Command, RecordCommand
from ravenpy.new_config.commands import LinearTransform, RainSnowTransition, CustomOutput, SoilProfile, SoilProfiles, \
    EvaluationPeriod, SubBasin, SubBasins, HRU, HRUsCommand, Reservoir, ReadFromNetCDF, SubBasinGroup, \
    ChannelProfile, Data


def test_linear_transform():
    class TestRV(RV):
        lt: LinearTransform = Field(None, alias="LinearTransform")

    t = TestRV(LinearTransform={"scale": 2, "offset": 5})
    assert t.to_rv() == ":LinearTransform 2.000000000000000 5.000000000000000\n"


def test_rain_snow_transition():
    class TestRV(RV):
        rst: RainSnowTransition = Field(None, alias="RainSnowTransition")

    t = TestRV(RainSnowTransition={"temp": 2, "delta": 5})
    assert t.to_rv() == ":RainSnowTransition 2.0 5.0\n"


def test_evaluation_period():
    class TestRV(RV):
        ep: EvaluationPeriod = Field(None, alias="EvaluationPeriod")

    t = TestRV(EvaluationPeriod={"name": "WET", "start": dt.datetime(2001, 1, 1), "end": dt.datetime(2002, 1, 1)})

    assert t.to_rv() == ":EvaluationPeriod WET 2001-01-01 2002-01-01"


def test_custom_output():
    class TestRV(RV):
        co: CustomOutput = Field(None, alias="CustomOutput")

    t = TestRV(CustomOutput={"time_per": "DAILY", "stat": "AVERAGE", "variable": "RAIN", "space_agg": "BY_HRU"})

    assert t.to_rv() == ":CustomOutput DAILY AVERAGE RAIN BY_HRU \n"


def test_soil_profiles():

    class TestRV(RV):
        soil_profiles: SoilProfiles = Field(None, alias="SoilProfiles")

    t = TestRV(SoilProfiles=(SoilProfile(name="DEFAULT", soil_classes=["FAST", "SLOW"], thicknesses=[10.0, 100.0]),))
    assert dedent(t.to_rv()) == dedent(
        """
        :SoilProfiles
            DEFAULT         ,   2,        FAST,  10.0,        SLOW, 100.0
        :EndSoilProfiles
        """)


def test_subbasins():
    class TestRV(RV):
        sb: SubBasins = Field(None, alias="SubBasins")

    t = TestRV(SubBasins=[SubBasin(),])

    expected = """
        :SubBasins
          :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
          :Units      none none          none    none           km    none
          1          sub_001    -1         None       ZERO-      1
        :EndSubBasins
        """
    for (a, e) in zip(t.to_rv().splitlines(), expected.splitlines()):
        assert a.strip() == e.strip()



def test_hrus():
    class TestRV(RV):
        hrus: HRUsCommand = Field(..., alias="HRUs")

    t = TestRV(HRUs=[HRU(),])
    t.to_rv()


def test_reservoir():
    class TestRV(RV):
        r: Sequence[Reservoir] = Field(None, alias="Reservoir")

    t = TestRV(Reservoir=[Reservoir()])
    assert t.to_rv().strip().startswith(":Reservoir")

def test_subbasin_group():
    class TestRV(RV):
        group: SubBasinGroup = Field(None, alias="SubBasinGroup")

    t = TestRV(SubBasinGroup=SubBasinGroup(name="Group", sb_ids=[1,2,3]))
    assert dedent(t.to_rv()) == dedent(
        """
        :SubBasinGroup Group
          1, 2, 3
        :EndSubBasinGroup
        """
    )

def test_channel_profile():
    class TestRV(RV):
        cp: Sequence[ChannelProfile] = Field(None, alias="ChannelProfile")

    t = TestRV(ChannelProfile=[ChannelProfile()])
    t.to_rv()

def test_read_from_netcdf():
    nc = ReadFromNetCDF(FileNameNC="test.nc", VarNameNC="pr", DimNamesNC=("time",))
    print(nc.to_rv())

def no_test_data():
    class TestRV(RV):
        data: Sequence[Data] = Field(None, alias="Data")

    t = TestRV(Data=(Data(),))
    t.to_rv()
