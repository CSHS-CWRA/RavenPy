import datetime as dt
from pathlib import Path
from textwrap import dedent
from typing import Sequence

import pytest
from pydantic import Field

from ravenpy.new_config.base import RV, Command, RecordCommand
from ravenpy.new_config.commands import (
    HRU,
    PL,
    BasinIndex,
    BasinStateVariables,
    ChannelProfile,
    CustomOutput,
    Data,
    EvaluationPeriod,
    Gauge,
    GriddedForcing,
    GridWeights,
    HRUsCommand,
    HRUState,
    HRUStateVariableTable,
    LandUseParameterList,
    LinearTransform,
    ObservationData,
    RainSnowTransition,
    ReadFromNetCDF,
    RedirectToFile,
    Reservoir,
    SoilClasses,
    SoilProfile,
    SoilProfiles,
    StationForcing,
    SubBasin,
    SubBasinGroup,
    SubBasins,
)


def test_linear_transform():
    class TestRV(RV):
        lt: LinearTransform = Field(None, alias="LinearTransform")

    t = TestRV(LinearTransform={"scale": 2, "offset": 5})
    assert t.commands() == ":LinearTransform 2.000000000000000 5.000000000000000\n"


def test_rain_snow_transition():
    class TestRV(RV):
        rst: RainSnowTransition = Field(None, alias="RainSnowTransition")

    t = TestRV(RainSnowTransition={"temp": 2, "delta": 5})
    assert t.commands() == ":RainSnowTransition 2.0 5.0\n"


def test_evaluation_period():
    class TestRV(RV):
        ep: EvaluationPeriod = Field(None, alias="EvaluationPeriod")

    t = TestRV(
        EvaluationPeriod={
            "name": "WET",
            "start": dt.datetime(2001, 1, 1),
            "end": dt.datetime(2002, 1, 1),
        }
    )

    assert t.commands() == ":EvaluationPeriod WET 2001-01-01 2002-01-01"


def test_custom_output():
    class TestRV(RV):
        co: CustomOutput = Field(None, alias="CustomOutput")

    t = TestRV(
        CustomOutput={
            "time_per": "DAILY",
            "stat": "AVERAGE",
            "variable": "RAIN",
            "space_agg": "BY_HRU",
        }
    )

    assert t.commands() == ":CustomOutput DAILY AVERAGE RAIN BY_HRU \n"


def test_soil_profiles():
    class TestRV(RV):
        soil_profiles: SoilProfiles = Field(None, alias="SoilProfiles")

    t = TestRV(
        SoilProfiles=(
            SoilProfile(
                name="DEFAULT", soil_classes=["FAST", "SLOW"], thicknesses=[10.0, 100.0]
            ),
        )
    )
    assert dedent(t.commands()) == dedent(
        """
        :SoilProfiles
          DEFAULT         ,   2,        FAST,  10.0,        SLOW, 100.0
        :EndSoilProfiles
        """
    )


def test_subbasins():
    class TestRV(RV):
        sb: SubBasins = Field(None, alias="SubBasins")

    t = TestRV(
        SubBasins=[
            SubBasin(),
        ]
    )

    expected = """
        :SubBasins
          :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
          :Units      none none          none    none           km    none
          1          sub_001    -1         None       ZERO-      1
        :EndSubBasins
        """
    for (a, e) in zip(t.commands().splitlines(), expected.splitlines()):
        assert a.strip() == e.strip()


def test_soil_classes():
    class TestRV(RV):
        sc: SoilClasses = Field(None, alias="SoilClasses")

    t = TestRV(SoilClasses=["TOPSOIL", "FAST_RES", "SLOW_RES"])

    assert dedent(t.commands()) == dedent(
        """
    :SoilClasses
      TOPSOIL
      FAST_RES
      SLOW_RES
    :EndSoilClasses
    """
    )


def test_hrus():
    class TestRV(RV):
        hrus: HRUsCommand = Field(..., alias="HRUs")

    t = TestRV(
        HRUs=[
            HRU(),
        ]
    )
    t.commands()


def test_reservoir():
    class TestRV(RV):
        r: Sequence[Reservoir] = Field(None, alias="Reservoir")

    t = TestRV(Reservoir=[Reservoir()])
    assert t.commands().strip().startswith(":Reservoir")


def test_subbasin_group():
    class TestRV(RV):
        group: SubBasinGroup = Field(None, alias="SubBasinGroup")

    t = TestRV(SubBasinGroup=SubBasinGroup(name="Group", sb_ids=[1, 2, 3]))
    assert dedent(t.commands()) == dedent(
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
    t.commands()


def test_read_from_netcdf():
    nc = ReadFromNetCDF(FileNameNC="test.nc", VarNameNC="pr", DimNamesNC=("time",))
    print(nc.commands())


def test_data():
    class TestRV(RV):
        data: Sequence[Data] = Field(None, alias="Data")

    nc = ReadFromNetCDF(FileNameNC="test.nc", VarNameNC="pr", DimNamesNC=("time",))
    t = TestRV(
        Data=[
            Data(ReadFromNetCDF=nc),
        ]
    )
    t.to_rv()


def test_observation_data():
    class TestRV(RV):
        data: Sequence[ObservationData] = Field(None, alias="ObservationData")

    nc = ReadFromNetCDF(FileNameNC="test.nc", VarNameNC="pr", DimNamesNC=("time",))
    t = TestRV(
        ObservationData=[
            ObservationData(uid=12, ReadFromNetCDF=nc),
        ]
    )
    t.to_rv()


def test_gauge():
    class TestRV(RV):
        gauge: Gauge = Field(None, alias="Gauge")

    nc = ReadFromNetCDF(FileNameNC="test.nc", VarNameNC="pr", DimNamesNC=("time",))
    t = TestRV(
        Gauge=Gauge(name="my_gauge", Latitude=56, Data=[Data(ReadFromNetCDF=nc)])
    )
    t.to_rv()


def test_grid_weights():
    class TestRV(RV):
        gw: GridWeights = Field(None, alias="GW")

    t = TestRV(
        GW=GridWeights(
            NumberHRUs=3,
            NumberGridCells=1,
            data=((1, 0, 1.0), (2, 0, 1.0), (3, 0, 1.0)),
        )
    )
    print(t.to_rv())


def test_redirect_to_file(tmpdir):
    gw = GridWeights()
    gw_path = tmpdir / Path("grid_weights.rvt")
    gw_path.write_text(gw.to_rv() + "\n", "utf8")
    rtf = RedirectToFile(path=gw_path)
    rtf.to_rv()


def test_gridded_forcing():
    nc = dict(FileNameNC="test.nc", VarNameNC="pr", DimNamesNC=("time",))
    gf = GriddedForcing(**nc)
    print(gf.to_rv())


def test_station_forcings():
    nc = dict(FileNameNC="test.nc", VarNameNC="pr", DimNamesNC=("time",))
    sf = StationForcing(ForcingType="PRECIP", **nc)
    print(sf.to_rv())


def test_hru_state():
    s = HRUState(index=1, data={"SOIL[0]": 1, "SOIL[1]": 2.0})
    assert s.to_rv() == "1,1.0,2.0"


def test_hru_state_variable_table():
    class TestRV(RV):
        hru: HRUStateVariableTable = Field(None, alias="HRUStateVariableTable")

    s1 = HRUState(index=1, data={"SOIL[0]": 0.1, "SOIL[1]": 1.0})

    t = TestRV(
        HRUStateVariableTable=[
            s1,
        ]
    )
    rv = t.commands()
    assert dedent(rv) == dedent(
        """
            :HRUStateVariableTable
              :Attributes,SOIL[0],SOIL[1]
              :Units
              1,0.1,1.0
            :EndHRUStateVariableTable
            """
    )

    s2 = HRUStateVariableTable.parse(rv)
    assert s2.to_rv() == rv


@pytest.mark.xfail
def test_hru_state_variable_table_non_uniform():
    class TestRV(RV):
        hru: HRUStateVariableTable = Field(None, alias="HRUStateVariableTable")

    s1 = HRUState(index=1, data={"SOIL[0]": 0.1, "SOIL[1]": 1.0})
    s2 = HRUState(index=2, data={"SOIL[3]": 3, "SOIL[2]": 2.0})

    t = TestRV(HRUStateVariableTable=[s1, s2])
    assert dedent(t.commands()) == dedent(
        """
            :HRUStateVariableTable
              :Attributes,SOIL[0],SOIL[1],SOIL[2],SOIL[3]
              1,0.1,1.0,0.0,0.0
              2,0.0,0.0,2.0,3.0
            :EndHRUStateVariableTable
            """
    )


def test_land_use_parameter_list():
    class TestRV(RV):
        land_use_parameter_list: LandUseParameterList = Field(
            None, alias="LandUseParameterList"
        )

    t = TestRV(
        LandUseParameterList=LandUseParameterList(
            names=["MELT_FACTOR", "AET_COEFF", "FOREST_SPARSENESS", "DD_MELT_TEMP"],
            records=[
                PL(name="DEFAULT", vals=[1, 2, 3, 4]),
                PL(name="PATATE", vals=[10, 2, 3, 4]),
            ],
        ),
    )
    t.to_rv()


def test_basin_index():
    bi = BasinIndex()
    rv = bi.to_rv()

    bi2 = BasinIndex.parse(rv)
    assert bi2.to_rv() == rv


def test_basin_state_variables():
    bs = BasinStateVariables.parse_obj(
        [
            BasinIndex(),
        ]
    )
    rv = bs.to_rv()

    bs2 = BasinStateVariables.parse(rv)
    assert bs2.to_rv() == rv
