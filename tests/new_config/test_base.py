import pytest
from pydantic import Field, ValidationError
from typing import Sequence, Dict
from enum import Enum
import datetime as dt
from textwrap import dedent
from ravenpy.new_config.base import RV, Command


class TestRV():
    def test_bool(self):
        class Test(RV):
            so_true: bool = Field(None, alias="SuppressOutput")
            so_false: bool = Field(None, alias="ThisIsFalse")

        assert Test(SuppressOutput=True).content() == ":SuppressOutput\n"
        assert dedent(Test(SuppressOutput=True).to_rv()) == dedent(
            """
            :Test
              :SuppressOutput
            :EndTest
            """)

    def test_value(self):
        class Test(RV):
            avg: float = Field(None, alias="AvgAnnualSnow")
            start_date: dt.datetime = Field(None, alias="StartDate")

        assert Test(AvgAnnualSnow=45.4).content() == ":AvgAnnualSnow        45.4\n"
        assert Test(StartDate=dt.datetime(2000, 1, 1)).content() == ":StartDate            2000-01-01 00:00:00\n"

    def test_dict(self):
        class Test(RV):
            gp: Dict = Field(None, alias="GlobalParameters")

        t = Test(GlobalParameters={"AvgAnnualSnow": 34, "PrecipitationLapseRate": 12.})
        assert t.content() == ":GlobalParameters AvgAnnualSnow 34\n:GlobalParameters PrecipitationLapseRate 12.0"

    def test_option(self):
        class Opt(Enum):
            a = "A"
            b = "B"

        class Test(RV):
            opt: Opt = Field(None, alias="Option")

        assert Test(Option="A").content() == ":Option               A\n"
        assert Test(Option=Opt.a).content() == ":Option               A\n"

    def test_option_list(self):
        class Opt(Enum):
            A = "A"
            B = "B"

        class Test(RV):
            opts: Sequence[Opt] = Field(None, alias="Options")

        assert Test(Options=["A", "B"]).content() == ":Options              A, B\n"

        with pytest.raises(ValidationError):
            assert Test(Options=["C", "D"])


def test_command():
    class Test(Command):
        a: str = Field(None, alias="Alias")
        b: bool = Field(None, alias="Bool")

    t = Test(Alias="spam", Bool=True)
    assert t.encode() == {"Alias": ":Alias                spam\n",
                                     "Bool": ":Bool\n"}


