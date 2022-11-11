import datetime as dt
from enum import Enum
from textwrap import dedent
from typing import Dict, Sequence

import pytest
from pydantic import Field, ValidationError

from ravenpy.new_config.base import RV, Command


class TestRV:
    def test_bool(self):
        class Test(RV):
            so_true: bool = Field(None, alias="SuppressOutput")
            so_false: bool = Field(None, alias="ThisIsFalse")

        assert Test(SuppressOutput=True).commands() == ":SuppressOutput\n"
        assert dedent(Test(SuppressOutput=True).to_rv()) == dedent(
            """
            :Test
              :SuppressOutput
            :EndTest
            """
        )

    def test_value(self):
        class Test(RV):
            avg: float = Field(None, alias="AvgAnnualSnow")
            start_date: dt.datetime = Field(None, alias="StartDate")

        assert Test(AvgAnnualSnow=45.4).commands() == ":AvgAnnualSnow        45.4\n"
        assert (
            Test(StartDate=dt.datetime(2000, 1, 1)).commands()
            == ":StartDate            2000-01-01 00:00:00\n"
        )

    def test_dict(self):
        class Test(RV):
            gp: Dict = Field(None, alias="GlobalParameters")

        t = Test(GlobalParameters={"AvgAnnualSnow": 34, "PrecipitationLapseRate": 12.0})
        assert (
            t.commands()
            == ":GlobalParameters AvgAnnualSnow 34\n:GlobalParameters PrecipitationLapseRate 12.0"
        )

    def test_option(self):
        class Opt(Enum):
            a = "A"
            b = "B"

        class Test(RV):
            opt: Opt = Field(None, alias="Option")

        assert Test(Option="A").commands() == ":Option               A\n"
        assert Test(Option=Opt.a).commands() == ":Option               A\n"

    def test_option_list(self):
        class Opt(Enum):
            A = "A"
            B = "B"

        class Test(RV):
            opts: Sequence[Opt] = Field(None, alias="Options")

        assert Test(Options=["A", "B"]).commands() == ":Options              A, B\n"

        with pytest.raises(ValidationError):
            assert Test(Options=["C", "D"])


def test_command():
    class Test(Command):
        a: str = Field(None, alias="Alias")
        b: bool = Field(None, alias="Bool")

    t = Test(Alias="spam", Bool=True)
    assert t.command_json() == {
        "Alias": ":Alias                spam\n",
        "Bool": ":Bool\n",
    }
