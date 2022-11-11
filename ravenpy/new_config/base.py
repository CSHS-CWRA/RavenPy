import datetime as dt
from enum import Enum
from textwrap import dedent
from typing import Any, ClassVar, Dict, List, Literal, Sequence, Tuple, TypeVar, Union

import pytest
from pydantic import BaseModel, Extra, Field, ValidationError, validator


def encoder(v):
    for cmd, obj in v.items():
        if obj is None:
            out = ""
        elif isinstance(obj, bool):
            out = f":{cmd}\n" if obj else ""
        elif isinstance(obj, dict):
            out = "\n".join([f":{cmd} {key} {value}" for (key, value) in obj.items()])
        elif isinstance(obj, Enum):
            out = f":{cmd:<20} {obj.value}\n"
        elif isinstance(obj, (Command, RecordCommand)):
            out = obj.to_rv()
        elif isinstance(obj, (list, tuple)):
            o0 = obj[0]
            if isinstance(o0, Command):
                out = "".join([o.to_rv() for o in obj])
            elif isinstance(o0, Enum):
                seq = ", ".join([o.value for o in obj])
                out = f":{cmd:<20} {seq}\n"
            else:
                seq = ", ".join(obj)
                out = f":{cmd:<20} {seq}\n"
        else:
            out = f":{cmd:<20} {obj}\n"

        v[cmd] = str(out)

    return v


class Command(BaseModel):
    """
    Base class for Raven commands.
    """

    def __str__(self):
        return self.to_rv()

    def fields(self):
        d = {}
        for key, field in self.__fields__.items():
            if field.has_alias and self.__dict__[key] is not None:
                d[field.alias] = self.__dict__[key]
        return d

    def encode(self):
        return encoder(self.fields())

    def content(self):
        return "".join(self.encode().values())

    def to_rv(self):
        template = """
                :{cmd}
                  {records}
                :End{cmd}
                """
        return dedent(template).format(
            records="\n  ".join(self.encode().values()).strip(),
            cmd=self.__class__.__name__,
        )

    class Config:
        extra = Extra.forbid


class RecordCommand(tuple):
    record: TypeVar
    template: str = """
        :{cmd}
            {records}
        :End{cmd}
        """

    class Config:
        extra = Extra.forbid

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, values):
        for v in values:
            if not isinstance(v, cls.record):
                raise TypeError
        return cls(values)

    def to_rv(self):
        return dedent(self.template).format(
            records="\n  ".join(map(str, self)), cmd=self.__class__.__name__
        )


class ParameterList(Command):
    name: str = ""
    vals: Sequence[Union[float, None]] = ()

    @validator("vals", pre=True)
    def no_none_in_default(cls, v, values):
        """Make sure that no values are None for the [DEFAULT] record."""
        if values["name"] == "[DEFAULT]" and None in v:
            raise ValueError("Default record can not contain None.")
        return v

    def to_rv(self, **kwds):
        fmt = "{name:<16}" + len(self.vals) * ",{:>18}"
        evals = []
        for v in self.vals:
            ev = "_DEFAULT" if v is None else v
            evals.append(ev)

        return fmt.format(name=self.name, *evals)


class ParameterListCommand(Command):
    names: Sequence[str] = ()
    records: Sequence[ParameterList] = ()

    @validator("records")
    def num_values_equal_num_names(cls, val, values):
        n = len(values["names"])
        for v in val:
            if len(v.vals) != n:
                raise ValidationError(
                    "Number of values should match number of parameters."
                )
        return val

    def to_rv(self):
        template: str = """
            :{cmd}
                :Parameters    {name_list}
                :Units         {unit_list}
                {records}
            :End{cmd}
        """

        fmt = ",{:>18}" * len(self.names)
        units = ",              none" * len(self.names)

        return dedent(template).format(
            cmd=self.__class__.__name__,
            name_list=fmt.format(*self.names),
            unit_list=units,
            records="\n".join([rec.to_rv() for rec in self.records]),
        )


class RV(Command):
    """"""


# ---
