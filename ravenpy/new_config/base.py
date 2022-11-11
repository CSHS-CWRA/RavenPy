import datetime as dt
from enum import Enum
from textwrap import dedent, indent
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
                seq = ", ".join(map(str, obj))
                out = f":{cmd:<20} {seq}\n"
        else:
            out = f":{cmd:<20} {obj}\n"

        v[cmd] = str(out)

    return v


class Command(BaseModel):
    """
    Base class for Raven commands.
    """

    _template: str = """
            :{_cmd}
            {_commands}
            :End{_cmd}
            """

    def __str__(self):
        return self.to_rv()

    def command_objs(self):
        """Return attributes that are Raven commands."""
        d = {}
        for key, field in self.__fields__.items():
            if field.has_alias or field.alias == "__root__":
                if self.__dict__[key] is not None:
                    d[field.alias] = self.__dict__[key]
        return d

    def command_json(self):
        """Return dictionary of Raven commands."""
        return encoder(self.command_objs())

    def commands(self):
        """String of Raven commands."""
        return "".join(self.command_json().values())

    def to_rv(self):
        """Return Raven configuration string."""
        d = self.dict()
        d.update(self.command_json())
        d["_cmd"] = self.__class__.__name__
        d["_commands"] = indent(self.commands().strip(), "  ")
        return dedent(self._template).format(**d)

    class Config:
        extra = Extra.forbid
        arbitrary_types_allowed = True


class Record(Command):
    __root__: List


class RecordCommand(tuple):
    record: TypeVar
    _template: str = """
        :{_cmd}
        {_commands}
        :End{_cmd}
        """

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, values):
        for v in values:
            if not isinstance(v, cls.record):
                raise TypeError
        return cls(values)

    def command_objs(self):
        """Return attributes that are Raven commands."""
        d = {}
        for key, field in self.__dict__.items():
            if field.has_alias and self.__dict__[key] is not None:
                d[field.alias] = self.__dict__[key]

        for i, v in enumerate(self):
            d[i] = v

        return d

    def command_json(self):
        """Return dictionary of Raven commands."""
        return encoder(self.command_objs())

    def commands(self):
        """String of Raven commands."""
        return "".join(self.command_json().values())

    def to_rv(self):
        """Return Raven configuration string."""
        d = self.__dict__.copy()
        d.update(self.command_json())
        d["_cmd"] = self.__class__.__name__
        d["_commands"] = indent(self.commands().strip(), "  ")
        return dedent(self._template).format(**d)

    def to_rv_old(self):
        d = {
            "_cmd": self.__class__.__name__,
            "_commands": indent("\n".join(map(str, self)), "  "),
        }
        return dedent(self._template).format(**d)


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
            :{_cmd}
              :Parameters     {name_list}
              :Units          {unit_list}
            {_commands}
            :End{_cmd}
        """

        fmt = ",{:>18}" * len(self.names)
        units = [
            "none",
        ] * len(self.names)

        return dedent(template).format(
            _cmd=self.__class__.__name__,
            name_list=fmt.format(*self.names),
            unit_list=fmt.format(*units),
            _commands=indent("\n".join([rec.to_rv() for rec in self.records]), "  "),
        )


class RV(Command):
    """"""


# ---
