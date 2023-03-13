import datetime as dt
from enum import Enum
from textwrap import dedent, indent
from typing import Any, ClassVar, Dict, List, Literal, Sequence, Tuple, TypeVar, Union

from pymbolic.primitives import Expression, Variable

"""
Notes
-----
In pydantic, BaseModel and dataclasses behave slightly differently.
A dataclass object can be instantiated using positional arguments, while BaseModels cannot.
This is why Params is defined as a dataclass

Two cases
RV
    Command, alias="A1"
    Sequence[Command], alias="A2"

Sometimes, we want A2 to be printed as a command (HydrologicProcesses)
Sometime, we don't (Gauges)
In general, let's print the alias except if the parent is __root__, or if the parent's name is the same as the
children (Data holding multiple Data elements).

"""
import pytest
from pydantic import BaseModel, Extra, Field, ValidationError, validator


class SymConfig:
    arbitrary_types_allowed = True


Sym = Union[Variable, Expression, float, None]


class Params:
    pass


def encoder(v: dict) -> dict:
    """
    Return string representation of objects in dictionary.

    This is meant to be applied to BaseModel attributes that either have an `alias` defined,
    or have a `__root__` attribute. The objective is to avoid creating `Command` objects for every
    configuration option.

    - bool: ':{cmd}\n' if obj else ''
    - dict: ':{cmd} {key} {value}'
    - enum: ':{cmd} {obj.value}'
    - Command: obj.to_rv()
    - Sequence: complicated
    - Any other: ':{cmd} {obj}'

    """
    import warnings

    for cmd, obj in v.items():
        if obj is None:
            continue
        elif isinstance(obj, bool):
            # :Command\n
            out = f":{cmd}\n" if obj else ""
        elif isinstance(obj, dict):
            # :Command key value\n
            out = (
                "\n".join([f":{cmd} {key} {value}" for (key, value) in obj.items()])
                + "\n\n"
            )
        elif isinstance(obj, Enum):
            # :Command value\n
            out = f":{cmd:<20} {obj.value}\n"
        elif hasattr(obj, "to_rv"):
            # Custom
            out = obj.to_rv()
        elif isinstance(obj, (list, tuple)):
            s = map(str, obj)
            o0 = obj[0]
            if cmd == "__root__" or issubclass(o0.__class__, FlatCommand):
                out = "\n".join(s)
            elif issubclass(o0.__class__, (Command, Record)):
                rec = indent("\n".join(s), Command._indent)
                out = f":{cmd}\n{rec}\n:End{cmd}\n"
            elif isinstance(o0, Enum):
                seq = " ".join([o.value for o in obj])
                out = f":{cmd:<20} {seq}\n"
            else:
                seq = " ".join(s)
                out = f":{cmd:<20} {seq}\n"
        else:
            out = f":{cmd:<20} {obj}\n"

        v[cmd] = str(out)

    return v


class Record(BaseModel):
    class Config:
        extra = Extra.forbid
        arbitrary_types_allowed = True
        allow_population_by_field_name = True


class Command(BaseModel):
    """
    Base class for Raven commands.
    """

    _template: str = """
            :{_cmd}
            {_commands}{_records}
            :End{_cmd}
            """
    _indent: str = "  "

    def __str__(self):
        return self.to_rv()

    def _models(self) -> Dict[str, str]:
        """Return dictionary of RV strings for class attributes that are Raven models."""
        d = {}
        for key, field in self.__fields__.items():
            if field.has_alias or field.alias == "__root__":
                obj = self.__dict__[key]
                if obj is not None:
                    d[field.alias] = self.__dict__[key]
        return encoder(d)

    def _records(self) -> List[str]:
        """Return list of RV strings for records."""
        return [
            str(o)
            for o in self.__dict__.get("__root__", [])
            if issubclass(o.__class__, Record)
        ]

    def to_rv(self):
        """Return Raven configuration string."""
        d = self.dict()
        cmds = self._models()
        d.update(cmds)

        recs = "\n".join(self._records())
        if recs:
            recs = indent(recs, self._indent)

        d["_cmd"] = self.__class__.__name__
        d["_commands"] = indent("".join(cmds.values()).strip(), self._indent)
        if cmds and recs and "_commands" in self._template:
            recs = "\n" + recs
        d["_records"] = recs
        return dedent(self._template).format(**d)

    class Config:
        extra = Extra.forbid
        arbitrary_types_allowed = True
        allow_population_by_field_name = True


class FlatCommand(Command):
    """Only used to discriminate Commands that should not be nested."""


class ParameterList(Command):
    name: str = ""
    values: Sequence[Union[Sym, None]] = ()

    @validator("values", pre=True)
    def no_none_in_default(cls, v, values):
        """Make sure that no values are None for the [DEFAULT] record."""
        if values["name"] == "[DEFAULT]" and None in v:
            raise ValueError("Default record can not contain None.")
        return v

    def to_rv(self, **kwds):
        fmt = "{name:<16}" + len(self.values) * ",{:>18}"
        evals = []
        for v in self.values:
            ev = "_DEFAULT" if v is None else v
            evals.append(ev)

        return fmt.format(name=self.name, *evals)


class GenericParameterList(Command):
    names: Sequence[str] = Field(None, description="Parameter names")
    pl: Sequence[ParameterList] = ()

    @validator("pl")
    def num_values_equal_num_names(cls, val, values):
        n = len(values["names"])
        for v in val:
            if len(v.values) != n:
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
            _commands=indent("\n".join([rec.to_rv() for rec in self.pl]), "  "),
        )


class RV(Command):
    """"""

    _template: str = "{_commands}"
    _indent: str = ""

    class Config:
        allow_population_by_field_name = True
        validate_all = True


def parse_symbolic(value, **kwds):
    """Inject values of symbolic variables into object and return object."""
    from pymbolic.mapper.evaluator import EvaluationMapper as EM
    from pymbolic.primitives import Expression, Variable

    if isinstance(value, dict):
        return {k: parse_symbolic(v, **kwds) for k, v in value.items()}

    elif isinstance(value, (list, tuple)):
        return [parse_symbolic(v, **kwds) for v in value]

    elif isinstance(value, (Command, Record)):
        attrs = value.dict()
        return value.__class__(**parse_symbolic(attrs, **kwds))

    elif isinstance(value, (Variable, Expression)):
        # Inject numerical values numerical value
        return EM(context=kwds)(value)

    else:
        return value


# ---
