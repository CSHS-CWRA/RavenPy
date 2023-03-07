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


def encoder(v):
    import warnings

    for cmd, obj in v.items():
        if obj is None:
            out = ""
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
            o0 = obj[0]
            if hasattr(o0, "to_rv"):
                c = "\n".join([o.to_rv() for o in obj])
                if cmd == "__root__" or cmd == o0.__class__.__name__:
                    out = c
                else:
                    out = f":{cmd}\n{indent(c, '  ')}\n:End{cmd}\n"
            elif isinstance(o0, Enum):
                seq = " ".join([o.value for o in obj])
                out = f":{cmd:<20} {seq}\n"
            else:
                seq = " ".join(map(str, obj))
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
    _indent: str = "  "

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
        d["_commands"] = indent(self.commands().strip(), self._indent)
        return dedent(self._template).format(**d)

    class Config:
        extra = Extra.forbid
        arbitrary_types_allowed = True
        allow_population_by_field_name = True


class ParameterList(Command):
    name: str = ""
    vals: Sequence[Union[Sym, None]] = ()

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

    elif isinstance(value, Command):
        attrs = value.dict()
        return value.__class__(**parse_symbolic(attrs, **kwds))

    elif isinstance(value, (Variable, Expression)):
        # Inject numerical values numerical value
        return EM(context=kwds)(value)

    else:
        return value


# ---
