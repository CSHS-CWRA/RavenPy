from enum import Enum
from textwrap import dedent, indent
from typing import Dict, Sequence, Tuple, Union

from pydantic import BaseModel, Extra, Field, root_validator, validator
from pymbolic.primitives import Expression, Variable

"""
Notes
-----
In pydantic, BaseModel and dataclasses behave slightly differently.
A dataclass object can be instantiated using positional arguments, while BaseModels cannot.
This is why Params is defined as a dataclass::

    Two cases
    RV
        Command, alias="A1"
        Sequence[Command], alias="A2"

Sometimes, we want A2 to be printed as a command (HydrologicProcesses)
Sometime, we don't (Gauges)
In general, let's print the alias except if the parent is __root__, or if the parent's name is the same as the
children (Data holding multiple Data elements).

"""


class SymConfig:
    arbitrary_types_allowed = True


Sym = Union[Variable, Expression, float, None]


class Params:
    pass


def encoder(v: dict) -> dict:
    r"""Return string representation of objects in dictionary.

    This is meant to be applied to BaseModel attributes that either have an `alias` defined,
    or have a `__root__` attribute. The objective is to avoid creating `Command` objects for every
    configuration option:

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
        elif issubclass(obj.__class__, Record):
            out = str(obj)
        elif isinstance(obj, (list, tuple)):
            s = map(str, obj)
            o0 = obj[0]
            if cmd == "__root__" or issubclass(o0.__class__, FlatCommand):
                out = "".join(s)
            # elif cmd == "Attributes":
            #     out = f":{cmd}," + ",".join(s) + "\n"
            elif cmd in ["Parameters", "Attributes", "Units"]:
                fmt = ":{cmd:<15}" + len(obj) * ",{:>18}" + "\n"
                out = fmt.format(cmd=cmd, *obj)
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
    """Base class for Raven commands."""

    _template: str = """
            :{_cmd}
            {_commands}{_records}
            :End{_cmd}
            """
    _indent: str = "  "

    def __str__(self):
        return self.to_rv()

    def __subcommands__(self) -> Tuple[Dict[str, str], list]:
        """Return dictionary of class attributes that are Raven models."""
        cmds = {}
        recs = []
        for key, field in self.__fields__.items():
            obj = self.__dict__[key]
            if obj is not None:
                if issubclass(obj.__class__, Record):
                    recs = [obj]
                elif issubclass(obj.__class__, Command):
                    cmds[field.alias] = self.__dict__[key]
                else:
                    try:
                        o = obj[0]
                    except (TypeError, IndexError, KeyError):
                        o = None
                    if issubclass(o.__class__, Record):
                        recs = obj
                    elif field.has_alias or field.alias == "__root__":
                        cmds[field.alias] = self.__dict__[key]
                # elif (
                #     getattr(field.annotation, "_name", "") == "Sequence"
                #     and len(obj)
                #     and issubclass(obj[0].__class__, Record)
                # ):
                #     recs = obj
                # elif field.has_alias or field.alias == "__root__":
                #     cmds[field.alias] = self.__dict__[key]

        for key, field in self.__private_attributes__.items():
            cmds[key.strip("_")] = getattr(self, key)

        return cmds, recs

    def to_rv(self):
        """Return Raven configuration string."""
        d = self.dict()
        cmds, recs = self.__subcommands__()

        # Write command strings
        d.update(encoder(cmds))

        # Write record strings
        recs = "\n".join(map(str, recs))
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


class ListCommand(Command):
    """Use so that commands with __root__: Sequence[Command] behave like a list."""

    def __iter__(self):
        return iter(self.__root__)

    def __getitem__(self, item):
        return self.__root__[item]

    def __len__(self):
        return len(self.__root__)


class ParameterList(Record):
    name: str = ""
    values: Sequence[Union[Sym, str, None]] = ()

    @validator("values", pre=True)
    def no_none_in_default(cls, v, values):
        """Make sure that no values are None for the [DEFAULT] record."""
        if values["name"] == "[DEFAULT]" and None in v:
            raise ValueError("Default record can not contain None.")
        return v

    def __str__(self):
        fmt = "{name:<16}" + len(self.values) * ",{:>18}"
        evals = []
        for v in self.values:
            ev = "_DEFAULT" if v is None else v
            evals.append(ev)

        return fmt.format(name=self.name, *evals)


class GenericParameterList(Command):
    parameters: Sequence[str] = Field(
        None, alias="Parameters", description="Parameter names"
    )
    units: Sequence[str] = Field(None, alias="Units")
    pl: Sequence[ParameterList] = Field(None)

    @root_validator(pre=True)
    def num_values_equal_num_names(cls, values):
        n = len(values["parameters"])
        pl = values["pl"]
        for v in pl:
            # FIXME: assertions should not be found outside of testing code. Replace with conditional logic.
            assert (
                len(v.values) == n
            ), "Number of values should match number of parameters."

        values["units"] = ("none",) * n
        return values


class RV(Command):
    """Base class for RV configuration objects."""

    _template: str = "{_commands}\n"
    _indent: str = ""

    class Config:
        allow_population_by_field_name = True
        validate_all = True
        validate_assignment = True


def parse_symbolic(value, **kwds):
    """Inject values of symbolic variables into object and return object."""
    from pymbolic.mapper.evaluator import EvaluationMapper as EM
    from pymbolic.primitives import Expression, Variable

    if isinstance(value, dict):
        return {k: parse_symbolic(v, **kwds) for k, v in value.items()}

    elif isinstance(value, (list, tuple)):
        return type(value)(parse_symbolic(v, **kwds) for v in value)

    elif isinstance(value, (Command, Record)):
        attrs = value.__dict__
        return value.__class__(**parse_symbolic(attrs, **kwds))

    elif isinstance(value, (Variable, Expression)):
        # Inject numerical values numerical value
        return EM(context=kwds)(value)

    else:
        return value
