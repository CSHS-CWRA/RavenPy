import re
from collections.abc import Sequence
from enum import Enum
from textwrap import dedent, indent
from typing import Any, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, RootModel, model_validator
from pymbolic.primitives import ExpressionNode, Variable

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


SymConfig = ConfigDict(arbitrary_types_allowed=True)


Sym = Union[Variable, ExpressionNode, float, None]


def optfield(**kwds):
    """Shortcut to create an optional field with an alias."""
    return Field(**kwds, default_factory=lambda: None, validate_default=False)


class Params:
    pass


def encoder(v: dict) -> dict:
    r"""
    Return string representation of objects in dictionary.

    This is meant to be applied to BaseModel attributes that either have an `alias` defined,
    or have a `root` attribute. The objective is to avoid creating `Command` objects for every
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
        elif issubclass(obj.__class__, _Record):
            out = str(obj)
        elif isinstance(obj, (list, tuple)):
            s = map(str, obj)
            o0 = obj[0]
            if cmd == "root" or issubclass(o0.__class__, FlatCommand):
                out = "".join(s)
            # elif cmd == "Attributes":
            #     out = f":{cmd}," + ",".join(s) + "\n"
            elif cmd in ["Parameters", "Attributes", "Units"]:
                fmt = ":{cmd:<15}" + len(obj) * ",{:>18}" + "\n"
                out = fmt.format(cmd=cmd, *obj)
            elif issubclass(o0.__class__, (_Command, _Record)):
                rec = indent("\n".join(s), o0._indent)
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


class _Record(BaseModel):
    pass


class Record(_Record):
    """A Record has no nested Command or Record objects. It is typically a list of named
    values on a single line.

    For example, SubBasins is a ListCommand, whose root is a list of `SubBasin` Records.
    """

    model_config = ConfigDict(
        extra="forbid", arbitrary_types_allowed=True, populate_by_name=True
    )


class RootRecord(RootModel, _Record):
    """A Record with a root attribute. This is typically used for an unnamed list of
    values on a single line.

    For example, the list of HRUs in an HRUGroup is a RootRecord, and the weights of a GridWeights
    command are a sequence of records.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True, populate_by_name=True)

    def __iter__(self):
        return iter(self.root)

    def __getitem__(self, item):
        return self.root[item]

    def __str__(self):
        return " ".join(map(str, self.root))


class _Command(BaseModel):
    """Base class for Raven commands."""

    @property
    def _indent(self):
        return "  "

    @property
    def _template(self):
        return """
               :{_cmd}
               {_commands}{_records}
               :End{_cmd}
               """

    def __str__(self):
        return self.to_rv()

    def __subcommands__(self) -> tuple[dict[str, str], list]:
        """Return dictionary of class attributes that are Raven models."""
        cmds = {}
        recs = []
        cls = self.__class__
        for key, field in cls.model_fields.items():
            obj = self.__dict__[key]
            if obj is not None:
                if issubclass(obj.__class__, _Record):
                    recs = [obj]
                elif issubclass(obj.__class__, _Command):
                    cmds[field.alias] = self.__dict__[key]
                else:
                    try:
                        o = obj[0]
                    except (TypeError, IndexError, KeyError):
                        o = None
                    if issubclass(o.__class__, _Record):
                        recs = obj
                    elif field.alias is not None:
                        cmds[field.alias] = self.__dict__[key]
                    elif key == "root":
                        cmds[key] = self.__dict__[key]
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
        cmds, recs = self.__subcommands__()

        d = self.model_dump()
        if not isinstance(d, dict):
            d = dict(root=d)

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


class Command(_Command):
    model_config = ConfigDict(
        extra="forbid", arbitrary_types_allowed=True, populate_by_name=True
    )


class RootCommand(RootModel, _Command):
    """Generic Command for root models."""

    root: Any
    model_config = ConfigDict(arbitrary_types_allowed=True, populate_by_name=True)


class FlatCommand(Command):
    """Only used to discriminate Commands that should not be nested.

    HRUGroup, ReadFromNetCDF, Reservoir are examples of FlatCommand.
    """


class LineCommand(FlatCommand):
    r"""
    A non-nested Command on a single line.

    :CommandName {field_1} {field_2} ... {field_n}\n

    EvaluationPeriod is a FlatCommand.
    """

    def to_rv(self):
        cls = self.__class__
        out = [f":{cls.__name__:<20}"]
        for field in cls.model_fields.keys():
            out.append(str(getattr(self, field)))  # noqa: PERF401

        return " ".join(out) + "\n"

    @classmethod
    def parse(cls, s):
        """Parse the command and return an instance of LineCommand."""
        pat = rf":{cls.__name__}\s+(?P<args>.+)"
        fields = cls.model_fields.keys()
        if match := re.match(pat, s):
            args = match.group("args").split()
            return cls(**dict(zip(fields, args)))


class ListCommand(RootModel, _Command):
    """Use so that commands with __root__: Sequence[Record] behave like a list."""

    root: Sequence[Any]

    def __iter__(self):
        return iter(self.root)

    def __getitem__(self, item):
        return self.root[item]

    def __len__(self):
        return len(self.root)

    model_config = ConfigDict(arbitrary_types_allowed=True, populate_by_name=True)


class ParameterList(Record):
    name: str = ""
    values: Sequence[Union[Sym, str, None]] = ()

    @model_validator(mode="before")
    @classmethod
    def no_none_in_default(cls, data):
        """Make sure that no values are None for the [DEFAULT] record."""
        if data["name"] == "[DEFAULT]" and None in data["values"]:
            raise ValueError("Default record can not contain None.")
        return data

    def __str__(self):
        fmt = "{name:<16}" + len(self.values) * ",{:>18}"
        evals = []
        for v in self.values:
            ev = "_DEFAULT" if v is None else v
            evals.append(ev)

        return fmt.format(name=self.name, *evals)


class GenericParameterList(Command):
    parameters: Sequence[str] = Field(alias="Parameters", description="Parameter names")
    units: Optional[Sequence[str]] = Field(None, alias="Units")
    pl: Sequence[ParameterList]

    @model_validator(mode="after")
    def num_values_equal_num_names(self):
        """Check that the length of the parameter list equals the number of given parameter names."""
        n = len(self.parameters)

        for pl in self.pl:
            if len(pl.values) != n:
                raise ValueError("Number of values should match number of parameters.")

        return self

    @model_validator(mode="after")
    def set_default_units(self):
        n = len(self.parameters)
        if self.units is None:
            self.units = ("none",) * n
        return self


class RV(Command):
    """Base class for RV configuration objects."""

    __rv__ = ""

    @property
    def _template(self):
        return "{_commands}\n"

    @property
    def _indent(self):
        return ""

    model_config = ConfigDict(
        populate_by_name=True, validate_default=True, validate_assignment=True
    )


def parse_symbolic(value, **kwds):
    """
    Inject values of symbolic variables into object and return object.

    Note that parsing the output of `model_dump` can cause problems because there is not always enough information in the
    dictionary to recreate the correct model.
    """
    from pymbolic.mapper.evaluator import EvaluationMapper
    from pymbolic.primitives import ExpressionNode, Variable

    if isinstance(value, dict):
        return {k: parse_symbolic(v, **kwds) for k, v in value.items()}

    elif isinstance(value, (list, tuple)):
        return type(value)(parse_symbolic(v, **kwds) for v in value)

    elif isinstance(value, (_Record, _Command)):
        attrs = value.__dict__
        if isinstance(value, RootModel):
            attrs = attrs["root"]
        return value.model_validate(parse_symbolic(attrs, **kwds))

    elif isinstance(value, (Variable, ExpressionNode)):
        # Inject numerical values numerical value
        return EvaluationMapper(context=kwds)(value)

    else:
        return value
