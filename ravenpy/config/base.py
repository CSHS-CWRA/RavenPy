from abc import ABC, abstractmethod
from dataclasses import asdict, fields
from enum import Enum
from typing import Any, OrderedDict, Sequence, Union

from pydantic.dataclasses import dataclass
from pymbolic.mapper.evaluator import EvaluationMapper as EM
from pymbolic.primitives import Expression, Variable

Sym = Union[Variable, Expression, float, None]


class SymConfig:
    arbitrary_types_allowed = True


class RavenCommand(ABC):
    """Base class used for all Raven commands that are dataclasses which must implement some specialized rendering logic."""

    @abstractmethod
    def to_rv(self):
        pass

    def __str__(self):
        return self.to_rv()


@dataclass
class RavenOption(RavenCommand):
    option: Union[Enum, str]

    def to_rv(self):
        return f":{self.__class__.__name__:<20} {self.option.value}\n"


@dataclass
class RavenOptionList(RavenCommand):
    options: Union[Sequence[Enum], Sequence[str]]

    def to_rv(self):
        options = ", ".join([str(option.value) for option in self.options])
        return f":{self.__class__.__name__:<20} {options}\n"


@dataclass
class RavenValue(RavenCommand):
    value: str = None

    def to_rv(self):
        if self.value is not None:
            return f":{self.__class__.__name__:<20} {self.value}\n"
        return ""


@dataclass(config=SymConfig)
class RavenCoefficient(RavenValue):
    value: Sym = None


@dataclass
class RavenSwitch(RavenCommand):
    value: bool = False

    def to_rv(self):
        return f":{self.__class__.__name__}\n" if self.value is True else ""


@dataclass
class RavenMapping(RavenCommand):
    value: OrderedDict[str, Any]

    def to_rv(self):
        return (
            "\n".join(
                [f":{self.__class__.__name__} {k} {v}" for (k, v) in self.value.items()]
            )
            + "\n"
        )


@dataclass(config=SymConfig)
class Params:
    """Model parameters.

    Define type as `Sym` to support symbolic expressions.
    """

    def to_rv(self):
        return ""


def parse_symbolic(value, **kwds):
    """Inject values of symbolic variables into object and return object."""

    if isinstance(value, dict):
        return {k: parse_symbolic(v, **kwds) for k, v in value.items()}

    elif isinstance(value, (list, tuple)):
        return [parse_symbolic(v, **kwds) for v in value]

    elif isinstance(value, RavenCommand):
        # Cannot use asdict here as it recurses down nested classes
        attrs = {f.name: value.__dict__[f.name] for f in fields(value)}
        return value.__class__(**parse_symbolic(attrs, **kwds))

    elif isinstance(value, (Variable, Expression)):
        # Inject numerical values numerical value
        return EM(context=kwds)(value)

    else:
        return value
