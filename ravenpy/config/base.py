from abc import ABC, abstractmethod
from enum import Enum
from typing import Any, OrderedDict, Sequence, Union

from pydantic.dataclasses import dataclass


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


@dataclass
class RavenCoefficient(RavenValue):
    value: float = None


@dataclass
class RavenSwitch(RavenCommand):
    value: bool = False

    def to_rv(self):
        return f":{self.__class__.__name__}\n" if self.value is True else ""


@dataclass
class RavenMapping(RavenCommand):
    value: OrderedDict[str, Any]

    def to_rv(self):
        return "\n".join(
            [f":{self.__class__.__name__} {k} {v}" for (k, v) in self.value]
        )
