"""
Raven model options that require non-trivial formatting for exportation to RV files.
"""

from dataclasses import dataclass
from typing import Any


REGISTRY = dict()


def register(rv, name=None):
    """Register a Raven option."""
    def wrapper(cls):
        REGISTRY[name] = cls
        return cls
    return wrapper


@dataclass
class Option:
    value: Any

    def to_rv(self):
        return self.value

    def __str__(self):
        return self.to_rv()

# Example of simple option
@register("suppress_output")
class SuppressOutput(Option):
    value: bool

    def to_rv(self):
        tag = ":SuppressOutput\n:DontWriteWatershedStorage"
        return tag if self.value else ""


