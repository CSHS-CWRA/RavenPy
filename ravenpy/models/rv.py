import collections
import datetime as dt
from collections import namedtuple
from dataclasses import dataclass
from pathlib import Path
from textwrap import dedent
from typing import Dict, Tuple

import cftime
import six

from ravenpy.config.commands import (
    BasinIndexCommand,
    BasinStateVariablesCommand,
    ChannelProfileCommand,
    DataCommand,
    GaugeCommand,
    GriddedForcingCommand,
    HRUsCommand,
    HRUStateVariableTableCommand,
    LandUseClassesCommand,
    RavenConfig,
    ReservoirCommand,
    RoutingCommand,
    SBGroupPropertyMultiplierCommand,
    SoilClassesCommand,
    SoilProfilesCommand,
    StationForcingCommand,
    SubBasinGroupCommand,
    SubBasinsCommand,
    VegetationClassesCommand,
)

HRU = HRUsCommand.Record
HRUState = HRUStateVariableTableCommand.Record
LU = LandUseClassesCommand.Record
Sub = SubBasinsCommand.Record

"""
Raven configuration
-------------------

The RV class is used to store Raven parameters for emulated models.

Each model should subclass RV to define the parameters it expects using a namedtuple class. For example::

    class MyModel(RV):
        params = namedtuple('ModelParams', 'x1, x2, x3')
        init = namedtuple('ModelInit', 'i1, i2')
        hru = namedtuple('ModelHRU', 'hru1', hru2')

It can then be instantiated by passing values that will set as default values for each parameter::

    rv = MyModel(params=MyModel.params(1,2,3), init=MyModel.init(0,0), hru=MyModel.hru(4,5), name='basin')

values can then be modified either using attributes or properties::

    rv.name = 'LacVert'
    rv['evaluation_metrics'] = 'LOG_NASH'


Simulation end date and duration are updated automatically when duration, start date or end date are changed.

"""


class RVFile:
    def __init__(self, fn):
        """Read the content."""
        fn = Path(fn)

        self.stem = fn.with_suffix("").with_suffix("").stem
        self.suffixes = "".join(fn.suffixes)

        self.ext = ""
        self._store_ext(fn)

        # Whether extension indicates an Ostrich template file.
        self.is_tpl = fn.suffix in [".tpl", ".txt"]

        self.content = ""
        self.content = fn.read_text()

    def _store_ext(self, fn):
        try:
            self.ext = fn.suffixes[0][1:]
        except IndexError as e:
            msg = "\nFile {} does not look like a valid Raven/Ostrich config file.".format(
                fn
            )
            raise ValueError(msg) from e

    def rename(self, name):
        self.stem = name

    def write(self, path, **kwds):
        fn = (path / self.stem).with_suffix(self.suffixes)

        content = self.content
        if kwds:
            content = content.format(**kwds)

        fn.write_text(content)
        return fn

    @property
    def tags(self):
        """Return a list of tags within the templates."""
        import re

        pattern = re.compile(r"{([\.\w]+)}")

        return pattern.findall(self.content)


class RV(collections.abc.Mapping, RavenConfig):
    """Generic configuration class.

    RV provides two mechanisms to set values, a dictionary-like interface and an object-like interface::

        rv = RV(a=None)
        rv['a'] = 1
        rv.a = 2

    The dictionary like interface only allows the modification of values for existing items, while the object interface
    allows the creation of new attributes::

      rv['c'] = 1

    will raise an AttributeError, while::

      rv.c = 1

    will create a new `c` attribute and assign it the value 1.

    """

    def __init__(self, **kwargs):
        # Set initial default values
        for key, val in kwargs.items():
            setattr(self, key, val)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        if not hasattr(self, key):
            raise AttributeError("Trying to assign unrecognized object: {}".format(key))

        setattr(self, key, value)

    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        return iter(self.keys())

    def keys(self):
        # Attributes
        a = list(filter(lambda x: not x.startswith("_"), self.__dict__))

        # Properties
        p = list(
            filter(
                lambda x: isinstance(getattr(self.__class__, x, None), property),
                dir(self),
            )
        )
        return a + p

    def items(self):
        for attribute in self.keys():
            yield attribute, getattr(self, attribute)

    def update(self, items, force=False):
        """Update values from dictionary items.

        Parameters
        ----------
        items : dict
          Dictionary of values.
        force : bool
          If True, un-initialized keys can be set.
        """
        if force:
            for key, val in items.items():
                setattr(self, key, val)
        else:
            for key, val in items.items():
                self[key] = val


# class Ost(RV):
#     def __init__(self, **kwargs):
#         self._max_iterations = None
#         self._random_seed = None

#         super(Ost, self).__init__(**kwargs)

#     @property
#     def max_iterations(self):
#         return self._max_iterations

#     @max_iterations.setter
#     def max_iterations(self, x):
#         if x < 1:
#             raise ValueError("Max iteration should be a positive integer: {}".format(x))
#         else:
#             self._max_iterations = x

#     @property
#     def random_seed(self):
#         if self._random_seed is not None:
#             return "RandomSeed {}".format(self._random_seed)
#         return ""

#     @random_seed.setter
#     def random_seed(self, value):
#         if value >= 0:
#             self._random_seed = value
#         else:
#             self._random_seed = None


# def isinstance_namedtuple(x):
#     a = isinstance(x, tuple)
#     b = getattr(x, "_fields", None) is not None
#     return a and b


def guess_linear_transform(actual, expected):
    """Return RVT compatible dictionary for variable unit transformations.

    Parameters
    ----------
    actual : dict
      The units of each variable.
    expected : dict
      The units expected by Raven.

    Returns
    -------
    dict
      Dictionary keyed by <variable_name>_linear_transform, storing "<scale> <offset>"
      strings used by Raven to transform units.

    """
    # TODO : For precip we also need the frequency to sum over one day.
