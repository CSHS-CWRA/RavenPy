import collections
import datetime as dt
from collections import namedtuple
from pathlib import Path
from textwrap import dedent
from typing import Dict, List, Tuple

import cftime
import six
import xarray as xr
from dataclasses import dataclass, replace
from xclim.core.units import units2pint

from .commands import (
    BasinIndexCommand,
    BasinStateVariablesCommand,
    ChannelProfileCommand,
    DataCommand,
    GaugeCommand,
    GriddedForcingCommand,
    GridWeightsCommand,
    HRUsCommand,
    HRUStateVariableTableCommand,
    LandUseClassesCommand,
    MonthlyAverageCommand,
    ObservationDataCommand,
    RainCorrection,
    RavenConfig,
    ReservoirCommand,
    Routing,
    SBGroupPropertyMultiplierCommand,
    SnowCorrection,
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

# CF standard names used for mappings
default_input_variables = (
    "pr",
    "rainfall",
    "prsn",
    "tasmin",
    "tasmax",
    "tas",
    "evspsbl",
    "water_volume_transport_in_river_channel",
)

# Map CF-Convention standard name to Raven Forcing name
forcing_names = {
    "tasmin": "TEMP_MIN",
    "tasmax": "TEMP_MAX",
    "tas": "TEMP_AVE",
    "rainfall": "RAINFALL",
    "pr": "PRECIP",
    "prsn": "SNOWFALL",
    "evspsbl": "PET",
    "water_volume_transport_in_river_channel": "HYDROGRAPH",
}


# Alternate typical variable names found in netCDF files, keyed by CF standard name
alternate_nc_names = {
    "tasmin": ["tasmin", "tmin"],
    "tasmax": ["tasmax", "tmax"],
    "tas": ["tas", "t2m"],
    "rainfall": ["rainfall", "rain"],
    "pr": ["pr", "precip", "prec", "precipitation", "tp"],
    "prsn": ["prsn", "snow", "snowfall", "solid_precip"],
    "evspsbl": ["pet", "evap", "evapotranspiration"],
    "water_volume_transport_in_river_channel": [
        "qobs",
        "discharge",
        "streamflow",
        "dis",
    ],
}

rain_snow_fraction_options = (
    "RAINSNOW_DATA",
    "RAINSNOW_DINGMAN",
    "RAINSNOW_UBC",
    "RAINSNOW_HBV",
    "RAINSNOW_HARDER",
    "RAINSNOW_HSPF",
)

evaporation_options = (
    "PET_CONSTANT",
    "PET_PENMAN_MONTEITH",
    "PET_PENMAN_COMBINATION",
    "PET_PRIESTLEY_TAYLOR",
    "PET_HARGREAVES",
    "PET_HARGREAVES_1985",
    "PET_FROMMONTHLY",
    "PET_DATA",
    "PET_HAMON_1961",
    "PET_TURC_1961",
    "PET_MAKKINK_1957",
    "PET_MONTHLY_FACTOR",
    "PET_MOHYSE",
    "PET_OUDIN",
)

calendar_options = (
    "PROLEPTIC_GREGORIAN",
    "JULIAN",
    "GREGORIAN",
    "STANDARD",
    "NOLEAP",
    "365_DAY",
    "ALL_LEAP",
    "366_DAY",
)


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


class RVT(RV):
    def __init__(self, **kwargs):
        self.pr = {}
        self.rainfall = {}
        self.prsn = {}
        self.tasmin = {}
        self.tasmax = {}
        self.tas = {}
        self.evspsbl = {}
        self.water_volume_transport_in_river_channel = {}

        self.nc_index = None
        self.gauge_latitude = None
        self.gauge_longitude = None
        self.gauge_elevation = None

        self.raincorrection = 1
        self.snowcorrection = 1

        self.monthly_ave_evaporation = ()
        self.monthly_ave_temperature = ()

        self.gridded_forcings = ()

        # For a distributed model these weights will be shared among all the StationForcing commands
        self.grid_weights = None

        self.var_cmds = {}
        # Dictionary of potential variable names, keyed by CF standard name.
        # http://cfconventions.org/Data/cf-standard-names/60/build/cf-standard-name-table.html
        # PET is the potential evapotranspiration, while evspsbl is the actual evap.
        # TODO: Check we're not mixing precip and rainfall.

        super(RVT, self).__init__(**kwargs)

    @property
    def variables(self):
        return (getattr(self, name) for name in default_input_variables)

    @property
    def gauge(self):
        data = [o for o in self.var_cmds.values() if isinstance(o, DataCommand)]
        if data:
            return GaugeCommand(
                latitude=self.gauge_latitude,
                longitude=self.gauge_longitude,
                elevation=self.gauge_elevation,
                raincorrection=self.raincorrection,
                snowcorrection=self.snowcorrection,
                monthly_ave_evaporation=self.monthly_ave_evaporation,
                monthly_ave_temperature=self.monthly_ave_temperature,
                data=data,
            )
        else:
            return ""

    @property
    def station_forcing_list(self):
        data = [
            o for o in self.var_cmds.values() if isinstance(o, StationForcingCommand)
        ]
        return "\n\n".join(map(str, data))

    @property
    def gridded_forcing_list(self):
        data = [
            o for o in self.var_cmds.values() if isinstance(o, GriddedForcingCommand)
        ]
        # This is really a hack for now, as model.rvt.gridded_forcings are set
        # directly by the user in TestRouting.test_lievre_tutorial
        data += self.gridded_forcings
        return "\n\n".join(map(str, data))

    @property
    def observed_data_cmd(self):
        return self.var_cmds.get("water_volume_transport_in_river_channel", "")

    @property
    def raincorrection_cmd(self):
        return RainCorrection(self.raincorrection)

    @property
    def snowcorrection_cmd(self):
        return SnowCorrection(self.snowcorrection)


class RVI(RV):
    def __init__(self, **kwargs):
        self.name = None
        self.area = None
        self.elevation = None
        self.latitude = None
        self.longitude = None
        self.run_index = 0
        self.raven_version = "3.0.1 rev#275"
        self.routing = "ROUTE_NONE"

        self._run_name = "run"
        self._start_date = None
        self._end_date = None
        self._now = None
        self._rain_snow_fraction = "RAINSNOW_DATA"
        self._evaporation = None
        self._ow_evaporation = None
        self._duration = 1
        self._time_step = 1.0
        self._evaluation_metrics = "NASH_SUTCLIFFE RMSE"
        self._suppress_output = False
        self._calendar = "standard"

        super(RVI, self).__init__(**kwargs)

    @property
    def run_name(self):
        return self._run_name

    @run_name.setter
    def run_name(self, x):
        if isinstance(x, six.string_types):
            self._run_name = x
        else:
            raise ValueError("Must be string")

    @property
    def start_date(self):
        return self._start_date

    @start_date.setter
    def start_date(self, x):
        if isinstance(x, dt.datetime):
            self._start_date = self._dt2cf(x)
        else:
            raise ValueError("Must be datetime")

        if x == dt.datetime(1, 1, 1):
            return

        if self._duration is None:
            self._update_duration()
        else:
            self._update_end_date()

    @property
    def end_date(self):
        return self._end_date

    @end_date.setter
    def end_date(self, x):
        if isinstance(x, dt.datetime):
            self._end_date = self._dt2cf(x)
        else:
            raise ValueError("Must be datetime")

        if x != dt.datetime(1, 1, 1):
            self._update_duration()

    @property
    def duration(self):
        return self._duration

    @duration.setter
    def duration(self, x):
        if isinstance(x, int):
            if x > 0:
                self._duration = x
        else:
            raise ValueError("Must be int")

        if x > 0:
            self._update_end_date()

    @property
    def time_step(self):
        return self._time_step

    @time_step.setter
    def time_step(self, x):
        self._time_step = x

    @property
    def evaluation_metrics(self):
        return self._evaluation_metrics

    @evaluation_metrics.setter
    def evaluation_metrics(self, x):
        if not isinstance(x, six.string_types):
            raise ValueError("Evaluation metrics must be string.")

        for metric in x.split():
            if metric not in {
                "NASH_SUTCLIFFE",
                "LOG_NASH",
                "RMSE",
                "PCT_BIAS",
                "ABSERR",
                "ABSMAX",
                "PDIFF",
                "TMVOL",
                "RCOEFF",
                "NSC",
                "KLING_GUPTA",
            }:
                raise ValueError("{} is not a metric recognized by Raven.")

        self._evaluation_metrics = x

    def _update_duration(self):
        if self.end_date is not None and self.start_date is not None:
            self._duration = (self.end_date - self.start_date).days

    def _update_end_date(self):
        if self.start_date is not None and self.duration is not None:
            self._end_date = self.start_date + dt.timedelta(days=self.duration)

    @property
    def now(self):
        return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    @property
    def suppress_output(self):
        tag = ":SuppressOutput\n:DontWriteWatershedStorage"
        return tag if self._suppress_output else ""

    @suppress_output.setter
    def suppress_output(self, value):
        if not isinstance(value, bool):
            raise ValueError
        self._suppress_output = value

    @property
    def routing_cmd(self):
        return Routing(value=self.routing)

    @property
    def rain_snow_fraction(self):
        """Rain snow partitioning."""
        return self._rain_snow_fraction

    @rain_snow_fraction.setter
    def rain_snow_fraction(self, value):
        """Can be one of

        - RAINSNOW_DATA
        - RAINSNOW_DINGMAN
        - RAINSNOW_UBC
        - RAINSNOW_HBV
        - RAINSNOW_HARDER
        - RAINSNOW_HSPF
        """
        v = value.upper()

        if v in rain_snow_fraction_options:
            self._rain_snow_fraction = v
        else:
            raise ValueError(f"Value should be one of {rain_snow_fraction_options}.")

    @property
    def evaporation(self):
        """Evaporation scheme"""
        return self._evaporation

    @evaporation.setter
    def evaporation(self, value):
        v = value.upper()
        if v in evaporation_options:
            self._evaporation = v
        else:
            raise ValueError(f"Value {v} should be one of {evaporation_options}.")

    @property
    def ow_evaporation(self):
        """Open-water evaporation scheme"""
        return self._ow_evaporation

    @ow_evaporation.setter
    def ow_evaporation(self, value):
        v = value.upper()
        if v in evaporation_options:
            self._ow_evaporation = v
        else:
            raise ValueError(f"Value {v} should be one of {evaporation_options}.")

    @property
    def calendar(self):
        """Calendar"""
        return self._calendar.upper()

    @calendar.setter
    def calendar(self, value):
        if value.upper() in calendar_options:
            self._calendar = value
        else:
            raise ValueError(f"Value should be one of {calendar_options}.")

    def _dt2cf(self, date):
        """Convert datetime to cftime datetime."""
        return cftime._cftime.DATE_TYPES[self._calendar.lower()](*date.timetuple()[:6])


class RVC(RV):
    def __init__(
        self,
        hru_states: Dict[int, HRUState] = None,
        basin_states: Dict[int, BasinIndexCommand] = None,
        **kwds,
    ):
        self.hru_states = hru_states or {}
        self.basin_states = basin_states or {}
        super().__init__(**kwds)

    def parse(self, rvc):
        """Set initial conditions based on *solution* output file.

        Parameters
        ----------
        path : string
          `solution.rvc` content.
        """
        objs = parse_solution(rvc)
        self.hru_states, self.basin_states = get_states(objs)

    @property
    def hru_state(self):
        return self.hru_states.get(1, None)

    @hru_state.setter
    def hru_state(self, value):
        self.hru_states[1] = value

    @property
    def basin_state(self):
        return self.basin_states.get(1, None)

    @basin_state.setter
    def basin_state(self, value):
        self.basin_states[1] = value

    @property
    def hru_states_cmd(self):
        """Return HRU state values."""
        return HRUStateVariableTableCommand(self.hru_states)

    @property
    def basin_states_cmd(self):
        """Return basin state variables."""
        return BasinStateVariablesCommand(self.basin_states)


@dataclass
class RVH(RV):
    subbasins: Tuple[SubBasinsCommand.Record] = (SubBasinsCommand.Record(),)
    land_subbasins: Tuple[int] = ()
    land_subbasin_property_multiplier: SBGroupPropertyMultiplierCommand = ""
    lake_subbasins: Tuple[int] = ()
    lake_subbasin_property_multiplier: SBGroupPropertyMultiplierCommand = ""
    reservoirs: Tuple[ReservoirCommand] = ()
    hrus: Tuple[HRUsCommand.Record] = ()

    template = """
    {subbasins_cmd}

    {hrus_cmd}

    {land_subbasin_group_cmd}

    {lake_subbasin_group_cmd}

    {reservoir_cmd_list}
    """

    @property
    def subbasins_cmd(self):
        return SubBasinsCommand(self.subbasins)

    @property
    def land_subbasin_group_cmd(self):
        return SubBasinGroupCommand("Land", self.land_subbasins)

    @property
    def lake_subbasin_group_cmd(self):
        return SubBasinGroupCommand("Lakes", self.lake_subbasins)

    @property
    def reservoir_cmd_list(self):
        return "\n\n".join(map(str, self.reservoirs))

    @property
    def hrus_cmd(self):
        return HRUsCommand(self.hrus)

    def to_rv(self):
        # params = self.items()
        return dedent(self.template).format(**dict(self.items()))


@dataclass
class RVP(RV):
    params: Tuple[namedtuple] = ()
    soil_classes: Tuple[SoilClassesCommand.Record] = ()
    soil_profiles: Tuple[SoilProfilesCommand.Record] = ()
    vegetation_classes: Tuple[VegetationClassesCommand.Record] = ()
    land_use_classes: Tuple[LandUseClassesCommand.Record] = ()
    channel_profiles: Tuple[ChannelProfileCommand] = ()
    avg_annual_runoff: float = 0

    @property
    def soil_classes_cmd(self):
        return SoilClassesCommand(self.soil_classes)

    @property
    def soil_profiles_cmd(self):
        return SoilProfilesCommand(self.soil_profiles)

    @property
    def vegetation_classes_cmd(self):
        return VegetationClassesCommand(self.vegetation_classes)

    @property
    def land_use_classes_cmd(self):
        return LandUseClassesCommand(self.land_use_classes)

    @property
    def channel_profile_cmd_list(self):
        return "\n\n".join(map(str, self.channel_profiles))

    # Note sure about this!
    # def to_rv(self):
    #     params = self.items()
    #     return dedent(self.template).format(**dict(self.items()))


class Ost(RV):
    def __init__(self, **kwargs):
        self._max_iterations = None
        self._random_seed = None

        super(Ost, self).__init__(**kwargs)

    @property
    def max_iterations(self):
        return self._max_iterations

    @max_iterations.setter
    def max_iterations(self, x):
        if x < 1:
            raise ValueError("Max iteration should be a positive integer: {}".format(x))
        else:
            self._max_iterations = x

    @property
    def random_seed(self):
        if self._random_seed is not None:
            return "RandomSeed {}".format(self._random_seed)
        return ""

    @random_seed.setter
    def random_seed(self, value):
        if value >= 0:
            self._random_seed = value
        else:
            self._random_seed = None


def isinstance_namedtuple(x):
    a = isinstance(x, tuple)
    b = getattr(x, "_fields", None) is not None
    return a and b


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


def get_states(solution, hru_index=None, basin_index=None):
    """Return state variables.

    Parameters
    ----------
    solution : dict
      `solution.rvc` parsed content.
    """
    hru_state = {}
    basin_state = {}

    for index, params in solution["HRUStateVariableTable"]["data"].items():
        hru_state[index] = HRUState(*params)

    for index, raw in solution["BasinStateVariables"]["BasinIndex"].items():
        params = {k.lower(): v for (k, v) in raw.items()}
        basin_state[index] = BasinIndexCommand(**params)

    if hru_index is not None:
        hru_state = hru_state[hru_index]

    if basin_index is not None:
        basin_state = basin_state[basin_index]

    return hru_state, basin_state


def parse_solution(rvc):
    """Parse solution file and return dictionary of parameters that can then be used to reinitialize the model.

    Parameters
    ----------
    rvc : str
      Content of a solution.rvc file.
    """
    # Create a generator that will consume lines one by one.
    # tags = ['{i.' + a.lower().replace('[', '').replace(']', '') + '}' for a in atts]
    lines = iter(rvc.splitlines())
    return _parser(lines)


def _parser(lines, indent="", fmt=str):
    import itertools
    import re

    header_pat = re.compile(r"(\s*):(\w+)\s?,?\s*(.*)")

    out = collections.defaultdict(dict)
    old_key = None
    for line in lines:
        header = header_pat.match(line)
        if header:
            new_indent, key, value = header.groups()
            if new_indent > indent:
                out[old_key] = _parser(itertools.chain([line], lines), new_indent)
            elif new_indent < indent:
                return out
            else:
                if key == "BasinIndex":
                    i, name = value.split(",")
                    i = int(i)
                    out[key][i] = dict(
                        index=i, name=name, **_parser(lines, new_indent + "  ", float)
                    )
                # elif key in ["Qlat", "Qout"]:
                #     n, *values, last = value.split(",")
                #     out[key] = list(map(float, values))
                #     out[key + "Last"] = float(last)
                elif key in ["Qin", "Qout", "Qlat"]:
                    n, *values = value.split(",")
                    out[key] = (int(n),) + tuple(map(float, values))
                else:
                    out[key] = (
                        list(map(fmt, value.split(","))) if "," in value else fmt(value)
                    )

            old_key = key
        else:
            data = line.split(",")
            i = int(data.pop(0))
            out["data"][i] = [i] + list(map(float, data))

    return out
