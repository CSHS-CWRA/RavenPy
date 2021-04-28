import collections
import datetime as dt
from abc import ABC, abstractmethod
from dataclasses import is_dataclass, replace
from enum import Enum
from pathlib import Path
from textwrap import dedent
from typing import Dict, Optional, Tuple

import cf_xarray
import cftime
import numpy as np
import xarray as xr

from ravenpy.config.commands import (
    BasinIndexCommand,
    BasinStateVariablesCommand,
    ChannelProfileCommand,
    DataCommand,
    GaugeCommand,
    GriddedForcingCommand,
    GridWeightsCommand,
    HRUsCommand,
    HRUState,
    HRUStateVariableTableCommand,
    LandUseClassesCommand,
    ObservationDataCommand,
    SoilClassesCommand,
    SoilProfilesCommand,
    StationForcingCommand,
    SubBasinGroupCommand,
    SubBasinsCommand,
    VegetationClassesCommand,
)


class RV(ABC):
    def __init__(self, config, **kwds):
        # Each RV has a reference to their parent object in order to access sibling RVs.
        self._config = config
        self.is_ostrich_tmpl = False
        self.content = None

    def update(self, key, value):
        if hasattr(self, key):
            # Special case: the user might be trying to update `params` or `derived_params`,
            # (which are dataclasses) by supplying a list of values (either a list, tuple of
            # numpy array); if that's the case, cast those values into a new instance of the
            # corresponding dataclass.
            attr = getattr(self, key)
            if is_dataclass(attr) and isinstance(value, (list, tuple, np.ndarray)):
                value = attr.__class__(*value)
            setattr(self, key, value)
            return True
        return False

    def get_extra_attributes(self, d):
        """
        Not sure about this: for the moment I use it only for certain params that must
        be injected in the RVH by the MOHYSE emulator. The idea is to complete the `d`
        dict used in the `to_rv` method for the template with extra attributes that have
        been added in the emulator.
        """
        e = {}
        for k, v in self.__dict__.items():
            if k not in d:
                e[k] = v
        return e

    def set_tmpl(self, tmpl, is_ostrich=False):
        self.tmpl = tmpl
        self.is_ostrich_tmpl = is_ostrich

    @property
    @abstractmethod
    def tmpl(self):
        pass

    @abstractmethod
    def to_rv(self) -> str:
        pass


#########
# R V C #
#########


class RVC(RV):

    tmpl = """
    {hru_states}

    {basin_states}
    """

    def __init__(self, config):
        super().__init__(config)
        self.hru_states: Dict[int, HRUState] = {}
        self.basin_states: Dict[int, BasinIndexCommand] = {}

    def reset(self, **kwargs):
        self.hru_states = {}
        self.basin_states = {}

    def set_hru_state(self, hru_state: HRUState):
        self.hru_states[hru_state.index] = hru_state

    def set_basin_state(self, basin_state: BasinIndexCommand):
        self.basin_states[basin_state.index] = basin_state

    @classmethod
    def create_solution(cls, solution_str):
        rvc = RVC(None)
        rvc.parse_solution(solution_str)
        return rvc

    def parse_solution(self, solution_str):
        self.hru_states = HRUStateVariableTableCommand.parse(solution_str).hru_states
        self.basin_states = BasinStateVariablesCommand.parse(solution_str).basin_states

    def to_rv(self):
        d = {
            "hru_states": HRUStateVariableTableCommand(self.hru_states),
            "basin_states": BasinStateVariablesCommand(self.basin_states),
        }

        d.update(self.get_extra_attributes(d))

        return dedent(self.tmpl).format(**d)


#########
# R V H #
#########


class RVH(RV):

    tmpl = """
    {subbasins}

    {hrus}

    {land_subbasin_group}

    {land_subbasin_property_multiplier}

    {lake_subbasin_group}

    {lake_subbasin_property_multiplier}

    {reservoirs}
    """

    def __init__(self, config):
        super().__init__(config)
        self.hrus = ()
        self.subbasins = ()
        self.land_subbasin_ids = ()
        self.land_subbasin_property_multiplier = None
        self.lake_subbasin_ids = ()
        self.lake_subbasin_property_multiplier = None
        self.reservoirs = ()

    def to_rv(self):
        d = {
            "subbasins": SubBasinsCommand(self.subbasins),
            "hrus": HRUsCommand(self.hrus),
            "land_subbasin_group": SubBasinGroupCommand("Land", self.land_subbasin_ids),
            "land_subbasin_property_multiplier": self.land_subbasin_property_multiplier
            or "",
            "lake_subbasin_group": SubBasinGroupCommand(
                "Lakes", self.lake_subbasin_ids
            ),
            "lake_subbasin_property_multiplier": self.lake_subbasin_property_multiplier
            or "",
            "reservoirs": "\n\n".join(map(str, self.reservoirs)),
        }

        d.update(self.get_extra_attributes(d))

        return dedent(self.tmpl).format(**d)


#########
# R V I #
#########


class RVI(RV):

    _pre_tmpl = """
    :Calendar              {calendar}
    :RunName               {run_name}-{run_index}
    :StartDate             {start_date}
    :EndDate               {end_date}
    :TimeStep              {time_step}
    :Method                ORDERED_SERIES
    """

    tmpl = """
    """

    # This part must be at the end of the file in particular because `:EvaluationMetrics` must
    # come after `:SoilModel`
    _post_tmpl = """
    :EvaluationMetrics     {evaluation_metrics}
    :WriteNetcdfFormat     yes
    #:WriteForcingFunctions
    :SilentMode
    :PavicsMode
    {suppress_output}

    :NetCDFAttribute title Simulated river discharge
    :NetCDFAttribute history Created on {now} by Raven
    :NetCDFAttribute references  Craig, J.R., and the Raven Development Team, Raven user's and developer's manual (Version 2.8), URL: http://raven.uwaterloo.ca/ (2018).
    :NetCDFAttribute comment Raven Hydrological Framework version {raven_version}
    :NetCDFAttribute model_id {identifier}
    :NetCDFAttribute time_frequency day
    :NetCDFAttribute time_coverage_start {start_date}
    :NetCDFAttribute time_coverage_end {end_date}
    """

    class EvaporationOptions(Enum):
        PET_CONSTANT = "PET_CONSTANT"
        PET_PENMAN_MONTEITH = "PET_PENMAN_MONTEITH"
        PET_PENMAN_COMBINATION = "PET_PENMAN_COMBINATION"
        PET_PRIESTLEY_TAYLOR = "PET_PRIESTLEY_TAYLOR"
        PET_HARGREAVES = "PET_HARGREAVES"
        PET_HARGREAVES_1985 = "PET_HARGREAVES_1985"
        PET_FROMMONTHLY = "PET_FROMMONTHLY"
        PET_DATA = "PET_DATA"
        PET_HAMON_1961 = "PET_HAMON_1961"
        PET_TURC_1961 = "PET_TURC_1961"
        PET_MAKKINK_1957 = "PET_MAKKINK_1957"
        PET_MONTHLY_FACTOR = "PET_MONTHLY_FACTOR"
        PET_MOHYSE = "PET_MOHYSE"
        PET_OUDIN = "PET_OUDIN"

    class CalendarOptions(Enum):
        PROLEPTIC_GREGORIAN = "PROLEPTIC_GREGORIAN"
        JULIAN = "JULIAN"
        GREGORIAN = "GREGORIAN"
        STANDARD = "STANDARD"
        NOLEAP = "NOLEAP"
        _365_DAY = "365_DAY"
        ALL_LEAP = "ALL_LEAP"
        _366_DAY = "366_DAY"

    class EvaluationMetrics(Enum):
        NASH_SUTCLIFFE = "NASH_SUTCLIFFE"
        LOG_NASH = "LOG_NASH"
        RMSE = "RMSE"
        PCT_BIAS = "PCT_BIAS"
        ABSERR = "ABSERR"
        ABSMAX = "ABSMAX"
        PDIFF = "PDIFF"
        TMVOL = "TMVOL"
        RCOEFF = "RCOEFF"
        NSC = "NSC"
        KLING_GUPTA = "KLING_GUPTA"

    class RainSnowFractionOptions(Enum):
        DATA = "RAINSNOW_DATA"
        DINGMAN = "RAINSNOW_DINGMAN"
        UBC = "RAINSNOW_UBC"
        HBV = "RAINSNOW_HBV"
        HARDER = "RAINSNOW_HARDER"
        HSPF = "RAINSNOW_HSPF"

    class RoutingOptions(Enum):
        DIFFUSIVE_WAVE = "ROUTE_DIFFUSIVE_WAVE"
        HYDROLOGIC = "ROUTE_HYDROLOGIC"
        NONE = "ROUTE_NONE"
        STORAGE_COEFF = "ROUTE_STORAGE_COEFF"
        PLUG_FLOW = "ROUTE_PLUG_FLOW"
        MUSKINGUM = "MUSKINGUM"

    def __init__(self, config):
        super().__init__(config)

        # These are attributes that can be modified/set directly
        self.run_name = "run"
        self.run_index = 0
        self.raven_version = "3.X.X"
        self.time_step = 1.0

        # These correspond to properties whose setters will pass their value through
        # an Enum cast (triggering a ValueError at runtime for unknown values) and
        # getters will be used when rendering the template in the `to_rv` method
        self._calendar = RVI.CalendarOptions.STANDARD
        self._routing = RVI.RoutingOptions.NONE
        self._start_date = None
        self._end_date = None
        self._rain_snow_fraction = RVI.RainSnowFractionOptions.DATA
        self._evaporation = None
        self._ow_evaporation = None
        self._duration = 1
        self._evaluation_metrics = [
            RVI.EvaluationMetrics.NASH_SUTCLIFFE,
            RVI.EvaluationMetrics.RMSE,
        ]
        self._suppress_output = False

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
    def evaluation_metrics(self):
        if self._evaluation_metrics:
            return ",".join(m.value for m in self._evaluation_metrics)
        return None

    @evaluation_metrics.setter
    def evaluation_metrics(self, values):
        if not isinstance(values, (list, set, tuple)):
            values = [values]
        ms = []
        for v in values:
            v = v.upper() if isinstance(v, str) else v.value
            ms.append(RVI.EvaluationMetrics(v))
        self._evaluation_metrics = ms

    def _update_duration(self):
        if self.end_date is not None and self.start_date is not None:
            self._duration = (self.end_date - self.start_date).days

    def _update_end_date(self):
        if self.start_date is not None and self.duration is not None:
            self._end_date = self.start_date + dt.timedelta(days=self.duration)

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
    def routing(self):
        return self._routing.value

    @routing.setter
    def routing(self, value):
        v = value.upper() if isinstance(value, str) else value.value
        self._routing = RVI.RoutingOptions(v)

    @property
    def rain_snow_fraction(self):
        """Rain snow partitioning."""
        return self._rain_snow_fraction.value

    @rain_snow_fraction.setter
    def rain_snow_fraction(self, value):
        v = value.upper() if isinstance(value, str) else value.value
        self._rain_snow_fraction = RVI.RainSnowFractionOptions(v)

    @property
    def evaporation(self):
        """Evaporation scheme"""
        return self._evaporation.value if self._evaporation else None

    @evaporation.setter
    def evaporation(self, value):
        v = value.upper() if isinstance(value, str) else value.value
        self._evaporation = RVI.EvaporationOptions(v)

    @property
    def ow_evaporation(self):
        """Open-water evaporation scheme"""
        return self._ow_evaporation.value if self._ow_evaporation else None

    @ow_evaporation.setter
    def ow_evaporation(self, value):
        v = value.upper() if isinstance(value, str) else value.value
        self._ow_evaporation = RVI.EvaporationOptions(v)

    @property
    def calendar(self):
        """Calendar"""
        return self._calendar.value

    @calendar.setter
    def calendar(self, value):
        v = value.upper() if isinstance(value, str) else value.value
        self._calendar = RVI.CalendarOptions(v)

    def _dt2cf(self, date):
        """Convert datetime to cftime datetime."""
        return cftime._cftime.DATE_TYPES[self.calendar.lower()](*date.timetuple()[:6])

    def to_rv(self):

        # Attributes (not starting with "_")
        a = list(filter(lambda x: not x.startswith("_"), self.__dict__))

        # Properties (computing values corresponding to attributes starting with "_')
        p = list(
            filter(
                lambda x: isinstance(getattr(self.__class__, x, None), property),
                dir(self),
            )
        )

        d = {attr: getattr(self, attr) for attr in a + p}

        d["identifier"] = self._config.identifier
        d["now"] = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        t = dedent(self._pre_tmpl) + dedent(self.tmpl) + dedent(self._post_tmpl)
        return t.format(**d)


##########
# R V P  #
##########


class RVP(RV):

    # This is expected to be defined by the emulators.
    tmpl = """
    """

    def __init__(self, config):
        super().__init__(config)

        # Model specific params and derived params
        self.params = None
        self.derived_params = None

        self.soil_classes: Tuple[SoilClassesCommand.Record, ...] = ()
        self.soil_profiles: Tuple[SoilProfilesCommand.Record, ...] = ()
        self.vegetation_classes: Tuple[VegetationClassesCommand.Record, ...] = ()
        self.land_use_classes: Tuple[LandUseClassesCommand.Record, ...] = ()
        self.channel_profiles: Tuple[ChannelProfileCommand, ...] = ()
        self.avg_annual_runoff: Optional[float] = None

    def to_rv(self):
        d = {
            "params": self.params,
            "derived_params": self.derived_params,
            "soil_classes": SoilClassesCommand(self.soil_classes),
            "soil_profiles": SoilProfilesCommand(self.soil_profiles),
            "vegetation_classes": VegetationClassesCommand(self.vegetation_classes),
            "land_use_classes": LandUseClassesCommand(self.land_use_classes),
            "channel_profiles": "\n\n".join(map(str, self.channel_profiles)),
            "avg_annual_runoff": f":AvgAnnualRunoff {self.avg_annual_runoff}"
            if self.avg_annual_runoff
            else "",
        }
        return dedent(self.tmpl).format(**d)


#########
# R V T #
#########


class RVT(RV):

    tmpl = """
    {gauge}

    {forcing_list}

    {observed_data}
    """

    # Keys are standard names
    NC_VARS = {
        "tasmin": {"raven": "TEMP_MIN", "alts": ["tmin"]},
        "tasmax": {"raven": "TEMP_MAX", "alts": ["tmax"]},
        "tas": {"raven": "TEMP_AVE", "alts": ["t2m"]},
        "rainfall": {"raven": "RAINFALL", "alts": ["rain"]},
        "pr": {
            "raven": "PRECIP",
            "alts": ["precip", "prec", "precipitation", "tp"],
        },
        "prsn": {
            "raven": "SNOWFALL",
            "alts": ["snow", "snowfall", "solid_precip"],
        },
        "evspsbl": {"raven": "PET", "alts": ["pet", "evap", "evapotranspiration"]},
        "water_volume_transport_in_river_channel": {
            "raven": "HYDROGRAPH",
            "alts": [
                "qobs",
                "discharge",
                "streamflow",
                "dis",
            ],
        },
    }

    def __init__(self, config):
        super().__init__(config)

        self._var_cmds = {k: {} for k in RVT.NC_VARS.keys()}

        self.nc_index = 0
        self.grid_weights = None
        self.rain_correction = None
        self.snow_correction = None
        self.monthly_ave_evaporation = None
        self.monthly_ave_temperature = None

        self._nc_latitude = []
        self._nc_longitude = []
        self._nc_elevation = []
        self._number_grid_cells = 0

    def add_nc_variable(self, **kwargs):
        std_name = kwargs.get("name", kwargs["var_name_nc"])
        is_obs_var = kwargs.pop("is_observation", False)
        if len(kwargs["dim_names_nc"]) == 1:
            if std_name == "water_volume_transport_in_river_channel" or is_obs_var:
                cmd = ObservationDataCommand(**kwargs)
            else:
                cmd = DataCommand(**kwargs)
        elif len(kwargs["dim_names_nc"]) == 2:
            if std_name == "water_volume_transport_in_river_channel" or is_obs_var:
                cmd = ObservationDataCommand(**kwargs)
            else:
                cmd = StationForcingCommand(**kwargs)
        else:
            cmd = GriddedForcingCommand(**kwargs)

        if isinstance(self._var_cmds.get(std_name, None), dict):
            self._var_cmds[std_name] = replace(cmd, **self._var_cmds[std_name])
        else:
            self._var_cmds[std_name] = cmd

    def configure_from_nc_data(self, fns):

        # Important note: if the object at key `k` is a dict (as opposed to a `Command`),
        # don't reset it because it contains initial user-defined config (for the future Command
        # object at that particular key)
        for std_name in RVT.NC_VARS:
            v = self._var_cmds[std_name]
            if not isinstance(v, dict):
                self._var_cmds[std_name] = {}

        for fn in fns:
            with xr.open_dataset(fn) as ds:
                try:
                    self.nc_latitude = ds.cf["latitude"]
                    self.nc_longitude = ds.cf["longitude"]
                    self.nc_elevation = ds.cf["vertical"]
                except KeyError:
                    # Will try to compute values later from first HRU (in self.to_rv)
                    pass

                # Check if any alternate variable name is in the file.
                for std_name in RVT.NC_VARS:
                    for var_name in [std_name] + RVT.NC_VARS[std_name]["alts"]:
                        if var_name not in ds.data_vars:
                            continue
                        nc_var = ds[var_name]
                        self.add_nc_variable(
                            name=std_name,
                            file_name_nc=fn,
                            data_type=RVT.NC_VARS[std_name]["raven"],
                            var_name_nc=var_name,
                            dim_names_nc=nc_var.dims,
                            units=nc_var.attrs.get("units"),
                        )
                        self._number_grid_cells = int(nc_var.size / len(ds["time"]))
                        break

    def update(self, key, value):
        if key in self._var_cmds:
            self._var_cmds[key].update(value)
            return True
        elif key == "nc_index":
            self.nc_index = value
            return True
        return False

    def to_rv(self):
        """
        IMPORTANT NOTE: as this method is called at the last moment in the model lifecycle,
        we can take the occasion to inject in the data structure some values that are guaranteed
        to be there (for instance we can assume that the RVH data is fully specified).
        """

        d = {
            "gauge": "",
            "forcing_list": "",
            "observed_data": "",
        }

        use_gauge = any(type(cmd) is DataCommand for cmd in self._var_cmds.values())
        if use_gauge:
            data = []
            for var, cmd in self._var_cmds.items():
                if cmd and not isinstance(cmd, ObservationDataCommand):
                    data.append(cmd)
            lat = (
                self._nc_latitude[self.nc_index]
                if self._nc_latitude
                else self._config.rvh.hrus[0].latitude
            )
            lon = (
                self._nc_longitude[self.nc_index]
                if self._nc_longitude
                else self._config.rvh.hrus[0].longitude
            )
            elev = (
                self._nc_elevation[self.nc_index]
                if self._nc_elevation
                else self._config.rvh.hrus[0].elevation
            )

            d["gauge"] = GaugeCommand(
                latitude=lat,
                longitude=lon,
                elevation=elev,
                rain_correction=self.rain_correction,
                snow_correction=self.snow_correction,
                monthly_ave_evaporation=self.monthly_ave_evaporation,
                monthly_ave_temperature=self.monthly_ave_temperature,
                data=data,
            )
        else:
            # Construct default grid weights applying equally to all HRUs
            data = [(hru.hru_id, self.nc_index, 1.0) for hru in self._config.rvh.hrus]
            gw = self.grid_weights or GridWeightsCommand(
                number_hrus=len(data),
                number_grid_cells=self._number_grid_cells,
                data=data,
            )
            cmds = []
            for var, cmd in self._var_cmds.items():
                if cmd and not isinstance(cmd, ObservationDataCommand):
                    # TODO: implement a RedirectToFile mechanism to avoid inlining the grid weights
                    # multiple times as we do here
                    if len(cmd.grid_weights.data) == 1:
                        cmd.grid_weights = gw
                    cmds.append(cmd)
            d["forcing_list"] = "\n".join(map(str, cmds))

        # QUESTION: is it possible to have (and if yes should we support) more than 1
        # observation variable? For now we don't.
        for cmd in self._var_cmds.values():
            if isinstance(cmd, ObservationDataCommand):
                # Search for the gauged SB, not sure what should happen when there are
                # more than one (should it be even supported?)
                for sb in self._config.rvh.subbasins:
                    if sb.gauged:
                        cmd.subbasin_id = sb.subbasin_id
                        break
                else:
                    raise Exception(
                        "Could not find an outlet subbasin for observation data"
                    )
                # Set the :StationxIdx (which starts at 1)
                cmd.index = self.nc_index + 1
                d["observed_data"] = cmd
                break

        return dedent(self.tmpl).format(**d)


#########
# O s t #
#########


class OST(RV):

    tmpl = """
    """

    def __init__(self, config):
        super().__init__(config)

        self._max_iterations = None
        self._random_seed = None
        self.lowerBounds = None
        self.upperBounds = None
        self.algorithm = None
        # If there's an OstRandomNumbers.txt file this is its path
        self.random_numbers_path = None

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

    @property
    def identifier(self):
        return self._config.identifier

    def to_rv(self):
        # Get those from RVI (there's probably a better way to do this!)
        self.run_name = self._config.rvi.run_name
        self.run_index = self._config.rvi.run_index

        # Attributes
        a = list(filter(lambda x: not x.startswith("_"), self.__dict__))

        # Properties
        p = list(
            filter(
                lambda x: isinstance(getattr(self.__class__, x, None), property),
                dir(self),
            )
        )

        d = {attr: getattr(self, attr) for attr in a + p}

        return dedent(self.tmpl).format(**d)


class Config:
    def __init__(self, **kwargs):
        self.rvc = RVC(self)
        self.rvh = RVH(self)
        self.rvi = RVI(self)
        self.rvp = RVP(self)
        self.rvt = RVT(self)
        self.ost = OST(self)
        self.identifier = None
        self.update(**kwargs)

    def update(self, key=None, value=None, **kwargs):
        def _update_single(key, value):
            updated = False
            for rv in [self.rvc, self.rvi, self.rvh, self.rvp, self.rvt, self.ost]:
                # Note that in certain cases we might need to update a key
                # for more than one rv object.
                if rv.update(key, value):
                    updated = True
            if not updated:
                raise AttributeError(
                    f"No field named `{key}` found in any RV* conf class"
                )

        identifier = kwargs.pop("identifier", None)
        if identifier:
            self.identifier = identifier

        if key is None and value is None:
            for k, v in kwargs.items():
                _update_single(k, v)
        else:
            _update_single(key, value)

    def set_rv_file(self, fn):
        fn = Path(fn)
        if fn.name == "OstRandomNumbers.txt":
            self.ost.random_numbers_path = fn
        else:
            rvx = fn.suffixes[0][1:]  # get first suffix: eg.g. .rvt[.tpl]
            rvo = getattr(self, rvx, None) or self.ost
            rvo.content = fn.read_text()
            rvo.is_ostrich_tmpl = fn.suffixes[-1] == ".tpl"
            if not self.identifier:
                # Get the "true" stem if there are more than one suffixes
                identifier = fn
                while identifier.suffixes:
                    identifier = Path(identifier.stem)
                self.identifier = identifier.as_posix()
            # This is a sorry hack: I want to have rvi.run_name have a default of "run"
            # because I don't want to burden the user with setting it.. but the problem
            # is that externally supplied rv files might not have it (to be discussed)
            self.rvi.run_name = None
