import collections
import datetime as dt
from abc import ABC, abstractmethod
from dataclasses import is_dataclass, replace
from pathlib import Path
from textwrap import dedent

import cf_xarray
import cftime
import numpy as np
import xarray as xr

from ravenpy.config.commands import (
    AvgAnnualRunoffCommand,
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
    MonthlyAverageCommand,
    ObservationDataCommand,
    RainCorrectionCommand,
    RoutingCommand,
    SnowCorrectionCommand,
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
        # TODO: find something better than this!
        self.is_ostrich_tmpl = False

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
        self.hru_states = {}
        self.basin_states = {}

    def reset(self, **kwargs):
        self.hru_states = {}
        self.basin_states = {}

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

    @staticmethod
    def parse_solution(solution_str):
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
                        out[old_key] = _parser(
                            itertools.chain([line], lines), new_indent
                        )
                    elif new_indent < indent:
                        return out
                    else:
                        if key == "BasinIndex":
                            i, name = value.split(",")
                            i = int(i)
                            out[key][i] = dict(
                                index=i,
                                name=name,
                                **_parser(lines, new_indent + "  ", float),
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
                                list(map(fmt, value.split(",")))
                                if "," in value
                                else fmt(value)
                            )

                    old_key = key
                else:
                    data = line.split(",")
                    i = int(data.pop(0))
                    out["data"][i] = [i] + list(map(float, data))

            return out

        lines = iter(solution_str.splitlines())
        return _parser(lines)

    @staticmethod
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

    def set_from_solution(self, solution_str):

        solution_objs = RVC.parse_solution(solution_str)

        self.hru_states = {}
        self.basin_states = {}

        for index, params in solution_objs["HRUStateVariableTable"]["data"].items():
            self.hru_states[index] = HRUStateVariableTableCommand.Record(*params)

        for index, raw in solution_objs["BasinStateVariables"]["BasinIndex"].items():
            params = {k.lower(): v for (k, v) in raw.items()}
            self.basin_states[index] = BasinIndexCommand(**params)

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

    tmpl = """
    """

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

    def __init__(self, config):
        super().__init__(config)
        self.name = None
        self.area = None
        self.elevation = None
        self.latitude = None
        self.longitude = None
        self.run_name = "run"
        self.run_index = 0
        self.raven_version = "3.0.1 rev#275"
        self.time_step = 1.0

        self._routing = "ROUTE_NONE"
        self._start_date = None
        self._end_date = None
        self._now = None
        self._rain_snow_fraction = "RAINSNOW_DATA"
        self._evaporation = None
        self._ow_evaporation = None
        self._duration = 1
        self._evaluation_metrics = "NASH_SUTCLIFFE RMSE"
        self._suppress_output = False
        self._calendar = "standard"

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
        return self._evaluation_metrics

    @evaluation_metrics.setter
    def evaluation_metrics(self, x):
        if not isinstance(x, str):
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
    def routing(self):
        return RoutingCommand(value=self._routing)

    @routing.setter
    def routing(self, value):
        self._routing = value

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

        if v in RVI.rain_snow_fraction_options:
            self._rain_snow_fraction = v
        else:
            raise ValueError(
                f"Value should be one of {RVI.rain_snow_fraction_options}."
            )

    @property
    def evaporation(self):
        """Evaporation scheme"""
        return self._evaporation

    @evaporation.setter
    def evaporation(self, value):
        v = value.upper()
        if v in RVI.evaporation_options:
            self._evaporation = v
        else:
            raise ValueError(f"Value {v} should be one of {RVI.evaporation_options}.")

    @property
    def ow_evaporation(self):
        """Open-water evaporation scheme"""
        return self._ow_evaporation

    @ow_evaporation.setter
    def ow_evaporation(self, value):
        v = value.upper()
        if v in RVI.evaporation_options:
            self._ow_evaporation = v
        else:
            raise ValueError(f"Value {v} should be one of {RVI.evaporation_options}.")

    @property
    def calendar(self):
        """Calendar"""
        return self._calendar.upper()

    @calendar.setter
    def calendar(self, value):
        if value.upper() in RVI.calendar_options:
            self._calendar = value
        else:
            raise ValueError(f"Value should be one of {RVI.calendar_options}.")

    def _dt2cf(self, date):
        """Convert datetime to cftime datetime."""
        return cftime._cftime.DATE_TYPES[self._calendar.lower()](*date.timetuple()[:6])

    def to_rv(self):

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


##########
# R V P  #
##########


class RVP(RV):

    tmpl = """
    """

    def __init__(self, config):
        super().__init__(config)

        # Model specific params and derived params
        self.params = None
        self.derived_params = None

        self.soil_classes: Tuple[SoilClassesCommand.Record] = ()
        self.soil_profiles: Tuple[SoilProfilesCommand.Record] = ()
        self.vegetation_classes: Tuple[VegetationClassesCommand.Record] = ()
        self.land_use_classes: Tuple[LandUseClassesCommand.Record] = ()
        self.channel_profiles: Tuple[ChannelProfileCommand] = ()
        self.avg_annual_runoff: float = None

    def to_rv(self):
        d = {
            "params": self.params,
            "derived_params": self.derived_params,
            "soil_classes": SoilClassesCommand(self.soil_classes),
            "soil_profiles": SoilProfilesCommand(self.soil_profiles),
            "vegetation_classes": VegetationClassesCommand(self.vegetation_classes),
            "land_use_classes": LandUseClassesCommand(self.land_use_classes),
            "channel_profiles": "\n\n".join(map(str, self.channel_profiles)),
            "avg_annual_runoff": AvgAnnualRunoffCommand(self.avg_annual_runoff),
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

    NC_VARS = {
        "tasmin": {"name": "TEMP_MIN", "alts": ["tasmin", "tmin"]},
        "tasmax": {"name": "TEMP_MAX", "alts": ["tasmax", "tmax"]},
        "tas": {"name": "TEMP_AVE", "alts": ["tas", "t2m"]},
        "rainfall": {"name": "RAINFALL", "alts": ["rainfall", "rain"]},
        "pr": {
            "name": "PRECIP",
            "alts": ["pr", "precip", "prec", "precipitation", "tp"],
        },
        "prsn": {
            "name": "SNOWFALL",
            "alts": ["prsn", "snow", "snowfall", "solid_precip"],
        },
        "evspsbl": {"name": "PET", "alts": ["pet", "evap", "evapotranspiration"]},
        "water_volume_transport_in_river_channel": {
            "name": "HYDROGRAPH",
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
        var_name = kwargs.get("name", kwargs["var_name_nc"])
        is_obs_var = kwargs.pop("is_observation", False)
        if len(kwargs["dim_names_nc"]) == 1:
            if var_name == "water_volume_transport_in_river_channel" or is_obs_var:
                cmd = ObservationDataCommand(**kwargs)
            else:
                cmd = DataCommand(**kwargs)
        elif len(kwargs["dim_names_nc"]) == 2:
            if var_name == "water_volume_transport_in_river_channel" or is_obs_var:
                cmd = ObservationDataCommand(**kwargs)
            else:
                cmd = StationForcingCommand(**kwargs)
        else:
            cmd = GriddedForcingCommand(**kwargs)

        if isinstance(self._var_cmds.get(var_name, None), dict):
            self._var_cmds[var_name] = replace(cmd, **self._var_cmds[var_name])
        else:
            self._var_cmds[var_name] = cmd

    def configure_from_nc_data(self, fns):

        self._var_cmds = {k: {} for k in RVT.NC_VARS.keys()}

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
                for var_name in RVT.NC_VARS:
                    for alt_name in RVT.NC_VARS[var_name]["alts"]:
                        if alt_name not in ds.data_vars:
                            continue
                        nc_var = ds[alt_name]
                        self.add_nc_variable(
                            name=var_name,
                            file_name_nc=fn,
                            data_type=RVT.NC_VARS[var_name]["name"],
                            var_name_nc=alt_name,
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
                d["observed_data"] = cmd
                break

        return dedent(self.tmpl).format(**d)


#########
# O s t #
#########


class Ost(RV):

    tmpl = """
    """

    def __init__(self, config, identifier=None):
        super().__init__(config)

        self._max_iterations = None
        self._random_seed = None
        self.lowerBounds = None
        self.upperBounds = None
        self.algorithm = None

        # TODO: find something better than this
        self.identifier = identifier

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
    def __init__(self, identifier, **kwargs):
        self.rvc = RVC(self)
        self.rvh = RVH(self)
        self.rvi = RVI(self)
        self.rvp = RVP(self)
        self.rvt = RVT(self)
        self.ost = Ost(self, identifier=identifier)
        self.identifier = identifier
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

        if key is None and value is None:
            for k, v in kwargs.items():
                _update_single(k, v)
        else:
            _update_single(key, value)
