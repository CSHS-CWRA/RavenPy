import collections
import datetime as dt
from abc import ABC, abstractmethod
from dataclasses import replace
from enum import Enum
from pathlib import Path
from textwrap import dedent
from typing import Any, Dict, List, Optional, Tuple, Union, cast

import cf_xarray
import cftime
import numpy as np
import xarray as xr
from numpy.distutils.misc_util import is_sequence

from ravenpy import __version__
from ravenpy.config.commands import (
    HRU,
    BaseDataCommand,
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
    ReservoirCommand,
    SBGroupPropertyMultiplierCommand,
    SoilClassesCommand,
    SoilProfilesCommand,
    StationForcingCommand,
    Sub,
    SubBasinGroupCommand,
    SubBasinsCommand,
    VegetationClassesCommand,
)


class RV(ABC):

    # This header will be prepended to all RV files when they are rendered
    tmpl_header = """
    ###########################################################################################################
    :FileType          {rv_type} ASCII Raven {raven_version}
    :WrittenBy         PAVICS RavenPy {ravenpy_version} based on setups provided by James Craig and Juliane Mai
    :CreationDate      {date}{model_and_description}
    #----------------------------------------------------------------------------------------------------------
    """

    def __init__(self, config, **kwds):
        # Each RV has a reference to their parent object in order to access sibling RVs.
        self._config = config

        self.is_ostrich_tmpl = False

        # This variable contains the RV file content when it was set from a file; if still
        # None at the moment Raven is called, it means the corresponding RV must be rendered
        # with the `to_rv` method.
        self.content = None

        # This contains extra attributes that might be used with a customized template
        # (currently used with HBVEC and MOHYSE emulators, for values in their RVH)
        self._extra_attributes = {}

    def update(self, key, value):
        if hasattr(self, key):
            setattr(self, key, value)
            return True
        return False

    def set_extra_attributes(self, **kwargs):
        for k, v in kwargs.items():
            self._extra_attributes[k] = v

    def get_extra_attribute(self, k):
        return self._extra_attributes[k]

    def set_tmpl(self, tmpl=None, is_ostrich=False):
        self.tmpl = tmpl or self.tmpl  # type: ignore
        self.is_ostrich_tmpl = is_ostrich

    @property
    @abstractmethod
    def tmpl(self):
        pass

    @abstractmethod
    def to_rv(self, s: str, rv_type: str) -> str:
        if not self._config:
            # In the case where the RV file has been created outside the context of a
            # Config object, don't include the header
            return s
        d = {
            "rv_type": rv_type,
            "raven_version": self._config.model.raven_version,
            "ravenpy_version": __version__,
            "date": dt.datetime.now().isoformat(),
            "model_and_description": "",
        }
        model = ""
        description = self._config.model.description or ""
        if self._config.model.__class__.__name__ not in ["Raven", "Ostrich"]:
            model = f"Emulation of {self._config.model.__class__.__name__}"
        model_and_description = list(filter(None, [model, description]))
        if model_and_description:
            model_and_description = ": ".join(model_and_description)
            d["model_and_description"] = f"\n#\n# {model_and_description}"
        return dedent(self.tmpl_header.lstrip("\n")).format(**d) + s


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

        d.update(self._extra_attributes)

        return super().to_rv(dedent(self.tmpl.lstrip("\n")).format(**d), "RVC")


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
        self.hrus: Tuple[HRU, ...] = ()
        self.subbasins: Tuple[Sub, ...] = ()
        self.land_subbasin_ids: Tuple[int, ...] = ()
        self.land_subbasin_property_multiplier: Optional[
            SBGroupPropertyMultiplierCommand
        ] = None
        self.lake_subbasin_ids: Tuple[int, ...] = ()
        self.lake_subbasin_property_multiplier: Optional[
            SBGroupPropertyMultiplierCommand
        ] = None
        self.reservoirs: Tuple[ReservoirCommand, ...] = ()

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

        d.update(self._extra_attributes)

        return super().to_rv(dedent(self.tmpl.lstrip("\n")).format(**d), "RVH")


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
    {evaluation_periods}
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
        self.run_name: Optional[str] = "run"
        self.run_index = 0
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
        self._evaluation_periods = []
        self._suppress_output = False

    def configure_from_nc_data(self, fns):
        with xr.open_mfdataset(fns, combine="by_coords") as ds:
            start, end = ds.indexes["time"][0], ds.indexes["time"][-1]
            cal = ds.time.encoding.get("calendar", "standard")

        if self.start_date in [None, dt.datetime(1, 1, 1)]:
            self.start_date = start

        if self.end_date in [None, dt.datetime(1, 1, 1)]:
            self.end_date = end

        self.calendar = RVI.CalendarOptions(cal.upper())

    @property
    def raven_version(self):
        return self._config.model.raven_version

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

        if not is_sequence(values):
            values = [values]
        ms = []
        for v in values:
            v = v.upper() if isinstance(v, str) else v.value
            ms.append(RVI.EvaluationMetrics(v))
        self._evaluation_metrics = ms

    @property
    def evaluation_periods(self):
        """:EvaluationPeriod option. Instantiate with list of EvaluationPeriod commands."""
        return "\n".join([str(p) for p in self._evaluation_periods])

    @evaluation_periods.setter
    def evaluation_periods(self, values):
        if not isinstance(values, (list, set, tuple)):
            values = [values]
        self._evaluation_periods = values

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

        d["identifier"] = self._config.model.identifier
        d["now"] = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        d.update(self._extra_attributes)

        t = (
            dedent(self._pre_tmpl.lstrip("\n"))
            + dedent(self.tmpl.lstrip("\n"))
            + dedent(self._post_tmpl.lstrip("\n"))
        )

        return super().to_rv(t.format(**d), "RVI")


##########
# R V P  #
##########


class RVP(RV):

    # This is expected to be defined by the emulators.
    tmpl = """
    """

    def __init__(self, config):
        super().__init__(config)

        # Model specific params
        self.params = None

        self.soil_classes: Tuple[SoilClassesCommand.Record, ...] = ()
        self.soil_profiles: Tuple[SoilProfilesCommand.Record, ...] = ()
        self.vegetation_classes: Tuple[VegetationClassesCommand.Record, ...] = ()
        self.land_use_classes: Tuple[LandUseClassesCommand.Record, ...] = ()
        self.channel_profiles: Tuple[ChannelProfileCommand, ...] = ()
        self.avg_annual_runoff: Optional[float] = None

    def update(self, key, value):
        if key == "params":
            if is_sequence(value):
                self.params = self._config.model.Params(*value)
            else:
                assert isinstance(value, self._config.model.Params)
                self.params = value
            return True
        else:
            return super().update(key, value)

    def to_rv(self):
        d = {
            "params": self.params,
            "soil_classes": SoilClassesCommand(self.soil_classes),
            "soil_profiles": SoilProfilesCommand(self.soil_profiles),
            "vegetation_classes": VegetationClassesCommand(self.vegetation_classes),
            "land_use_classes": LandUseClassesCommand(self.land_use_classes),
            "channel_profiles": "\n\n".join(map(str, self.channel_profiles)),
            "avg_annual_runoff": f":AvgAnnualRunoff {self.avg_annual_runoff}"
            if self.avg_annual_runoff
            else "",
        }

        d.update(self._extra_attributes)

        return super().to_rv(dedent(self.tmpl.lstrip("\n")).format(**d), "RVP")


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

        # These are customized variable attributes specified by the user
        self._var_specs: Dict[str, Dict[str, Any]] = {k: {} for k in RVT.NC_VARS.keys()}

        # These are the actual variable as `commands.BaseDataCommand` objects
        self._var_cmds: Dict[str, Optional[BaseDataCommand]] = {
            k: None for k in RVT.NC_VARS.keys()
        }

        # Specifies whether the variables must be configured using file NC data
        self._auto_nc_configure = True

        self.nc_index = 0
        self.grid_weights = None
        self.rain_correction = None
        self.snow_correction = None
        self.monthly_ave_evaporation = None
        self.monthly_ave_temperature = None

        self._nc_latitude: Optional[xr.DataArray] = None
        self._nc_longitude: Optional[xr.DataArray] = None
        self._nc_elevation: Optional[xr.DataArray] = None
        self._number_grid_cells = 0

    def _add_nc_variable(self, **kwargs):
        std_name = kwargs.get("name", kwargs["var_name_nc"])
        # If the name is not the standard one, search for it
        # TODO: reorganize NC_VARS so that the keys are the Raven names
        if std_name not in RVT.NC_VARS:
            for sn, rec in RVT.NC_VARS.items():
                if rec["raven"] == kwargs["data_type"] or std_name in rec["alts"]:
                    std_name = sn
                    break
            else:
                assert False, f"{std_name} not found in the list of standard names"
        is_obs_var = kwargs.pop("is_observation", False)
        cmd: BaseDataCommand
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

        spec = self._var_specs[std_name]
        self._var_cmds[std_name] = replace(cmd, **spec)

    def set_nc_variables(self, nc_variables):
        """
        This is meant for manually setting the variables, and should prevent
        automatic configuration from an nc file.
        """
        for nc_var in nc_variables:
            self._add_nc_variable(**nc_var)
        self._auto_nc_configure = False

    def configure_from_nc_data(self, fns):

        assert self._auto_nc_configure is True

        self._var_cmds = {k: None for k in RVT.NC_VARS.keys()}

        for fn in fns:
            with xr.open_dataset(fn) as ds:
                try:
                    self._nc_latitude = ds.cf["latitude"]
                    self._nc_longitude = ds.cf["longitude"]
                    self._nc_elevation = ds.cf["vertical"]
                except KeyError:
                    # Will try to compute values later from first HRU (in self.to_rv)
                    pass

                # Check if any alternate variable name is in the file.
                for std_name in RVT.NC_VARS:
                    for var_name in [std_name] + RVT.NC_VARS[std_name]["alts"]:  # type: ignore
                        if var_name not in ds.data_vars:
                            continue
                        nc_var = ds[var_name]
                        self._add_nc_variable(
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
        if key in self._var_specs:
            self._var_specs[key].update(value)
            return True
        return super().update(key, value)

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
            data_cmds = []
            for var, cmd in self._var_cmds.items():
                if cmd and not isinstance(cmd, ObservationDataCommand):
                    cmd = cast(DataCommand, cmd)
                    data_cmds.append(cmd)

            if (
                self._nc_latitude
                and self._nc_latitude.shape
                and len(self._nc_latitude) > self.nc_index
            ):
                lat = self._nc_latitude.values[self.nc_index]
            else:
                lat = self._config.rvh.hrus[0].latitude

            if (
                self._nc_longitude
                and self._nc_longitude.shape
                and len(self._nc_longitude) > self.nc_index
            ):
                lon = self._nc_longitude.values[self.nc_index]
            else:
                lon = self._config.rvh.hrus[0].longitude

            if (
                self._nc_elevation
                and self._nc_elevation
                and len(self._nc_elevation) > self.nc_index
            ):
                elev = self._nc_elevation.values[self.nc_index]
            else:
                elev = self._config.rvh.hrus[0].elevation

            d["gauge"] = GaugeCommand(
                latitude=lat,
                longitude=lon,
                elevation=elev,
                rain_correction=self.rain_correction,
                snow_correction=self.snow_correction,
                monthly_ave_evaporation=self.monthly_ave_evaporation,
                monthly_ave_temperature=self.monthly_ave_temperature,
                data_cmds=tuple(data_cmds),
            )  # type: ignore
        else:
            # Construct default grid weights applying equally to all HRUs
            data = [(hru.hru_id, self.nc_index, 1.0) for hru in self._config.rvh.hrus]
            gw = self.grid_weights or GridWeightsCommand(
                number_hrus=len(data),
                number_grid_cells=self._number_grid_cells,
                data=tuple(data),
            )
            cmds = []
            for var, cmd in self._var_cmds.items():
                if cmd and not isinstance(cmd, ObservationDataCommand):
                    # TODO: implement a RedirectToFile mechanism to avoid inlining the grid weights
                    # multiple times as we do here
                    cmd = cast(Union[GriddedForcingCommand, StationForcingCommand], cmd)
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
                d["observed_data"] = cmd  # type: ignore
                break

        return super().to_rv(dedent(self.tmpl.lstrip("\n")).format(**d), "RVT")


#########
# O s t #
#########


class OST(RV):

    tmpl = """
    """

    # Multiplier applied to metric before passing to minimization algorithm.
    _evaluation_metrics_multiplier = dict(
        NASH_SUTCLIFFE=-1,
        LOG_NASH=-1,
        RMSE=1,
        PCT_BIAS="Not Supported",
        ABSERR=1,
        ABSMAX=1,
        PDIFF=1,
        TMVOL=1,
        RCOEFF=1,
        NSC=-1,
        KLING_GUPTA=-1,
    )

    def __init__(self, config):
        super().__init__(config)

        self._max_iterations = None
        self._random_seed = None
        self.lowerBounds = None
        self.upperBounds = None
        self.algorithm = None
        # If there's an OstRandomNumbers.txt file this is its path
        self.random_numbers_path = None

    def update(self, key, value):
        if key in ["lowerBounds", "upperBounds"]:
            if is_sequence(value):
                setattr(self, key, self._config.model.Params(*value))
            else:
                assert isinstance(value, self._config.model.Params)
                setattr(self, key, value)
            return True
        else:
            return super().update(key, value)

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
    def evaluation_metric_multiplier(self):
        """For Ostrich."""
        m = self._config.rvi._evaluation_metrics[0]
        if m.name == "PCT_BIAS":
            raise ValueError(
                "PCT_BIAS cannot be properly minimized. Select another evaluation metric."
            )
        return self._evaluation_metrics_multiplier[m.value]

    @property
    def identifier(self):
        return self._config.model.identifier

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

        d.update(self._extra_attributes)

        return super().to_rv(dedent(self.tmpl.lstrip("\n")).format(**d), "OST")


class Config:
    def __init__(self, model, **kwargs):
        self.model = model
        self.rvc = RVC(self)
        self.rvh = RVH(self)
        self.rvi = RVI(self)
        self.rvp = RVP(self)
        self.rvt = RVT(self)
        self.ost = OST(self)
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

    def set_rv_file(self, fn):
        fn = Path(fn)
        if fn.name == "OstRandomNumbers.txt":
            self.ost.random_numbers_path = fn
        else:
            rvx = fn.suffixes[0][1:]  # get first suffix: eg.g. .rvt[.tpl]
            rvo = getattr(self, rvx, None) or self.ost
            rvo.content = fn.read_text()
            rvo.is_ostrich_tmpl = fn.suffixes[-1] == ".tpl"
            # This is a sorry hack: I want to have rvi.run_name have a default of "run"
            # because I don't want to burden the user with setting it.. but the problem
            # is that externally supplied rv files might not have it (to be discussed)
            self.rvi.run_name = None
