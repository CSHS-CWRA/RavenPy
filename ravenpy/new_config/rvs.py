import datetime as dt
from dataclasses import asdict
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Sequence, Type, TypeVar, Union

import cftime
from pydantic import BaseModel, Extra, Field, validator
from pydantic.dataclasses import dataclass

from . import commands as rc
from . import options as o
from .base import RV, parse_symbolic


class RVI(RV):
    # Run parameters
    run_name: str = Field("run", alias="RunName")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    start_date: cftime.datetime = Field(None, alias="StartDate")
    end_date: cftime.datetime = Field(None, alias="EndDate")
    duration: float = Field(None, alias="Duration")
    time_step: float = Field(1.0, alias="TimeStep")
    evaluation_period: Sequence[rc.EvaluationPeriod] = Field(
        None, alias="EvaluationPeriod"
    )
    evaluation_metrics: Sequence[o.EvaluationMetrics] = Field(
        None, alias="EvaluationMetrics"
    )

    # Model description
    routing: o.Routing = Field(None, alias="Routing")
    # method: rc.Method = None
    # interpolation: rc.Interpolation = None
    catchment_route: o.CatchmentRoute = Field(None, alias="CatchmentRoute")
    evaporation: o.Evaporation = Field(None, alias="Evaporation")
    ow_evaporation: o.Evaporation = Field(None, alias="OW_Evaporation")
    sw_radiation_method: o.SWRadiationMethod = Field(None, alias="SWRadiationMethod")
    sw_cloud_correct: o.SWCloudCorrect = Field(None, alias="SWCloudCorrect")
    sw_canopy_correct: o.SWCanopyCorrect = Field(None, alias="SWCanopyCorrect")
    # lw_radiation_method: rc.LWRadiationMethod = None
    windspeed_method: o.WindspeedMethod = Field(None, alias="WindspeedMethod")
    rain_snow_fraction: o.RainSnowFraction = Field(None, alias="RainSnowFraction")
    potential_melt_method: o.PotentialMeltMethod = Field(
        None, alias="PotentialMeltMethod"
    )
    oro_temp_correct: o.OroTempCorrect = Field(None, alias="OroTempCorrect")
    oro_precip_correct: o.OroPrecipCorrect = Field(None, alias="OroPrecipCorrect")
    oro_pet_correct: o.OroPETCorrect = Field(None, alias="OroPETCorrect")
    cloud_cover_method: o.CloudCoverMethod = Field(None, alias="CloudCoverMethod")
    precip_icept_frac: o.PrecipIceptFract = Field(None, alias="PrecipIceptFract")
    subdaily_method: o.SubdailyMethod = Field(None, alias="SubdailyMethod")
    mim: o.MonthlyInterpolationMethod = Field(None, alias="MonthlyInterpolationMethod")
    soil_model: Union[int, rc.SoilModel] = Field(None, alias="SoilModel")
    lake_storage: o.StateVariables = Field(None, alias="LakeStorage")
    # alias: Dict[str, str] = Field(None, alias="Alias")
    hydrologic_processes: Sequence[rc.Process] = Field(
        None, alias="HydrologicProcesses"
    )

    # Options
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    netcdf_attribute: Dict[str, str] = Field(None, alias="NetCDFAttribute")

    custom_output: rc.CustomOutput = Field(None, alias="CustomOutput")
    direct_evaporation: bool = Field(
        None,
        alias="DirectEvaporation",
        description="Rainfall is automatically reduced through evapotranspiration up to the limit of the calculated PET.",
    )
    deltares_fews_mode: bool = Field(None, alias="DeltaresFEWSMode")
    debug_mode: bool = Field(None, alias="DebugMode")
    dont_write_watershed_storage: bool = Field(
        None,
        alias="DontWriteWatershedStorage",
        description="Do not write watershed storage variables to disk.",
    )
    pavics_mode: bool = Field(None, alias="PavicsMode")
    silent_mode: bool = Field(None, alias="SilentMode")
    suppress_output: bool = Field(
        None,
        alias="SuppressOutput",
        description="Write minimal output to disk when enabled.",
    )
    write_forcing_functions: bool = Field(
        None,
        alias="WriteForcingFunctions",
        description="Write watershed averaged forcing functions (e.g. rainfall, radiation, PET, etc).",
    )
    write_subbasin_file: bool = Field(None, alias="WriteSubbasinFile")  # Undocumented

    @validator("soil_model", pre=True)
    def init_soil_model(cls, v):
        if isinstance(v, int):
            return rc.SoilModel(n=v)
        return v

    @validator("start_date", "end_date", pre=True)
    def dates2cf(cls, val, values):
        """Convert dates to cftime dates."""
        if val is not None:
            return cftime._cftime.DATE_TYPES[values["calendar"].value.lower()](
                *val.timetuple()[:6]
            )

    class Config:
        arbitrary_types_allowed = True


class RVT(RV):
    gauge: Sequence[rc.Gauge] = Field(None, alias="Gauge")
    # data: Sequence[rc.Data] = ()
    station_forcing: Sequence[rc.StationForcing] = Field(None, alias="StationForcing")
    gridded_forcing: Sequence[rc.GriddedForcing] = Field(None, alias="GriddedForcing")
    observation_data: Sequence[rc.ObservationData] = Field(
        None, alias="ObservationData"
    )


class RVP(RV):
    params: Any
    soil_classes: rc.SoilClasses = Field(None, alias="SoilClasses")
    soil_profiles: Sequence[rc.SoilProfile] = Field(None, alias="SoilProfiles")
    vegetation_classes: rc.VegetationClasses = Field(None, alias="VegetationClasses")
    land_use_classes: rc.LandUseClasses = Field(None, alias="LandUseClasses")
    soil_parameter_list: rc.SoilParameterList = Field(None, alias="SoilParameterList")
    land_use_parameter_list: rc.LandUseParameterList = Field(
        None, alias="LandUseParameterList"
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        None, alias="VegetationParameterList"
    )

    channel_profiles: Sequence[rc.ChannelProfile] = ()

    # TODO: create list of all available parameters to constrain key
    global_parameter: Dict[str, str] = Field({}, alias="GlobalParameter")
    rain_snow_transition: rc.RainSnowTransition = Field(
        None, alias="RainSnowTransition"
    )


class RVC(RV):
    hru_states: rc.HRUStateVariableTable = ()
    basin_states: rc.BasinStateVariables = ()
    uniform_initial_conditions: Dict[str, float] = Field(
        None, alias="UniformInitialConditions"
    )

    @classmethod
    def from_solution(cls, solution: str):
        hru_states = rc.HRUStateVariableTable.parse(solution)
        basin_states = rc.BasinStateVariables.parse(solution)
        return cls(hru_states=hru_states, basin_states=basin_states)


class RVH(RV):
    subbasins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    hrus: rc.HRUs = Field(None, alias="HRUs")
    land_subbasin_group: Sequence[rc.SubBasinGroup] = ()
    lake_subbasin_group: Sequence[rc.SubBasinGroup] = ()
    land_subbasin_property_multiplier: rc.SBGroupPropertyMultiplierCommand = None
    lake_subbasin_property_multiplier: rc.SBGroupPropertyMultiplierCommand = None
    reservoirs: Sequence[rc.Reservoir] = ()


class Config(RVI, RVC, RVH, RVT, RVP):
    @validator("*", pre=True)
    def assign_symbolic(cls, v, values, config, field):
        if field.name != "params":
            return parse_symbolic(v, **asdict(values["params"]))
        return v

    @validator("global_parameter", pre=True)
    def update_defaults(cls, v, values, config, field):
        """Some configuration parameters should be updated with user given arguments, not overwritten."""
        return {**cls.__fields__[field.name].default, **v}

    def to_rv(self, rv: str):
        """Return RV configuration text."""
        # Get RV class
        rvs = {b.__name__: b for b in Config.__bases__}
        cls = rvs[rv.upper()]

        # Instantiate RV class
        attrs = dict(self)
        p = {f: attrs[f] for f in cls.__fields__}
        rv = cls(**p)
        return rv.to_rv()

    def write(self, workdir: Union[str, Path], overwrite=False):
        """Write configuration files to disk.

        Parameters
        ----------
        workdir: str, Path
          A directory where rv files will be written to disk.
        overwrite: bool
          If True, overwrite existing configuration files.
        """
        workdir = Path(workdir)
        if not workdir.exists():
            workdir.mkdir()

        for rv in ["rvi", "rvp", "rvc", "rvh", "rvt"]:
            fn = workdir / f"{self.run_name}.{rv}"
            if fn.exists() and not overwrite:
                raise OSError(f"{fn} already exists and would be overwritten.")
            fn.write_text(self.to_rv(rv))
