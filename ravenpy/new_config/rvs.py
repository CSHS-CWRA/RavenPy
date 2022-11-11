import datetime as dt
from enum import Enum
from typing import Dict, Sequence, TypeVar, Union

import cftime
from pydantic import BaseModel, Extra, Field, validator
from pydantic.dataclasses import dataclass

from ravenpy.config import options as o

from . import commands as rc
from .base import RV


class RVI(RV):
    # Run parameters
    run_name: str = Field("run", alias="RunName")
    calendar: o.Calendar = Field("STANDARD", alias="Calendar")
    start_date: cftime.datetime = Field(None, alias="StartDate")
    end_date: cftime.datetime = Field(None, alias="EndDate")
    duration: float = Field(1, alias="Duration")
    time_step: float = Field(1.0, alias="TimeStep")
    evaluation_period: Sequence[rc.EvaluationPeriod] = Field(
        None, alias="EvaluationPeriod"
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
    # oro_temp_correct: rc.OroTempCorrect = None
    # oro_precip_correct: rc.OroPrecipCorrect = None
    # oro_pet_correct: rc.OroPETCorrect = None
    cloud_cover_method: o.CloudCoverMethod = Field(None, alias="CloudCoverMethod")
    precip_icept_frac: o.PrecipIceptFract = Field(None, alias="PrecipIceptFract")
    subdaily_method: o.SubdailyMethod = Field(None, alias="SubdailyMethod")
    mim: o.MonthlyInterpolationMethod = Field(None, alias="MonthlyInterpolationMethod")
    soil_model: Union[int, rc.SoilModel] = Field(None, alias="SoilModel")
    lake_storage: o.StateVariables = Field(None, alias="LakeStorage")

    # Options
    alias: Dict[str, str] = Field(None, alias="Alias")
    netcdf_attribute: Dict[str, str] = Field(None, alias="NetCDFAttribute")
    # hydrologic_processes: rc.HydrologicProcesses

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
        return cftime._cftime.DATE_TYPES[values["calendar"].lower()](
            *val.timetuple()[:6]
        )

    class Config:
        arbitrary_types_allowed = True


class RVT(RV):
    gauge: Sequence[rc.Gauge] = ()
    data: Sequence[rc.Data] = ()
    observation: Sequence[rc.ObservationData] = ()


class RVP(RV):
    # params = rc.Params
    soil_classes: rc.SoilClasses = Field(None, alias="SoilClasses")
    soil_profiles: rc.SoilProfiles = Field(None, alias="SoilProfiles")
    vegetation_classes: rc.VegetationClassesCommand = Field(
        None, alias="VegetationClasses"
    )
    land_use_classes: rc.LandUseClassesCommand = Field(None, alias="LandUseClasses")
    soil_parameter_list: rc.SoilParameterList = Field(None, alias="SoilParameterList")
    land_use_parameter_list: rc.LandUseParameterList = Field(
        None, alias="LandUseParameterList"
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        None, alias="VegetationParameterList"
    )

    channel_profiles: Sequence[rc.ChannelProfile] = ()

    # TODO: create list of all available parameters to constrain key
    global_parameter: Dict[str, str] = Field(None, alias="GlobalParameter")
    rain_snow_transition: rc.RainSnowTransition = Field(
        None, alias="RainSnowTransition"
    )


class RVC(RV):
    hru_states: rc.HRUStateVariableTable = ()
    basin_states: rc.BasinStateVariables = ()


class RVH(RV):
    subbasins: Sequence[rc.SubBasin] = (rc.SubBasin(),)
    hrus: Sequence[rc.HRU] = ()
    land_subbasin_group: Sequence[rc.SubBasinGroup] = ()
    lake_subbasin_group: Sequence[rc.SubBasinGroup] = ()
    land_subbasin_property_multiplier: rc.SBGroupPropertyMultiplierCommand = None
    lake_subbasin_property_multiplier: rc.SBGroupPropertyMultiplierCommand = None
    reservoirs: Sequence[rc.Reservoir] = ()


"""
@dataclass
class Emulator:
    rvi: RVI = RVI()
    rvp: RVP = RVP()
    rvc: RVC = RVC()
    rvh: RVH = RVH()
    rvt: RVT = RVT()

    def write(self, stem="raven", path="."):
        p = Path(path)
        for k, v in asdict(self).items():
            fn = p / f"{stem}.{k}"
            fn.write_text(v.to_rv())
"""
