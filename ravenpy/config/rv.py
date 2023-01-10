from dataclasses import asdict, fields
from enum import Enum
from typing import Any

from pydantic.dataclasses import dataclass

from . import options
from .base import RavenCommand, parse_symbolic
from .commands import *

"""
# Notes

Either all attributes have defaults, or none have. Otherwise, it's bound to cause issues when creating emulators from multiple subclasses.


"""


class RV:
    """"""

    _header = """
    ###########################################################################################################
    :FileType          {rv_type} ASCII Raven {raven_version}
    :WrittenBy         PAVICS RavenPy {ravenpy_version} based on setups provided by James Craig and Juliane Mai
    :CreationDate      {date}{model_and_description}
    #----------------------------------------------------------------------------------------------------------
    """

    def to_rv(self):
        out = ""
        for key, val in self.__dict__.items():
            if val is not None and not key.startswith("_"):
                out += val.to_rv()
        return out


@dataclass
class RVC(RV):
    hru_states: HRUStateVariableTableCommand = HRUStateVariableTableCommand()
    basin_states: BasinStateVariablesCommand = BasinStateVariablesCommand()


@dataclass
class RVH(RV):
    subbasins: SubBasinsCommand = None
    hrus: HRUsCommand = None
    land_subbasin_group: SubBasinGroupCommand = None
    land_subbasin_property_multiplier: SBGroupPropertyMultiplierCommand = None
    lake_subbasin_group: SubBasinGroupCommand = None
    lake_subbasin_property_multiplier: SBGroupPropertyMultiplierCommand = None
    reservoirs: Tuple[ReservoirCommand] = None


@dataclass
class RVI(RV):
    calendar: Calendar = Calendar("STANDARD")
    run_name: RunName = RunName("run")
    start_date: StartDate = None
    end_date: EndDate = None
    time_step: TimeStep = TimeStep(1.0)
    evaluation_metrics: EvaluationMetrics = EvaluationMetrics(
        ("NASH_SUTCLIFFE", "RMSE")
    )
    evaluation_periods: Tuple[EvaluationPeriod] = None
    write_netcdf_format: WriteNetcdfFormat = WriteNetcdfFormat(True)
    silent_mode: SilentMode = SilentMode(True)
    pavics_mode: PavicsMode = PavicsMode(True)
    suppress_output: SuppressOutput = SuppressOutput(False)
    write_forcing_functions: WriteForcingFunctions = WriteForcingFunctions(False)
    custom_output: Tuple[CustomOutput] = None
    netcdf_attribute: NetCDFAttribute = None
    rain_snow_fraction: RainSnowFraction = RainSnowFraction("RAINSNOW_DATA")
    routing: Routing = Routing("ROUTE_NONE")
    evaporation: Evaporation = None
    ow_evaporation: OW_Evaporation = None


@dataclass
class RVP(RV):
    params: Any = None
    soil_classes: SoilClasses = None
    soil_profiles: SoilProfiles = SoilProfiles()
    vegetation_classes: VegetationClasses = VegetationClasses()
    land_use_classes: LandUseClassesCommand = LandUseClassesCommand()
    channel_profiles: ChannelProfileCommand = ChannelProfileCommand()
    soil_parameter_list: SoilParameterListCommand = SoilParameterListCommand()
    avg_annual_runoff: AvgAnnualRunoff = None
    rain_snow_transition: RainSnowTransition = None
    air_snow_coeff: AirSnowCoeff = None
    avg_annual_snow: AvgAnnualSnow = None
    precipitation_lapse_rate: PrecipitationLapseRate = None
    adiabatic_lapse_rate: AdiabaticLapseRate = None


@dataclass
class RVT(RV):
    gauge: Tuple[GaugeCommand] = None
    forcing: Tuple[Union[GriddedForcingCommand, StationForcingCommand]] = None
    observed_data: Tuple[ObservationDataCommand] = None


@dataclass
class Config(RVI, RVC, RVH, RVT, RVP):
    @validator("*", pre=True)
    def assign_symbolic(cls, v, values, config, field):
        print(field.name)
        if field.name != "params":
            return parse_symbolic(v, **asdict(values["params"]))
        return v

    def to_rv(self, rv: str):
        """Return RV configuration text."""
        rvs = {b.__name__: b for b in Config.__bases__}
        cls = rvs[rv.upper()]

        p = {f.name: self.__dict__[f.name] for f in fields(cls)}
        rv = cls(**p)
        return rv.to_rv()


@dataclass
class TestConfig(Config):
    calendar: Calendar = Calendar("JULIAN")
