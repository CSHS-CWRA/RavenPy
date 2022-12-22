from dataclasses import asdict, fields
from enum import Enum

from pydantic.dataclasses import dataclass

from . import options
from .base import RavenCommand
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
    params: Enum = None
    soil_classes: SoilClasses = SoilClassesCommand()
    soil_profiles: SoilProfilesCommand = SoilProfilesCommand()
    vegetation_classes: VegetationClassesCommand = VegetationClassesCommand()
    land_use_classes: LandUseClassesCommand = LandUseClassesCommand()
    channel_profiles: ChannelProfileCommand = ChannelProfileCommand()
    soil_parameter_list: SoilParameterListCommand = SoilParameterListCommand()
    avg_annual_runoff: AvgAnnualRunoff = None
    rain_snow_transition: RainSnowTransition = None
    air_snow_coeff: AirSnowCoeff = None
    avg_annual_snow: AvgAnnualSnow = None
    precipitation_lapse_rate: PrecipitationLapseRate = None
    adiabatic_lapse_rate: AdiabaticLapseRate


@dataclass
class RVT(RV):
    gauge: Tuple[GaugeCommand] = None
    forcing: Tuple[Union[GriddedForcingCommand, StationForcingCommand]] = None
    observed_data: Tuple[ObservationDataCommand] = None


class Config(RVP, RVI, RVC, RVH, RVT):
    def to_rv(self, rv: str):
        """Return RV configuration text."""
        rvs = {b.__name__: b for b in self.__class__.__bases__}
        cls = rvs[rv.upper()]

        p = {f.name: self.__dict__[f.name] for f in fields(cls)}
        rv = cls(**p)
        return rv.to_rv()


class GR4JCN:
    air_snow_coeff: AirSnowCoeff = None  # 1 - CEMANEIGE_X2
    avg_annual_snow: AvgAnnualSnow = None  # CEMANEIGE_X1
    rain_snow_transition: RainSnowTransition = RainSnowTransition(0, 1.0)
    precipitation_lapse_rate: PrecipitationLapseRate = PrecipitationLapseRate(0.0004)
    adiabatic_lapse_rate: AdiabaticLapseRate = AdiabaticLapseRate(0.0065)
    soil_classes: SoilClasses(
        names=("SOIL_PROD", "SOIL_ROUT", "SOIL_TEMP", "SOIL_GW", "AQUIFER")
    )
    soil_parameter_list: SoilParameterListCommand = SoilParameterListCommand(
        names=("[DEFAULT]",),
        records=[PL(name="[DEFAULT]", vals=(1, p.GR4J_X3, p.GR4J_X2))],
    )
