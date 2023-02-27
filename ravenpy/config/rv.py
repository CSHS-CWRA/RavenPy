from dataclasses import asdict, fields
from enum import Enum
from typing import Any

from pydantic.dataclasses import dataclass
import ravenpy
from . import options
from .base import RavenCommand, parse_symbolic
from .commands import *
import datetime as dt

"""
# Notes

Either all attributes have defaults, or none have. Otherwise, it's bound to cause issues when creating emulators from multiple subclasses.


- Ditch dataclasses and subclass everything from BaseModel
- Add extra='forbid' and allow_mutation=False.

"""


class RV:
    """"""

    _header = """
    ###########################################################################################################
    :FileType          {rv_type} ASCII Raven
    :WrittenBy         PAVICS RavenPy {ravenpy_version} based on setups provided by James Craig and Juliane Mai
    :CreationDate      {date}
    #----------------------------------------------------------------------------------------------------------

    """

    def to_rv(self):
        out = dedent(self._header.format(rv_type=self.__class__.__name__,
                                  ravenpy_version=ravenpy.__version__,
                                  date=dt.datetime.now().date()))

        for key, val in self.__dict__.items():
            if val is not None and not key.startswith("_"):
                out += val.to_rv()
        return out


@dataclass
class RVC(RV):
    hru_states: HRUStateVariableTableCommand = HRUStateVariableTableCommand()
    basin_states: BasinStateVariablesCommand = BasinStateVariablesCommand()

    @classmethod
    def from_solution(cls, solution: str):
        hru_states = HRUStateVariableTableCommand.parse(solution).hru_states
        basin_states = BasinStateVariablesCommand.parse(solution).basin_states
        return cls(hru_states=hru_states, basin_states=basin_states)


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
    soil_model: SoilModel = None
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
    rain_snow_fraction: RainSnowFraction = RainSnowFraction("RAINSNOW_DATA")
    routing: Routing = Routing("ROUTE_NONE")
    evaporation: Evaporation = None
    ow_evaporation: OW_Evaporation = None
    catchment_route: CatchmentRoute = None
    potential_melt: PotentialMeltMethod = None
    oro_temp_correct: OroTempCorrect = OroTempCorrect("OROCORR_SIMPLELAPSE")
    oro_precip_correct: OroPrecipCorrect = OroPrecipCorrect("OROCORR_SIMPLELAPSE")
    alias: Alias = None
    lake_storage: LakeStorage = None
    hydrologic_processes: HydrologicProcesses = ()
    netcdf_attribute: NetCDFAttribute = None



@dataclass
class RVP(RV):
    params: Any
    soil_classes: SoilClasses = None
    soil_profiles: SoilProfiles = SoilProfiles()
    vegetation_classes: VegetationClasses = VegetationClasses()
    land_use_classes: LandUseClasses = LandUseClasses()
    channel_profiles: ChannelProfiles = ChannelProfiles()
    soil_parameter_list: SoilParameterListCommand = SoilParameterListCommand()
    land_use_parameter_list: LandUseParameterListCommand = LandUseParameterListCommand()
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

    def write(self, workdir: Union[str, Path]):
        """Write configuration files to disk.

        Parameters
        ----------
        workdir: str, Path
          An existing directory where rv files will be written to disk.
        """
        # Check that params has been set
        workdir = Path(workdir)

        for rv in ["rvi", "rvp", "rvc", "rvh", "rvt"]:
            fn = workdir / f"{self.run_name.value}.{rv}"
            fn.write_text(self.to_rv(rv))



@dataclass
class TestConfig(Config):
    calendar: Calendar = Calendar("JULIAN")
