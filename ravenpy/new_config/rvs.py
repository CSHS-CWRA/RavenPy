import datetime as dt
from dataclasses import asdict
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Sequence, Type, TypeVar, Union

import cftime
from pydantic import BaseModel, Extra, Field, root_validator, validator
from pydantic.dataclasses import dataclass

from . import commands as rc
from . import options as o
from .base import RV, Sym, parse_symbolic

"""
Generic Raven model configuration

Note that alias are set to identify class attributes as Raven commands.
"""
date = Union[dt.date, dt.datetime, cftime.datetime]


class RVI(RV):
    # Run parameters
    silent_mode: bool = Field(None, alias="SilentMode")
    noisy_mode: bool = Field(None, alias="NoisyMode")

    run_name: str = Field(None, alias="RunName")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    start_date: date = Field(None, alias="StartDate")
    assimilation_start_time: date = Field(None, alias="AssimilationStartTime")
    end_date: date = Field(None, alias="EndDate")
    duration: float = Field(None, alias="Duration")
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
    lw_radiation_method: o.LWRadiationMethod = Field(None, alias="LWRadiationMethod")
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
    monthly_interpolation_method: o.MonthlyInterpolationMethod = Field(
        None, alias="MonthlyInterpolationMethod"
    )
    soil_model: rc.SoilModel = Field(None, alias="SoilModel")
    lake_storage: o.StateVariables = Field(None, alias="LakeStorage")
    # alias: Dict[str, str] = Field(None, alias="Alias")

    define_hru_groups: Sequence[str] = Field(None, alias="DefineHRUGroups")

    hydrologic_processes: Sequence[
        Union[rc.Process, rc.Conditional, rc.ProcessGroup]
    ] = Field(None, alias="HydrologicProcesses")
    evaluation_metrics: Sequence[o.EvaluationMetrics] = Field(
        None, alias="EvaluationMetrics"
    )

    ensemble_mode: rc.EnsembleMode = Field(None, alias="EnsembleMode")

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
            return rc.SoilModel.parse_obj(v)
        return v

    @validator("start_date", "end_date", "assimilation_start_time")
    def dates2cf(cls, val, values):
        """Convert dates to cftime dates."""
        if val is not None:
            return cftime._cftime.DATE_TYPES[values["calendar"].value.lower()](
                *val.timetuple()[:6]
            )

    class Config:
        arbitrary_types_allowed = True


class RVT(RV):
    gauge: Sequence[rc.Gauge] = Field(None, alias="Gauge", flat=True)
    station_forcing: Sequence[rc.StationForcing] = Field(None, alias="StationForcing")
    gridded_forcing: Sequence[rc.GriddedForcing] = Field(None, alias="GriddedForcing")
    observation_data: Sequence[rc.ObservationData] = Field(
        None, alias="ObservationData"
    )


class RVP(RV):
    params: Any
    soil_classes: rc.SoilClasses = Field(None, alias="SoilClasses")
    soil_profiles: rc.SoilProfiles = Field(None, alias="SoilProfiles")
    vegetation_classes: rc.VegetationClasses = Field(None, alias="VegetationClasses")
    land_use_classes: rc.LandUseClasses = Field(None, alias="LandUseClasses")
    terrain_classes: rc.TerrainClasses = Field(None, alias="TerrainClasses")
    soil_parameter_list: rc.SoilParameterList = Field(None, alias="SoilParameterList")
    land_use_parameter_list: rc.LandUseParameterList = Field(
        None, alias="LandUseParameterList"
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        None, alias="VegetationParameterList"
    )

    channel_profiles: Sequence[rc.ChannelProfile] = ()

    # TODO: create list of all available parameters to constrain key
    global_parameter: Dict[str, Sym] = Field({}, alias="GlobalParameter")
    rain_snow_transition: rc.RainSnowTransition = Field(
        None, alias="RainSnowTransition"
    )
    seasonal_relative_lai: rc.SeasonalRelativeLAI = Field(
        None, alias="SeasonalRelativeLAI"
    )
    seasonal_relative_height: rc.SeasonalRelativeHeight = Field(
        None, alias="SeasonalRelativeLAI"
    )


class RVC(RV):
    hru_states: rc.HRUStateVariableTable = Field(None, alias="HRUStateVariableTable")
    basin_states: rc.BasinStateVariables = Field(None, alias="BasinStateVariables")
    uniform_initial_conditions: Dict[str, Sym] = Field(
        None, alias="UniformInitialConditions"
    )


class RVH(RV):
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    sub_basin_group: Sequence[rc.SubBasinGroup] = ()
    sub_basin_properties: rc.SubBasinProperties = Field(
        None, alias="SubBasinProperties"
    )
    sub_basin_property_multiplier: rc.SBGroupPropertyMultiplierCommand = None
    hrus: rc.HRUs = Field(None, alias="HRUs")
    hru_group: Sequence[rc.HRUGroup] = Field(None, alias="HRUGroup")
    reservoirs: Sequence[rc.Reservoir] = ()


class RVE(RV):
    """Ensemble Kalman filter configuration"""

    output_directory_format: Union[str, Path] = Field(
        None, alias="OutputDirectoryFormat"
    )
    warm_ensemble: Union[str, Path] = Field(
        None,
        alias="WarmEnsemble",
        description="RunName of the simulation whose initial states is used.",
    )
    forecast_rvt_filename: str = Field(None, alias="ForecastRVTFilename")
    truncate_hindcasts: bool = Field(None, alias="TruncateHindcasts")
    forcing_perturbation: Sequence[rc.ForcingPerturbation] = Field(
        None, alias="ForcingPerturbation"
    )
    assimilated_state: Sequence[rc.AssimilatedState] = Field(
        None, alias="AssimilatedState"
    )
    assimilate_streamflow: Sequence[rc.AssimilateStreamflow] = Field(
        None, alias="AssimilateStreamflow"
    )
    observational_error_model: Sequence[rc.ObservationalErrorModel] = Field(
        None, alias="ObservationalErrorModel"
    )


class Config(RVI, RVC, RVH, RVT, RVP, RVE):
    @root_validator
    def assign_symbolic(cls, values):
        """If params is numerical, convert symbolic expressions from other fields."""

        if values["params"] is not None:
            p = asdict(values["params"])

            if not is_symbolic(p):
                return parse_symbolic(values, **p)

        return values

    @validator("global_parameter", pre=True)
    def update_defaults(cls, v, values, config, field):
        """Some configuration parameters should be updated with user given arguments, not overwritten."""
        return {**cls.__fields__[field.name].default, **v}

    def set_params(self, params) -> "Config":
        """Return a new instance of Config with params set to their numerical values."""
        # Get the configuration fields
        dc = self.__dict__.copy()

        # Create params with numerical values
        sym_p = dc.pop("params")
        if not is_symbolic(sym_p):
            raise ValueError(
                "Setting `params` on a configuration without symbolic expressions has no effect."
                "Leave `params` to its default value when instantiating the emulator configuration."
            )

        num_p = sym_p.__class__(*params)

        # Parse symbolic expressions using numerical params values
        out = parse_symbolic(dc, **asdict(num_p))

        # Instantiate config class
        # Note: `construct` skips validation. benchmark to see if it speeds things up.
        return self.__class__.construct(params=num_p, **out)

    def set_solution(self, fn: Path):
        """Return a new instance of Config with hru, basin states and start date set from an existing solution.

        Note that EndDate is set to None to avoid inconsistencies.
        """
        from .parsers import parse_solution

        out = self.__dict__.copy()
        out.update(
            **parse_solution(fn, calendar=self.calendar.value),
            uniform_initial_conditions={},
        )

        return self.__class__(**out)

    def _rv(self, rv: str):
        """Return RV configuration."""

        if self.params and is_symbolic(self.params):
            raise ValueError(
                "Cannot write RV files if `params` has symbolic variables. Use `set_params` method to set numerical "
                "values for `params`."
            )

        # Get RV class
        rvs = {b.__name__: b for b in Config.__bases__}
        cls = rvs[rv.upper()]

        # Instantiate RV class
        attrs = dict(self)
        p = {f: attrs[f] for f in cls.__fields__}
        rv = cls(**p)
        return rv.to_rv()

    @property
    def is_symbolic(self):
        """Return True if configuration contains symbolic expressions."""
        if self.params is not None:
            p = asdict(self.params)
            return is_symbolic(p)

        return False

    @property
    def rvi(self):
        return self._rv("rvi")

    @property
    def rvt(self):
        return self._rv("rvt")

    @property
    def rvp(self):
        return self._rv("rvp")

    @property
    def rvc(self):
        return self._rv("rvc")

    @property
    def rvh(self):
        return self._rv("rvh")

    @property
    def rve(self):
        return self._rv("rve")

    def write_rv(self, workdir: Union[str, Path], modelname=None, overwrite=False):
        """Write configuration files to disk.

        Parameters
        ----------
        workdir: str, Path
          A directory where rv files will be written to disk.
        modelname: str
          File name stem for rv files. If not given, defaults to `RunName` if set, otherwise `raven`.
        overwrite: bool
          If True, overwrite existing configuration files.
        """
        workdir = Path(workdir)
        if not workdir.exists():
            workdir.mkdir(parents=True)

        mandatory = ["rvi", "rvp", "rvc", "rvh", "rvt"]
        optional = ["rve"]

        if modelname is None:
            modelname = self.run_name or "raven"

        out = {}
        for rv in mandatory + optional:
            fn = workdir / f"{modelname}.{rv}"
            if fn.exists() and not overwrite:
                raise OSError(f"{fn} already exists and would be overwritten.")

            text = self._rv(rv)
            if rv in mandatory or text.strip():
                fn.write_text(text)
                out[rv] = fn

        return out

    def zip(self, path):
        """Write configuration to zip file.

        Parameters
        ----------
        path: Path, str
          Path to zip archive storing RV files.
        """
        import zipfile

        with zipfile.ZipFile(path, "w") as fh:
            for rv in ["rvi", "rvp", "rvc", "rvh", "rvt"]:
                fn = self.rv_fn(rv)
                fh.write(fn, self._rv(rv))
        return zip


def is_symbolic(params: Dict) -> bool:
    """Return True if parameters include a symbolic variable."""
    from dataclasses import is_dataclass

    from pymbolic.primitives import Variable

    if is_dataclass(params):
        params = asdict(params)

    return any([isinstance(v, Variable) for v in params.values()])
