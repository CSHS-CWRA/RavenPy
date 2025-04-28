import datetime as dt
import zipfile
from collections.abc import Sequence
from dataclasses import asdict, fields, is_dataclass
from pathlib import Path
from textwrap import dedent
from typing import Any, Optional, Union

import cftime
from pydantic import ConfigDict, Field, ValidationInfo, field_validator

from ..config import commands as rc
from ..config import options as o
from ..config import processes as rp
from .base import RV, Sym, optfield, parse_symbolic

try:
    from raven_hydro import __raven_version__
except ImportError:
    __raven_version__ = "0.0.0"


"""
Generic Raven model configuration.

Note that aliases are set to identify class attributes as Raven commands.
"""
date = Union[dt.date, dt.datetime, cftime.datetime]


class RVI(RV):
    __rv__ = "RVI"

    # Run parameters
    silent_mode: Optional[bool] = optfield(alias="SilentMode")
    noisy_mode: Optional[bool] = optfield(alias="NoisyMode")

    run_name: Optional[str] = optfield(alias="RunName")
    calendar: Optional[o.Calendar] = optfield(alias="Calendar")
    start_date: Optional[Union[str, date]] = optfield(alias="StartDate")
    assimilation_start_time: Optional[date] = optfield(alias="AssimilationStartTime")
    end_date: Optional[Union[str, date]] = optfield(alias="EndDate")
    duration: Optional[float] = optfield(alias="Duration")
    time_step: Optional[Union[float, str]] = optfield(alias="TimeStep")
    interpolation: Optional[o.Interpolation] = optfield(alias="Interpolation")

    # Model description
    routing: Optional[o.Routing] = optfield(alias="Routing")
    catchment_route: Optional[o.CatchmentRoute] = optfield(alias="CatchmentRoute")
    evaporation: Optional[o.Evaporation] = optfield(alias="Evaporation")
    ow_evaporation: Optional[o.Evaporation] = optfield(alias="OW_Evaporation")
    sw_radiation_method: Optional[o.SWRadiationMethod] = optfield(
        alias="SWRadiationMethod"
    )
    sw_cloud_correct: Optional[o.SWCloudCorrect] = optfield(alias="SWCloudCorrect")
    sw_canopy_correct: Optional[o.SWCanopyCorrect] = optfield(alias="SWCanopyCorrect")
    lw_radiation_method: Optional[o.LWRadiationMethod] = optfield(
        alias="LWRadiationMethod"
    )
    windspeed_method: Optional[o.WindspeedMethod] = optfield(alias="WindspeedMethod")
    rain_snow_fraction: Optional[o.RainSnowFraction] = optfield(
        alias="RainSnowFraction"
    )
    potential_melt_method: Optional[o.PotentialMeltMethod] = optfield(
        alias="PotentialMeltMethod"
    )
    oro_temp_correct: Optional[o.OroTempCorrect] = optfield(alias="OroTempCorrect")
    oro_precip_correct: Optional[o.OroPrecipCorrect] = optfield(
        alias="OroPrecipCorrect"
    )
    oro_pet_correct: Optional[o.OroPETCorrect] = optfield(alias="OroPETCorrect")
    cloud_cover_method: Optional[o.CloudCoverMethod] = optfield(
        alias="CloudCoverMethod"
    )
    precip_icept_frac: Optional[o.PrecipIceptFract] = optfield(alias="PrecipIceptFract")
    subdaily_method: Optional[o.SubdailyMethod] = optfield(alias="SubdailyMethod")
    monthly_interpolation_method: Optional[o.MonthlyInterpolationMethod] = optfield(
        alias="MonthlyInterpolationMethod"
    )
    soil_model: Optional[rc.SoilModel] = optfield(alias="SoilModel")
    temperature_correction: Optional[bool] = optfield(
        alias="TemperatureCorrection",
        description="Gridded or gauged temperature bias correction.",
    )
    lake_storage: Optional[o.StateVariables] = optfield(alias="LakeStorage")
    relative_humidity_method: Optional[o.RelativeHumidityMethod] = optfield(
        alias="RelativeHumidityMethod"
    )

    define_hru_groups: Optional[Sequence[str]] = optfield(alias="DefineHRUGroups")

    hydrologic_processes: Optional[
        Sequence[Union[rc.Process, rp.Conditional, rp.ProcessGroup]]
    ] = optfield(alias="HydrologicProcesses")
    evaluation_metrics: Optional[Sequence[o.EvaluationMetrics]] = optfield(
        alias="EvaluationMetrics"
    )
    evaluation_period: Optional[Sequence[rc.EvaluationPeriod]] = optfield(
        alias="EvaluationPeriod"
    )
    ensemble_mode: Optional[rc.EnsembleMode] = optfield(alias="EnsembleMode")

    # Options
    write_netcdf_format: Optional[bool] = optfield(alias="WriteNetcdfFormat")
    netcdf_attribute: Optional[dict[str, str]] = optfield(alias="NetCDFAttribute")

    custom_output: Optional[Sequence[rc.CustomOutput]] = optfield(alias="CustomOutput")
    direct_evaporation: Optional[bool] = optfield(
        alias="DirectEvaporation",
        description="Rainfall is automatically reduced through evapotranspiration up to the limit of the calculated PET.",
    )
    deltares_fews_mode: Optional[bool] = optfield(alias="DeltaresFEWSMode")
    debug_mode: Optional[bool] = optfield(alias="DebugMode")
    dont_write_watershed_storage: Optional[bool] = optfield(
        alias="DontWriteWatershedStorage",
        description="Do not write watershed storage variables to disk.",
    )
    pavics_mode: Optional[bool] = optfield(alias="PavicsMode")

    suppress_output: Optional[bool] = optfield(
        alias="SuppressOutput", description="Write minimal output to disk when enabled."
    )
    write_forcing_functions: Optional[bool] = optfield(
        alias="WriteForcingFunctions",
        description="Write watershed averaged forcing functions (e.g. rainfall, radiation, PET, etc).",
    )
    write_subbasin_file: Optional[bool] = optfield(
        alias="WriteSubbasinFile"
    )  # Undocumented
    write_local_flows: Optional[bool] = optfield(
        alias="WriteLocalFlows",
        description="Write local contribution to hydrograph in hydrograph.csv",
    )

    @field_validator("soil_model", mode="before")
    @classmethod
    def init_soil_model(cls, v):
        if isinstance(v, int):
            return rc.SoilModel(v)
        return v

    @field_validator("start_date", "end_date", "assimilation_start_time")
    @classmethod
    def dates2cf(cls, v, info):
        """Convert dates to cftime dates."""
        if v is not None:
            calendar = (
                info.data.get("calendar") or o.Calendar.PROLEPTIC_GREGORIAN
            ).value.lower()

            obj = cftime._cftime.DATE_TYPES[calendar]
            if isinstance(v, str):
                v = dt.datetime.fromisoformat(v)

            return obj(*v.timetuple()[:6])

    model_config = ConfigDict(arbitrary_types_allowed=True)


class RVT(RV):
    __rv__ = "RVT"

    gauge: Optional[Sequence[rc.Gauge]] = optfield(alias="Gauge")
    station_forcing: Optional[Sequence[rc.StationForcing]] = optfield(
        alias="StationForcing"
    )
    gridded_forcing: Optional[Sequence[rc.GriddedForcing]] = optfield(
        alias="GriddedForcing"
    )
    observation_data: Optional[Sequence[rc.ObservationData]] = optfield(
        alias="ObservationData"
    )


class RVP(RV):
    __rv__ = "RVP"

    params: Any = None
    soil_classes: Optional[rc.SoilClasses] = optfield(alias="SoilClasses")
    soil_profiles: Optional[rc.SoilProfiles] = optfield(alias="SoilProfiles")
    vegetation_classes: Optional[rc.VegetationClasses] = optfield(
        alias="VegetationClasses"
    )
    land_use_classes: Optional[rc.LandUseClasses] = optfield(alias="LandUseClasses")
    terrain_classes: Optional[rc.TerrainClasses] = optfield(alias="TerrainClasses")
    soil_parameter_list: Optional[rc.SoilParameterList] = optfield(
        alias="SoilParameterList"
    )
    land_use_parameter_list: Optional[rc.LandUseParameterList] = optfield(
        alias="LandUseParameterList"
    )
    vegetation_parameter_list: Optional[rc.VegetationParameterList] = optfield(
        alias="VegetationParameterList"
    )

    channel_profile: Optional[Sequence[rc.ChannelProfile]] = optfield(
        alias="ChannelProfile"
    )

    # TODO: create list of all available parameters to constrain key
    global_parameter: Optional[dict[str, Sym]] = Field({}, alias="GlobalParameter")
    rain_snow_transition: Optional[rc.RainSnowTransition] = optfield(
        alias="RainSnowTransition"
    )
    seasonal_relative_lai: Optional[rc.SeasonalRelativeLAI] = optfield(
        alias="SeasonalRelativeLAI"
    )
    seasonal_relative_height: Optional[rc.SeasonalRelativeHeight] = optfield(
        alias="SeasonalRelativeHeight"
    )


class RVC(RV):
    __rv__ = "RVC"

    hru_state_variable_table: Optional[rc.HRUStateVariableTable] = optfield(
        alias="HRUStateVariableTable"
    )
    basin_state_variables: Optional[rc.BasinStateVariables] = optfield(
        alias="BasinStateVariables"
    )
    uniform_initial_conditions: Optional[dict[str, Sym]] = optfield(
        alias="UniformInitialConditions"
    )


class RVH(RV):
    __rv__ = "RVH"

    sub_basins: Optional[rc.SubBasins] = optfield(alias="SubBasins")
    sub_basin_group: Optional[Sequence[rc.SubBasinGroup]] = optfield(
        alias="SubBasinGroup"
    )
    sub_basin_properties: Optional[rc.SubBasinProperties] = optfield(
        alias="SubBasinProperties"
    )
    sb_group_property_multiplier: Optional[Sequence[rc.SBGroupPropertyMultiplier]] = (
        optfield(alias="SBGroupPropertyMultiplier")
    )
    hrus: Optional[rc.HRUs] = optfield(alias="HRUs")
    hru_group: Optional[Sequence[rc.HRUGroup]] = optfield(alias="HRUGroup")
    reservoirs: Optional[Sequence[rc.Reservoir]] = optfield(alias="Reservoirs")


class RVE(RV):
    __rv__ = "RVE"

    enkf_mode: Optional[o.EnKFMode] = optfield(alias="EnKFMode")
    window_size: Optional[int] = optfield(alias="WindowSize")
    solution_run_name: Optional[str] = optfield(alias="SolutionRunName")
    extra_rvt_filename: Optional[str] = optfield(alias="ExtraRVTFilename")
    output_directory_format: Optional[Union[str, Path]] = optfield(
        alias="OutputDirectoryFormat"
    )
    forecast_rvt_filename: Optional[str] = optfield(alias="ForecastRVTFilename")
    truncate_hindcasts: Optional[bool] = optfield(alias="TruncateHindcasts")
    forcing_perturbation: Optional[Sequence[rc.ForcingPerturbation]] = optfield(
        alias="ForcingPerturbation"
    )
    assimilated_state: Optional[Sequence[rc.AssimilatedState]] = optfield(
        alias="AssimilatedState"
    )
    assimilate_streamflow: Optional[Sequence[rc.AssimilateStreamflow]] = optfield(
        alias="AssimilateStreamflow"
    )
    observation_error_model: Optional[Sequence[rc.ObservationErrorModel]] = optfield(
        alias="ObservationErrorModel"
    )


class Config(RVI, RVC, RVH, RVT, RVP, RVE):
    __rv__ = None

    @staticmethod
    def header(rv):
        """Return the header to print at the top of each RV file."""
        import ravenpy

        return dedent(
            f"""
        ###########################################################################################################
        :FileType          {rv.upper()} Raven {__raven_version__}
        :WrittenBy         RavenPy {ravenpy.__version__} based on setups provided by James Craig and Juliane Mai
        :CreationDate      {dt.datetime.now().isoformat(timespec="seconds")}
        #----------------------------------------------------------------------------------------------------------
        \n\n
        """
        )

    @field_validator("params", mode="before")
    @classmethod
    def _cast_to_dataclass(cls, data: Union[dict, Sequence]):
        """Cast params to a dataclass."""
        # Needed because pydantic v2 does not cast tuples automatically.
        if data is not None:
            if is_dataclass(data):
                return data

            if isinstance(data, dict):
                return cls.model_fields["params"].annotation(**data)

            return cls.model_fields["params"].annotation(*data)

    @field_validator("global_parameter", mode="before")
    @classmethod
    def _update_defaults(cls, v, info: ValidationInfo):
        """Some configuration parameters should be updated with user given arguments, not overwritten."""
        return {**cls.model_fields[info.field_name].default, **v}

    def set_params(self, params: Union[dict, Sequence]) -> "Config":
        """Return a new instance of Config with params frozen to their numerical values."""
        # Create params with numerical values
        if not self.is_symbolic:
            raise ValueError(
                "Setting `params` on a configuration without symbolic expressions has no effect. "
                "Leave `params` to its default value when instantiating the emulator configuration."
            )

        p = self.__class__._cast_to_dataclass(params)

        # Parse symbolic expressions using numerical params values
        out = parse_symbolic(self.__dict__, **asdict(p))
        out["params"] = p

        # Instantiate config class
        return self.__class__.model_construct(**out)

    def set_solution(self, fn: Path, timestamp: bool = True) -> "Config":
        """
        Return a new instance of Config with hru, basin states and start date set from an existing solution.

        Parameters
        ----------
        fn : Path
           Path to solution file.
        timestamp : bool
           If False, ignore time stamp information in the solution. If True, the solution
           will set StartDate to the solution's timestamp.

        Returns
        -------
        Config
            Config with internal state set from the solution file.
        """
        from .defaults import CALENDAR
        from .parsers import parse_solution

        try:
            calendar = self.calendar.value
        except AttributeError:
            calendar = CALENDAR

        out = self.__dict__.copy()
        sol = parse_solution(fn, calendar=calendar)
        if timestamp is False:
            sol.pop("start_date")

        out.update(**sol, uniform_initial_conditions={})

        return self.__class__(**out)

    def duplicate(self, **kwds):
        """Duplicate this model, changing the values given in the keywords."""
        out = self.model_copy(deep=True)
        for key, val in (
            self.model_validate(kwds).model_dump(exclude_unset=True).items()
        ):
            setattr(out, key, val)
        return out

    def _rv(self, rv: str):
        """Return RV configuration."""
        import inspect

        # if self.is_symbolic:
        #     raise ValueError(
        #         "Cannot write RV files if `params` has symbolic variables. Use `set_params` method to set numerical "
        #         "values for `params`."
        #     )
        # Get RV class
        for cls in inspect.getmro(self.__class__):
            if getattr(cls, "__rv__") == rv.upper():
                break

        # Instantiate RV class
        attrs = dict(self)
        rv_attrs = {f: attrs[f] for f in cls.model_fields}

        # Get model parameters and convert symbolic expressions.
        if self.params is not None:
            p = asdict(self.params)
            rv_attrs = parse_symbolic(rv_attrs, **p)

        rv = cls.model_validate(rv_attrs)
        return rv.to_rv()

    @property
    def is_symbolic(self):
        """Return True if configuration contains symbolic expressions."""
        if self.params is not None:
            p = {
                field.name: getattr(self.params, field.name)
                for field in fields(self.params)
            }
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

    def write_rv(
        self,
        workdir: Union[str, Path],
        modelname: Optional[str] = None,
        overwrite: bool = False,
        header: bool = True,
    ):
        """
        Write configuration files to disk.

        Parameters
        ----------
        workdir : str, Path
            A directory where rv files will be written to disk.
        modelname : str
            File name stem for rv files. If not given, defaults to `RunName` if set, otherwise `raven`.
        overwrite : bool
            If True, overwrite existing configuration files.
        header : bool
            If True, write a header at the top of each RV file.
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

            # RV header
            text = self.header(rv) if header else ""

            # RV content
            text += self._rv(rv)

            # Write to disk
            if rv in mandatory or text.strip():
                fn.write_text(text)
                out[rv] = fn

        return out

    def zip(
        self,
        workdir: Union[str, Path],
        modelname: Optional[str] = None,
        overwrite: bool = False,
    ):
        """
        Write configuration to zip file.

        Parameters
        ----------
        workdir : Path, str
            Path to zip archive storing RV files.
        modelname : str, optional
            File name stem for rv files. If not given, defaults to `RunName` if set, otherwise `raven`.
        overwrite : bool
            If True, overwrite existing configuration zip file.
        """
        workdir = Path(workdir)
        if not workdir.exists():
            workdir.mkdir(parents=True)

        if modelname is None:
            modelname = self.run_name or "raven"

        fn = workdir / f"{modelname}.zip"
        if fn.exists() and not overwrite:
            raise OSError(f"{fn} already exists and would be overwritten.")

        with zipfile.ZipFile(fn, "w") as fh:
            for rv in ["rvi", "rvp", "rvc", "rvh", "rvt", "rve"]:
                f = f"{modelname}.{rv}"
                txt = self._rv(rv)
                if txt.strip():
                    fh.writestr(f, txt)
        return zip


def is_symbolic(params: Union[dict, Any]) -> bool:
    """Return True if parameters include a symbolic variable."""
    from pymbolic.primitives import Variable

    if is_dataclass(params) and not isinstance(params, type):
        params = {field.name: getattr(params, field.name) for field in fields(params)}

    return any([isinstance(v, Variable) for v in params.values()])
