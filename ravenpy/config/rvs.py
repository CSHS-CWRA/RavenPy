from __future__ import annotations

import datetime as dt
from collections import namedtuple
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Union

import cftime
from pydantic import ConfigDict, Field, ValidationInfo, field_validator, model_validator

from ..config import commands as rc
from ..config import options as o
from ..config import processes as rp
from .base import RV, Sym, optfield, parse_symbolic

"""
Generic Raven model configuration

Note that alias are set to identify class attributes as Raven commands.
"""
date = Union[dt.date, dt.datetime, cftime.datetime]


class RVI(RV):
    # Run parameters
    silent_mode: bool | None = optfield(alias="SilentMode")
    noisy_mode: bool | None = optfield(alias="NoisyMode")

    run_name: str | None = optfield(alias="RunName")
    calendar: o.Calendar | None = optfield(alias="Calendar")
    start_date: date | None = optfield(alias="StartDate")
    assimilation_start_time: date | None = optfield(alias="AssimilationStartTime")
    end_date: date | None = optfield(alias="EndDate")
    duration: float | None = optfield(alias="Duration")
    time_step: float | str | None = optfield(alias="TimeStep")

    # Model description
    routing: o.Routing | None = optfield(alias="Routing")
    catchment_route: o.CatchmentRoute | None = optfield(alias="CatchmentRoute")
    evaporation: o.Evaporation | None = optfield(alias="Evaporation")
    ow_evaporation: o.Evaporation | None = optfield(alias="OW_Evaporation")
    sw_radiation_method: o.SWRadiationMethod | None = optfield(
        alias="SWRadiationMethod"
    )
    sw_cloud_correct: o.SWCloudCorrect | None = optfield(alias="SWCloudCorrect")
    sw_canopy_correct: o.SWCanopyCorrect | None = optfield(alias="SWCanopyCorrect")
    lw_radiation_method: o.LWRadiationMethod | None = optfield(
        alias="LWRadiationMethod"
    )
    windspeed_method: o.WindspeedMethod | None = optfield(alias="WindspeedMethod")
    rain_snow_fraction: o.RainSnowFraction | None = optfield(alias="RainSnowFraction")
    potential_melt_method: o.PotentialMeltMethod | None = optfield(
        alias="PotentialMeltMethod"
    )
    oro_temp_correct: o.OroTempCorrect | None = optfield(alias="OroTempCorrect")
    oro_precip_correct: o.OroPrecipCorrect | None = optfield(alias="OroPrecipCorrect")
    oro_pet_correct: o.OroPETCorrect | None = optfield(alias="OroPETCorrect")
    cloud_cover_method: o.CloudCoverMethod | None = optfield(alias="CloudCoverMethod")
    precip_icept_frac: o.PrecipIceptFract | None = optfield(alias="PrecipIceptFract")
    subdaily_method: o.SubdailyMethod | None = optfield(alias="SubdailyMethod")
    monthly_interpolation_method: o.MonthlyInterpolationMethod | None = optfield(
        alias="MonthlyInterpolationMethod"
    )
    soil_model: rc.SoilModel | None = optfield(alias="SoilModel")
    lake_storage: o.StateVariables | None = optfield(alias="LakeStorage")

    define_hru_groups: Sequence[str] | None = optfield(alias="DefineHRUGroups")

    hydrologic_processes: Sequence[
        rc.Process | rp.Conditional | rp.ProcessGroup
    ] | None = optfield(alias="HydrologicProcesses")
    evaluation_metrics: Sequence[o.EvaluationMetrics] | None = optfield(
        alias="EvaluationMetrics"
    )
    evaluation_period: Sequence[rc.EvaluationPeriod] | None = optfield(
        alias="EvaluationPeriod"
    )
    ensemble_mode: rc.EnsembleMode | None = optfield(alias="EnsembleMode")

    # Options
    write_netcdf_format: bool | None = optfield(alias="WriteNetcdfFormat")
    netcdf_attribute: dict[str, str] | None = optfield(alias="NetCDFAttribute")

    custom_output: Sequence[rc.CustomOutput] | None = optfield(alias="CustomOutput")
    direct_evaporation: bool | None = optfield(
        alias="DirectEvaporation",
        description="Rainfall is automatically reduced through evapotranspiration up to the limit of the calculated PET.",
    )
    deltares_fews_mode: bool | None = optfield(alias="DeltaresFEWSMode")
    debug_mode: bool | None = optfield(alias="DebugMode")
    dont_write_watershed_storage: bool | None = optfield(
        alias="DontWriteWatershedStorage",
        description="Do not write watershed storage variables to disk.",
    )
    pavics_mode: bool | None = optfield(alias="PavicsMode")

    suppress_output: bool | None = optfield(
        alias="SuppressOutput", description="Write minimal output to disk when enabled."
    )
    write_forcing_functions: bool | None = optfield(
        alias="WriteForcingFunctions",
        description="Write watershed averaged forcing functions (e.g. rainfall, radiation, PET, etc).",
    )
    write_subbasin_file: bool | None = optfield(
        alias="WriteSubbasinFile"
    )  # Undocumented

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

            return obj(*v.timetuple()[:6])

    model_config = ConfigDict(arbitrary_types_allowed=True)


class RVT(RV):
    gauge: Sequence[rc.Gauge] | None = optfield(alias="Gauge")
    station_forcing: Sequence[rc.StationForcing] | None = optfield(
        alias="StationForcing"
    )
    gridded_forcing: Sequence[rc.GriddedForcing] | None = optfield(
        alias="GriddedForcing"
    )
    observation_data: Sequence[rc.ObservationData] | None = optfield(
        alias="ObservationData"
    )


class RVP(RV):
    params: Any = None
    soil_classes: rc.SoilClasses | None = optfield(alias="SoilClasses")
    soil_profiles: rc.SoilProfiles | None = optfield(alias="SoilProfiles")
    vegetation_classes: rc.VegetationClasses | None = optfield(
        alias="VegetationClasses"
    )
    land_use_classes: rc.LandUseClasses | None = optfield(alias="LandUseClasses")
    terrain_classes: rc.TerrainClasses | None = optfield(alias="TerrainClasses")
    soil_parameter_list: rc.SoilParameterList | None = optfield(
        alias="SoilParameterList"
    )
    land_use_parameter_list: rc.LandUseParameterList | None = optfield(
        alias="LandUseParameterList"
    )
    vegetation_parameter_list: rc.VegetationParameterList | None = optfield(
        alias="VegetationParameterList"
    )

    channel_profile: Sequence[rc.ChannelProfile] | None = optfield(
        alias="ChannelProfile"
    )

    # TODO: create list of all available parameters to constrain key
    global_parameter: dict[str, Sym] | None = Field({}, alias="GlobalParameter")
    rain_snow_transition: rc.RainSnowTransition | None = optfield(
        alias="RainSnowTransition"
    )
    seasonal_relative_lai: rc.SeasonalRelativeLAI | None = optfield(
        alias="SeasonalRelativeLAI"
    )
    seasonal_relative_height: rc.SeasonalRelativeHeight | None = optfield(
        alias="SeasonalRelativeHeight"
    )


class RVC(RV):
    hru_state_variable_table: rc.HRUStateVariableTable | None = optfield(
        alias="HRUStateVariableTable"
    )
    basin_state_variables: rc.BasinStateVariables | None = optfield(
        alias="BasinStateVariables"
    )
    uniform_initial_conditions: dict[str, Sym] | None = optfield(
        alias="UniformInitialConditions"
    )


class RVH(RV):
    sub_basins: rc.SubBasins | None = optfield(alias="SubBasins")
    sub_basin_group: Sequence[rc.SubBasinGroup] | None = optfield(alias="SubBasinGroup")
    sub_basin_properties: rc.SubBasinProperties | None = optfield(
        alias="SubBasinProperties"
    )
    sb_group_property_multiplier: Sequence[
        rc.SBGroupPropertyMultiplier
    ] | None = optfield(alias="SBGroupPropertyMultiplier")
    hrus: rc.HRUs | None = optfield(alias="HRUs")
    hru_group: Sequence[rc.HRUGroup] | None = optfield(alias="HRUGroup")
    reservoirs: Sequence[rc.Reservoir] | None = optfield(alias="Reservoirs")


class RVE(RV):
    enkf_mode: o.EnKFMode | None = optfield(alias="EnKFMode")
    window_size: int | None = optfield(alias="WindowSize")
    solution_run_name: str | None = optfield(alias="SolutionRunName")
    extra_rvt_filename: str | None = optfield(alias="ExtraRVTFilename")
    output_directory_format: str | Path | None = optfield(alias="OutputDirectoryFormat")
    forecast_rvt_filename: str | None = optfield(alias="ForecastRVTFilename")
    truncate_hindcasts: bool | None = optfield(alias="TruncateHindcasts")
    forcing_perturbation: Sequence[rc.ForcingPerturbation] | None = optfield(
        alias="ForcingPerturbation"
    )
    assimilated_state: Sequence[rc.AssimilatedState] | None = optfield(
        alias="AssimilatedState"
    )
    assimilate_streamflow: Sequence[rc.AssimilateStreamflow] | None = optfield(
        alias="AssimilateStreamflow"
    )
    observation_error_model: Sequence[rc.ObservationErrorModel] | None = optfield(
        alias="ObservationErrorModel"
    )


class Config(RVI, RVC, RVH, RVT, RVP, RVE):
    def header(self, rv):
        """Return the header to print at the top of each RV file."""
        import datetime as dt
        from textwrap import dedent

        import ravenpy

        # TODO: Better mechanism to fetch version
        version = "3.7"

        return dedent(
            f"""
        ###########################################################################################################
        :FileType          {rv.upper()} Raven {version}
        :WrittenBy         RavenPy {ravenpy.__version__} based on setups provided by James Craig and Juliane Mai
        :CreationDate      {dt.datetime.now().isoformat(timespec="seconds")}
        #----------------------------------------------------------------------------------------------------------
        \n\n
        """
        )

    @field_validator("params", mode="before")
    @classmethod
    def _cast_to_dataclass(cls, data):
        """Cast params to a dataclass."""
        # Needed because pydantic v2 does not cast tuples automatically.
        if data is not None and not is_dataclass(data):
            return cls.model_fields["params"].annotation(*data)
        return data

    @field_validator("global_parameter", mode="before")
    @classmethod
    def _update_defaults(cls, v, info: ValidationInfo):
        """Some configuration parameters should be updated with user given arguments, not overwritten."""
        return {**cls.model_fields[info.field_name].default, **v}

    # @model_validator(mode="after")
    # def _parse_symbolic(self):
    #     """If params is numerical, convert symbolic expressions from other fields.
    #     """
    #
    #     if self.params is not None:
    #         p = asdict(self.params)
    #
    #         if not is_symbolic(p):
    #             for key in self.model_fields.keys():
    #                 if key != "params":
    #                     setattr(self, key, parse_symbolic(getattr(self, key), **p))
    #
    #     return self

    def set_params(self, params) -> Config:
        """Return a new instance of Config with params set to their numerical values."""
        # Create params with numerical values
        if not self.is_symbolic:
            raise ValueError(
                "Setting `params` on a configuration without symbolic expressions has no effect."
                "Leave `params` to its default value when instantiating the emulator configuration."
            )

        num_p = self.model_fields["params"].annotation(*params)

        # Parse symbolic expressions using numerical params values
        out = parse_symbolic(self.__dict__, **asdict(num_p))
        out["params"] = num_p

        # Instantiate config class
        # Note: `construct` skips validation. benchmark to see if it speeds things up.
        return self.__class__.model_construct(**out)

    def set_solution(self, fn: Path, timestamp: bool = True) -> Config:
        """Return a new instance of Config with hru, basin states
        and start date set from an existing solution.

        Parameters
        ----------
        fn : Path
          Path to solution file.
        timestamp: bool
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

        # if self.is_symbolic:
        #     raise ValueError(
        #         "Cannot write RV files if `params` has symbolic variables. Use `set_params` method to set numerical "
        #         "values for `params`."
        #     )

        # Get RV class
        rvs = {b.__name__: b for b in Config.__bases__}
        cls = rvs[rv.upper()]

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

    def write_rv(
        self, workdir: str | Path, modelname=None, overwrite=False, header=True
    ):
        """Write configuration files to disk.

        Parameters
        ----------
        workdir: str, Path
          A directory where rv files will be written to disk.
        modelname: str
          File name stem for rv files. If not given, defaults to `RunName` if set, otherwise `raven`.
        overwrite: bool
          If True, overwrite existing configuration files.
        header: bool
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

    def zip(self, workdir: str | Path, modelname=None, overwrite=False):
        """Write configuration to zip file.

        Parameters
        ----------
        workdir: Path, str
          Path to zip archive storing RV files.
        modelname: str
          File name stem for rv files. If not given, defaults to `RunName` if set, otherwise `raven`.
        overwrite: bool
          If True, overwrite existing configuration zip file.
        """
        import zipfile

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


def is_symbolic(params: dict) -> bool:
    """Return True if parameters include a symbolic variable."""
    from dataclasses import is_dataclass

    from pymbolic.primitives import Variable

    if is_dataclass(params):
        params = asdict(params)

    return any([isinstance(v, Variable) for v in params.values()])
