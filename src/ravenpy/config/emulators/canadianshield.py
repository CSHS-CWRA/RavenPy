from collections.abc import Sequence
from dataclasses import field, make_dataclass
from typing import Union

from pydantic import Field, field_validator, model_validator
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

import ravenpy.config.processes as p
from ravenpy.config import commands as rc
from ravenpy.config import options as o
from ravenpy.config.base import Sym, SymConfig
from ravenpy.config.commands import (
    HRU,
    PL,
    LandUseClasses,
    LandUseParameterList,
    SoilClasses,
    SoilParameterList,
    SoilProfiles,
    VegetationClasses,
)
from ravenpy.config.defaults import nc_attrs
from ravenpy.config.rvs import Config

P = dataclass(
    make_dataclass(
        "Params",
        [(f"X{i:02}", Sym, field(default=Variable(f"X{i:02}"))) for i in range(1, 35)],
    ),
    config=SymConfig,
)

# Bug: RavenC tries to write two variables with the same name `Snow Melt (Liquid) [mm]` to netCDF.


class OrganicHRU(HRU):
    hru_id: int = 1
    land_use_class: str = "FOREST"
    veg_class: str = "FOREST"
    soil_profile: str = "SOILP_ORG"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"


class BedRockHRU(HRU):
    hru_id: int = 2
    land_use_class: str = "FOREST"
    veg_class: str = "FOREST"
    soil_profile: str = "SOILP_BEDROCK"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"


class HRUs(rc.HRUs):
    """HRUs command for CanadianShield."""

    root: tuple[OrganicHRU, BedRockHRU]


class CanadianShield(Config):
    """The Canadian Shield model is a useful configuration of Raven for
    Canadian shield basins characterized by shallow soils atop rock,
    with ample exposed rock and lakes. It is a custom model for this
    type of region, but there is no reference to it in the literature.
    """

    params: P = P()
    hrus: HRUs = Field([OrganicHRU(), BedRockHRU()], alias="HRUs")
    netcdf_attribute: dict[str, str] = {"model_id": "CanadianShield"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    land_use_classes: LandUseClasses = Field(
        [rc.LU(name="FOREST", impermeable_frac=0.0, forest_coverage=0.02345)],
        alias="LandUseClasses",
    )

    routing: o.Routing = Field("ROUTE_NONE", alias="Routing")
    catchment_route: o.CatchmentRoute = Field(
        "ROUTE_TRI_CONVOLUTION", alias="CatchmentRouting"
    )
    evaporation: o.Evaporation = Field(
        o.Evaporation.HARGREAVES_1985, alias="Evaporation"
    )
    ow_evaporation: o.Evaporation = Field(
        o.Evaporation.HARGREAVES_1985, alias="OW_Evaporation"
    )
    sw_canopy_correct: o.SWCanopyCorrect = Field(
        o.SWCanopyCorrect.STATIC, alias="SWCanopyCorrect"
    )
    rain_snow_fraction: o.RainSnowFraction = Field(
        o.RainSnowFraction.DINGMAN, alias="RainSnowFraction"
    )
    potential_melt_method: o.PotentialMeltMethod = Field(
        o.PotentialMeltMethod.DEGREE_DAY, alias="PotentialMeltMethod"
    )
    precip_icept_frac: o.PrecipIceptFract = Field(
        o.PrecipIceptFract.LAI, alias="PrecipIceptFrac"
    )
    soil_model: rc.SoilModel = Field(3, alias="SoilModel")
    monthly_interpolation_method: o.MonthlyInterpolationMethod = Field(
        o.MonthlyInterpolationMethod.LINEAR_MID, alias="MonthlyInterpolationMethod"
    )
    hydrologic_processes: Sequence[Union[rc.Process, p.Conditional]] = Field(
        [
            p.SnowRefreeze(algo="FREEZE_DEGREE_DAY", source="SNOW_LIQ", to="SNOW"),
            p.Precipitation(algo="PRECIP_RAVEN", source="ATMOS_PRECIP", to="MULTIPLE"),
            p.CanopyEvaporation(
                algo="CANEVP_MAXIMUM", source="CANOPY", to="ATMOSPHERE"
            ),
            p.CanopySublimation(
                algo="CANEVP_MAXIMUM", source="CANOPY_SNOW", to="ATMOSPHERE"
            ),
            p.SnowBalance(algo="SNOBAL_TWO_LAYER", source="MULTIPLE", to="MULTIPLE"),
            p.Abstraction(algo="ABST_FILL", source="PONDED_WATER", to="DEPRESSION"),
            p.OpenWaterEvaporation(
                algo="OPEN_WATER_EVAP", source="DEPRESSION", to="ATMOSPHERE"
            ),
            p.Infiltration(algo="INF_HBV", source="PONDED_WATER", to="MULTIPLE"),
            p.Interflow(algo="INTERFLOW_PRMS", source="SOIL[0]", to="SURFACE_WATER"),
            p.Baseflow(algo="BASE_POWER_LAW", source="SOIL[1]", to="SURFACE_WATER"),
            p.Conditional(kind="HRU_TYPE", op="IS", value="SOIL_ORG"),
            p.Baseflow(algo="BASE_POWER_LAW", source="SOIL[2]", to="SURFACE_WATER"),
            p.Percolation(algo="PERC_GAWSER", source="SOIL[0]", to="SOIL[1]"),
            p.Conditional(kind="HRU_TYPE", op="IS", value="SOIL_ORG"),
            p.Percolation(algo="PERC_GAWSER", source="SOIL[1]", to="SOIL[2]"),
            p.Conditional(kind="HRU_TYPE", op="IS", value="SOIL_ORG"),
            p.Percolation(algo="PERC_GAWSER", source="SOIL[0]", to="SOIL[2]"),
            p.Conditional(kind="HRU_TYPE", op="IS", value="SOIL_BEDROCK"),
            p.SoilEvaporation(algo="SOILEVAP_ROOT", source="SOIL[0]", to="ATMOSPHERE"),
        ],
        alias="HydrologicProcesses",
    )

    soil_classes: SoilClasses = Field(
        [{"name": "TOPSOIL"}, {"name": "VADOSE"}, {"name": "FRACBEDROCK"}],
        alias="SoilClasses",
    )
    soil_profiles: SoilProfiles = Field(
        [
            {"name": "LAKE"},
            {"name": "ROCK"},
            {
                "name": "SOILP_ORG",
                "soil_classes": ("TOPSOIL", "VADOSE", "FRACBEDROCK"),
                "thicknesses": (P.X01, P.X02, P.X03),
            },
            {
                "name": "SOILP_BEDROCK",
                "soil_classes": ("TOPSOIL", "VADOSE", "FRACBEDROCK"),
                "thicknesses": (P.X01, 0, P.X03),
            },
        ],
        alias="SoilProfiles",
    )
    vegetation_classes: VegetationClasses = Field(
        [{"name": "FOREST", "max_ht": 5, "max_lai": 5, "max_leaf_cond": 5}],
        alias="VegetationClasses",
    )

    global_parameter: dict = Field(
        {
            "SNOW_SWI": P.X15,
            "SNOW_SWI_MIN": P.X16,
            "SNOW_SWI_MAX": P.X17,
            "SWI_REDUCT_COEFF": P.X18,
            "RAINSNOW_TEMP": P.X19,
            "RAINSNOW_DELTA": P.X20,
            "MAX_SWE_SURFACE": P.X21,
            "TOC_MULTIPLIER": P.X22,
        }
    )
    soil_parameter_list: SoilParameterList = Field(
        {
            "parameters": (
                "POROSITY",
                "HBV_BETA",
                "BASEFLOW_COEFF",
                "BASEFLOW_N",
                "MAX_INTERFLOW_RATE",
                "FIELD_CAPACITY",
                "SAT_WILT",
                "MAX_PERC_RATE",
                "PET_CORRECTION",
            ),
            "pl": [
                PL(
                    name="TOPSOIL",
                    values=(1.0, P.X07, 0.0, 0.0, P.X12, P.X06, P.X05, P.X13, P.X04),
                ),
                PL(name="VADOSE", values=(1.0, 0.0, P.X08, P.X10, 0, 0, 0, P.X14, 0)),
                PL(name="FRACBEDROCK", values=(1.0, 0.0, P.X09, P.X11, 0, 0, 0, 0, 0)),
            ],
        },
        alias="SoilParameterList",
    )
    land_use_parameter_list: LandUseParameterList = Field(
        {
            "parameters": (
                "FOREST_SPARSENESS",
                "MELT_FACTOR",
                "DD_MELT_TEMP",
                "REFREEZE_FACTOR",
                "DEP_MAX",
                "OW_PET_CORR",
            ),
            "pl": [PL(name="FOREST", values=(0, P.X25, P.X24, P.X23, P.X26, P.X27))],
        },
        alias="LandUseParameterList",
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        {
            "parameters": (
                "SVF_EXTINCTION",
                "SAI_HT_RATIO",
                "RAIN_ICEPT_FACT",
                "SNOW_ICEPT_FACT",
                "MAX_CAPACITY",
                "MAX_SNOW_CAPACITY",
            ),
            "pl": [PL(name="FOREST", values=(0.5, 1, P.X28, P.X29, P.X30, P.X31))],
        },
        alias="VegetationParameterList",
    )
    seasonal_relative_lai: rc.SeasonalRelativeLAI = Field(
        rc.SeasonalRelativeLAI(), alias="SeasonalRelativeLAI"
    )
    seasonal_relative_height: rc.SeasonalRelativeHeight = Field(
        rc.SeasonalRelativeHeight(), alias="SeasonalRelativeHeight"
    )

    hru_state_variable_table: rc.HRUStateVariableTable = Field(
        [
            rc.HRUState(
                hru_id=1,
                data={
                    "SOIL[0]": P.X01 * 500,
                    "SOIL[1]": P.X02 * 500,
                    "SOIL[2]": P.X03 * 500,
                },
            ),
            rc.HRUState(
                hru_id=2,
                data={
                    "SOIL[0]": P.X01 * 500,
                    "SOIL[1]": P.X02 * 500,
                    "SOIL[2]": P.X03 * 500,
                },
            ),
        ],
        alias="HRUStateVariableTable",
    )

    @field_validator("hrus", mode="before")
    @classmethod
    def equal_area(cls, hrus):
        if len(hrus) != 2:
            raise ValueError(
                f"CanadianShield must have two HRUs, one organic and one bedrock, has {len(hrus)}."
            )

        areas = []
        for hru in hrus:
            if isinstance(hru, dict):
                areas.append(hru["area"])
            else:
                areas.append(hru.area)

        if len(set(areas)) > 1:
            raise ValueError(f"HRUs must have equal areas: {areas}")

        return hrus

    _nc_attrs = field_validator("netcdf_attribute")(nc_attrs)
    _initialized: bool = False

    @model_validator(mode="after")
    def _init(self):
        if not self._initialized:
            hrus = self.hrus.root
            hrus[0].area *= self.params.X34
            hrus[1].area *= 1 - self.params.X34

            if self.gauge:
                for gauge in self.gauge:
                    gauge.rain_correction = self.params.X32
                    gauge.snow_correction = self.params.X33

            self._initialized = True

        return self
