from collections.abc import Sequence
from dataclasses import field, make_dataclass
from typing import Union

from pydantic import Field, field_validator
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
        [(f"X{i:02}", Sym, field(default=Variable(f"X{i:02}"))) for i in range(1, 36)]
        + [(f"R{i:02}", Sym, field(default=Variable(f"R{i:02}"))) for i in range(1, 9)],
    ),
    config=SymConfig,
)

# Bug: RavenC tries to write two variables with the same name `Snow Melt (Liquid) [mm]` to netCDF.


class LandHRU(HRU):
    land_use_class: str = "FOREST"
    veg_class: str = "FOREST"
    soil_profile: str = "DEFAULT_P"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "DEFAULT_T"


class HRUs(rc.HRUs):
    """HRUs command for GR4J.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    root: Sequence[LandHRU]


class Blended(Config):
    """Blended hydrological model for blending outputs from different submodules.

    References
    ----------
    Mai, J., Craig, J. R., and Tolson, B. A.: Simultaneously determining global
    sensitivities of model parameters and model structure, Hydrol. Earth Syst.
    Sci., 24, 5835–5858, https://doi.org/10.5194/hess-24-5835-2020, 2020.

    Chlumsky, R., Mai, J., Craig, J. R., & Tolson, B. A. (2021). Simultaneous
    calibration of hydrologic model structure and parameters using a blended
    model. Water Resources Research, 57, e2020WR029229.
    https://doi.org/10.1029/2020WR029229
    """

    params: P = P()
    hrus: HRUs = Field(
        [
            LandHRU(),
        ],
        alias="HRUs",
    )
    netcdf_attribute: dict[str, str] = {"model_id": "Blended"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    routing: o.Routing = Field("ROUTE_NONE", alias="Routing")
    catchment_route: o.CatchmentRoute = Field("ROUTE_DUMP", alias="CatchmentRouting")
    evaporation: o.Evaporation = Field(o.Evaporation.OUDIN, alias="Evaporation")
    rain_snow_fraction: o.RainSnowFraction = Field(
        o.RainSnowFraction.HBV, alias="RainSnowFraction"
    )
    potential_melt_method: o.PotentialMeltMethod = Field(
        o.PotentialMeltMethod.HMETS, alias="PotentialMeltMethod"
    )
    soil_model: rc.SoilModel = Field(3, alias="SoilModel")

    hydrologic_processes: Sequence[Union[rc.Process, p.ProcessGroup]] = Field(
        [
            p.Precipitation(algo="RAVEN_DEFAULT", source="ATMOS_PRECIP", to="MULTIPLE"),
            p.ProcessGroup(
                p=[
                    p.Infiltration(
                        algo="INF_HMETS", source="PONDED_WATER", to="MULTIPLE"
                    ),
                    p.Infiltration(
                        algo="INF_VIC_ARNO", source="PONDED_WATER", to="MULTIPLE"
                    ),
                    p.Infiltration(
                        algo="INF_HBV", source="PONDED_WATER", to="MULTIPLE"
                    ),
                ],
                params=[P.R01, P.R02],
            ),
            p.Overflow(algo="OVERFLOW_RAVEN", source="SOIL[0]", to="CONVOLUTION[1]"),
            p.ProcessGroup(
                p=[
                    p.Baseflow(
                        algo="BASE_LINEAR_ANALYTIC",
                        source="SOIL[0]",
                        to="SURFACE_WATER",
                    ),
                    p.Baseflow(algo="BASE_VIC", source="SOIL[0]", to="SURFACE_WATER"),
                    p.Baseflow(
                        algo="BASE_TOPMODEL", source="SOIL[0]", to="SURFACE_WATER"
                    ),
                ],
                params=(P.R03, P.R04),
            ),
            p.Percolation(algo="PERC_LINEAR", source="SOIL[0]", to="SOIL[1]"),
            p.Overflow(algo="OVERFLOW_RAVEN", source="SOIL[1]", to="CONVOLUTION[1]"),
            p.Percolation(algo="PERC_LINEAR", source="SOIL[1]", to="SOIL[2]"),
            p.ProcessGroup(
                p=[
                    p.SoilEvaporation(
                        algo="SOILEVAP_ALL", source="SOIL[0]", to="ATMOSPHERE"
                    ),
                    p.SoilEvaporation(
                        algo="SOILEVAP_TOPMODEL", source="SOIL[0]", to="ATMOSPHERE"
                    ),
                ],
                params=(P.R05,),
            ),
            p.Convolve(
                algo="CONVOL_GAMMA", source="CONVOLUTION[0]", to="SURFACE_WATER"
            ),
            p.Convolve(
                algo="CONVOL_GAMMA_2", source="CONVOLUTION[1]", to="SURFACE_WATER"
            ),
            p.ProcessGroup(
                p=[
                    p.Baseflow(
                        algo="BASE_LINEAR_ANALYTIC",
                        source="SOIL[1]",
                        to="SURFACE_WATER",
                    ),
                    p.Baseflow(
                        algo="BASE_POWER_LAW", source="SOIL[1]", to="SURFACE_WATER"
                    ),
                ],
                params=[
                    P.R06,
                ],
            ),
            p.ProcessGroup(
                p=[
                    p.SnowBalance(
                        algo="SNOBAL_HMETS", source="MULTIPLE", to="MULTIPLE"
                    ),
                    p.SnowBalance(
                        algo="SNOBAL_SIMPLE_MELT", source="SNOW", to="PONDED_WATER"
                    ),
                    p.SnowBalance(algo="SNOBAL_HBV", source="MULTIPLE", to="MULTIPLE"),
                ],
                params=(P.R07, P.R08),
            ),
        ],
        alias="HydrologicProcesses",
    )

    soil_classes: SoilClasses = Field(
        [
            {"name": "TOPSOIL"},
            {"name": "PHREATIC"},
            {"name": "DEEP_GW"},
        ],
        alias="SoilClasses",
    )
    vegetation_classes: VegetationClasses = Field(
        [{"name": "FOREST", "max_ht": 4, "max_lai": 5, "max_leaf_cond": 5}],
        alias="VegetationClasses",
    )
    land_use_classes: LandUseClasses = Field(
        [rc.LU(name="FOREST", impermeable_frac=0, forest_coverage=0.02345)],
        alias="LandUseClasses",
    )
    terrain_classes: rc.TerrainClasses = Field(
        [
            rc.TC(
                name="DEFAULT_T",
                hillslope_length=1,
                drainage_density=1,
                topmodel_lambda=P.X07,
            )
        ],
        alias="TerrainClasses",
    )

    soil_profiles: SoilProfiles = Field(
        [
            {"name": "LAKE"},
            {"name": "ROCK"},
            {
                "name": "DEFAULT_P",
                "soil_classes": ("TOPSOIL", "PHREATIC", "DEEP_GW"),
                "thicknesses": (P.X29, P.X30, 1e6),
            },
        ],
        alias="SoilProfiles",
    )

    global_parameter: dict = Field(
        {
            "SNOW_SWI_MIN": P.X13,
            "SNOW_SWI_MAX": P.X14,
            "SWI_REDUCT_COEFF": P.X15,
            "SNOW_SWI": P.X19,
            "RAINSNOW_TEMP": P.X31,
            "RAINSNOW_DELTA": P.X32,
        }
    )

    soil_parameter_list: SoilParameterList = Field(
        {
            "parameters": (
                "POROSITY",
                "PERC_COEFF",
                "PET_CORRECTION",
                "BASEFLOW_COEFF",
                "B_EXP",
                "HBV_BETA",
                "MAX_BASEFLOW_RATE",
                "BASEFLOW_N",
                "FIELD_CAPACITY",
                "SAT_WILT",
            ),
            "pl": [
                PL(
                    name="TOPSOIL",
                    values=(
                        1.0,
                        P.X28,
                        P.X08,
                        P.X04,
                        P.X02,
                        P.X03,
                        P.X05,
                        P.X06,
                        P.X10,
                        P.X09,
                    ),
                ),
                PL(
                    name="PHREATIC",
                    values=(1.0, P.X35, 0, P.X11, 0, 0, 0, P.X12, 0, 0),
                ),
                PL(
                    name="DEEP_GW",
                    values=(1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                ),
            ],
        },
        alias="SoilParameterList",
    )
    land_use_parameter_list: LandUseParameterList = Field(
        {
            "parameters": (
                "MIN_MELT_FACTOR",
                "MAX_MELT_FACTOR",
                "DD_MELT_TEMP",
                "DD_AGGRADATION",
                "REFREEZE_FACTOR",
                "REFREEZE_EXP",
                "DD_REFREEZE_TEMP",
                "HMETS_RUNOFF_COEFF",
                "GAMMA_SHAPE",
                "GAMMA_SCALE",
                "GAMMA_SHAPE2",
                "GAMMA_SCALE2",
                "FOREST_SPARSENESS",
            ),
            "pl": [
                PL(
                    name="[DEFAULT]",
                    values=(
                        P.X24,
                        P.X25,
                        P.X26,
                        P.X27,
                        P.X18,
                        P.X17,
                        P.X16,
                        P.X01,
                        P.X20,
                        P.X21,
                        P.X22,
                        P.X23,
                        0,
                    ),
                )
            ],
        },
        alias="LandUseParameterList",
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        {
            "parameters": ("RAIN_ICEPT_PCT", "SNOW_ICEPT_PCT", "SAI_HT_RATIO"),
            "pl": [PL(name="[DEFAULT]", values=(0, 0, 0))],
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
                    "SOIL[0]": P.X29 * 500,
                    "SOIL[1]": P.X30 * 500,
                },
            ),
        ],
        alias="HRUStateVariableTable",
    )
    _nc_attrs = field_validator("netcdf_attribute")(nc_attrs)

    def __init__(self, **data):
        super().__init__(**data)

        if self.gauge:
            for gauge in self.gauge:
                gauge.rain_correction = self.params.X33
                gauge.snow_correction = self.params.X34
