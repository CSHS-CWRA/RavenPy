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
        [(f"X{i:02}", Sym, field(default=Variable(f"X{i:02}"))) for i in range(1, 22)],
    ),
    config=SymConfig,
)

# Bug: RavenC tries to write two variables with the same name `Snow Melt (Liquid) [mm]` to netCDF.


class LandHRU(HRU):
    land_use_class: str = "FOREST"
    veg_class: str = "FOREST"
    soil_profile: str = "DEFAULT_P"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"


class HRUs(rc.HRUs):
    """HRUs command for GR4J.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    root: Sequence[LandHRU]


class SACSMA(Config):
    """Sacramento - Soil Moisture Accounting model

    References
    ----------
    Sorooshian, S., Duan, Q., and Gupta, V. K. (1993), Calibration of rainfall-runoff models:
    Application of global optimization to the Sacramento Soil Moisture Accounting Model,
    Water Resour. Res., 29( 4), 1185– 1194, doi:10.1029/92WR02617.
    """

    params: P = P()
    hrus: HRUs = Field([LandHRU()], alias="HRUs")
    netcdf_attribute: dict[str, str] = {"model_id": "SACSMA"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    routing: o.Routing = Field("ROUTE_NONE", alias="Routing")
    catchment_route: o.CatchmentRoute = Field(
        "ROUTE_GAMMA_CONVOLUTION", alias="CatchmentRouting"
    )
    evaporation: o.Evaporation = Field(o.Evaporation.OUDIN, alias="Evaporation")
    rain_snow_fraction: o.RainSnowFraction = Field(
        o.RainSnowFraction.DATA, alias="RainSnowFraction"
    )
    potential_melt_method: o.PotentialMeltMethod = Field(
        o.PotentialMeltMethod.DEGREE_DAY, alias="PotentialMeltMethod"
    )
    soil_model: rc.SoilModel = Field(7, alias="SoilModel")

    hydrologic_processes: Sequence[Union[rc.Process, p.Conditional]] = Field(
        [
            p.SnowBalance(algo="SNOBAL_SIMPLE_MELT", source="SNOW", to="PONDED_WATER"),
            p.Precipitation(algo="RAVEN_DEFAULT", source="ATMOS_PRECIP", to="MULTIPLE"),
            p.SoilEvaporation(
                algo="SOILEVAP_SACSMA", source="MULTIPLE", to="ATMOSPHERE"
            ),
            p.SoilBalance(algo="SOILBAL_SACSMA", source="MULTIPLE", to="MULTIPLE"),
            p.OpenWaterEvaporation(
                algo="OPEN_WATER_RIPARIAN", source="SURFACE_WATER", to="ATMOSPHERE"
            ),
        ],
        alias="HydrologicProcesses",
    )

    soil_classes: SoilClasses = Field(
        [
            {"name": "UZT"},
            {"name": "UZF"},
            {"name": "LZT"},
            {"name": "LZFP"},
            {"name": "LZFS"},
            {"name": "ADIM"},
            {"name": "GW"},
        ],
        alias="SoilClasses",
    )
    vegetation_classes: VegetationClasses = Field(
        [{"name": "FOREST", "max_ht": 4, "max_lai": 5, "max_leaf_cond": 5}],
        alias="VegetationClasses",
    )
    land_use_classes: LandUseClasses = Field(
        [rc.LU(name="FOREST", impermeable_frac=P.X13, forest_coverage=1.0)],
        alias="LandUseClasses",
    )

    soil_profiles: SoilProfiles = Field(
        [
            {"name": "LAKE"},
            {
                "name": "DEFAULT_P",
                "soil_classes": ("UZT", "UZF", "LZT", "LZFP", "LZFS", "ADIM", "GW"),
                "thicknesses": (P.X04, P.X05, P.X06, P.X08, P.X07, 100, 100),
            },
        ],
        alias="SoilProfiles",
    )

    soil_parameter_list: SoilParameterList = Field(
        {
            "parameters": (
                "POROSITY",
                "SAC_PERC_ALPHA",
                "SAC_PERC_EXPON",
                "SAC_PERC_PFREE",
                "BASEFLOW_COEFF",
                "UNAVAIL_FRAC",
            ),
            "pl": [
                PL(
                    name="[DEFAULT]",
                    values=(1.0, P.X11, P.X10, P.X09, 0, 0),
                ),
                PL(
                    name="UZF",
                    values=(1.0, P.X11, P.X10, P.X09, P.X03, 0),
                ),
                PL(
                    name="LZFP",
                    values=(1.0, P.X11, P.X10, P.X09, P.X01, P.X16),
                ),
                PL(
                    name="LZFS",
                    values=(1.0, P.X11, P.X10, P.X09, P.X02, 0),
                ),
            ],
        },
        alias="SoilParameterList",
    )
    land_use_parameter_list: LandUseParameterList = Field(
        {
            "parameters": (
                "GAMMA_SHAPE",
                "GAMMA_SCALE",
                "MELT_FACTOR",
                "STREAM_FRACTION",
                "MAX_SAT_AREA_FRAC",
                "BF_LOSS_FRACTION",
            ),
            "pl": [
                PL(
                    name="[DEFAULT]",
                    values=(
                        P.X18,
                        P.X19,
                        P.X17,
                        P.X15,
                        P.X14,
                        P.X12,
                    ),
                )
            ],
        },
        alias="LandUseParameterList",
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        {
            "parameters": ("RAIN_ICEPT_PCT", "SNOW_ICEPT_PCT"),
            "pl": [PL(name="[DEFAULT]", values=(0, 0))],
        },
        alias="VegetationParameterList",
    )
    hru_state_variable_table: rc.HRUStateVariableTable = Field(
        [
            rc.HRUState(
                hru_id=1,
                data={
                    "SOIL[0]": P.X04 * 1000,
                    "SOIL[2]": P.X06 * 1000,
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
                gauge.rain_correction = self.params.X20
                gauge.snow_correction = self.params.X21
