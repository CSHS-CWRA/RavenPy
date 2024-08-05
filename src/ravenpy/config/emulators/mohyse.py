from collections.abc import Sequence
from typing import Literal, Union

from pydantic import Field, field_validator
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

import ravenpy.config.processes as p
from ravenpy.config import commands as rc
from ravenpy.config import options as o
from ravenpy.config.base import Params, Sym, SymConfig
from ravenpy.config.commands import HRU, PL, Process
from ravenpy.config.defaults import nc_attrs
from ravenpy.config.rvs import Config


@dataclass(config=SymConfig)
class P(Params):
    X01: Sym = Variable("X01")  #
    X02: Sym = Variable("X02")
    X03: Sym = Variable("X03")
    X04: Sym = Variable("X04")
    X05: Sym = Variable("X05")
    X06: Sym = Variable("X06")
    X07: Sym = Variable("X07")
    X08: Sym = Variable("X08")
    X09: Sym = Variable("X09")
    X10: Sym = Variable("X10")


class LandHRU(HRU):
    land_use_class: str = "LU_ALL"
    veg_class: str = "VEG_ALL"
    soil_profile: str = "DEFAULT_P"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    hru_type: Literal["land"] = "land"


class HRUs(rc.HRUs):
    """HRUs command for GR4J.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    root: Sequence[LandHRU]


class Mohyse(Config):
    """Modèle Hydrologique Simplifié à l'Extrême (MOHYSE)

    References
    ----------
    Fortin, V.; Turcotte, R. Le modèle hydrologique MOHYSE. In Note de Cours Pour SCA7420,
    Université du Québec à Montréal: Montréal, QC, Canada, 2007; p. 14.

    Troin, M., Arsenault, R. and Brissette, F., 2015. Performance and uncertainty evaluation
    of snow models on snowmelt flow simulations over a Nordic catchment (Mistassibi, Canada).
    Hydrology, 2(4), pp.289-317.
    """

    params: P = P()
    hrus: HRUs = Field([LandHRU()], alias="HRUs")
    netcdf_attribute: dict[str, str] = {"model_id": "Mohyse"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    soil_classes: rc.SoilClasses = Field(
        [{"name": "TOPSOIL"}, {"name": "GWSOIL"}], alias="SoilClasses"
    )
    land_use_classes: rc.LandUseClasses = Field(
        [{"name": "LU_ALL", "impermeable_frac": 0, "forest_coverage": 1}],
        alias="LandUseClasses",
    )
    vegetation_classes: rc.VegetationClasses = Field(
        [{"name": "VEG_ALL"}], alias="VegetationClasses"
    )
    soil_profiles: rc.SoilProfiles = Field(
        [
            {"name": "LAKE"},
            {"name": "ROCK"},
            {
                "name": "DEFAULT_P",
                "soil_classes": ["TOPSOIL", "GWSOIL"],
                "thicknesses": [P.X05, 10.0],
            },
        ],
        alias="SoilProfiles",
    )
    global_parameter: dict = Field(
        {"RAINSNOW_TEMP": -2, "TOC_MULTIPLIER": 1, "MOHYSE_PET_COEFF": P.X01},
        alias="GlobalParameter",
    )
    soil_parameter_list: rc.SoilParameterList = Field(
        {
            "parameters": [
                "POROSITY",
                "PET_CORRECTION",
                "HBV_BETA",
                "BASEFLOW_COEFF",
                "PERC_COEFF",
            ],
            "pl": [
                PL(name="TOPSOIL", values=(1.0, 1.0, 1.0, P.X07, P.X06)),
                PL(name="GWSOIL", values=(1.0, 1.0, 1.0, P.X08, 0.0)),
            ],
        },
        alias="SoilParameterList",
    )
    land_use_parameter_list: rc.LandUseParameterList = Field(
        {
            "parameters": [
                "MELT_FACTOR",
                "AET_COEFF",
                "FOREST_SPARSENESS",
                "DD_MELT_TEMP",
            ],
            "pl": [PL(name="[DEFAULT]", values=(P.X03, P.X02, 0, P.X04))],
        },
        alias="LandUseParameterList",
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        {
            "parameters": ["SAI_HT_RATIO", "RAIN_ICEPT_PCT", "SNOW_ICEPT_PCT"],
            "pl": [PL(name="[DEFAULT]", values=(0, 0, 0))],
        },
        alias="VegetationParameterList",
    )
    soil_model: rc.SoilModel = Field(2, alias="SoilModel")
    potential_melt_method: o.PotentialMeltMethod = Field(
        "POTMELT_DEGREE_DAY", alias="PotentialMeltMethod"
    )
    routing: o.Routing = Field("ROUTE_NONE", alias="Routing")
    catchment_route: o.CatchmentRoute = Field(
        "ROUTE_GAMMA_CONVOLUTION", alias="CatchmentRoute"
    )
    evaporation: o.Evaporation = Field("PET_MOHYSE", alias="Evaporation")
    direct_evaporation: bool = Field(True, alias="DirectEvaporation")
    rain_snow_fraction: o.RainSnowFraction = Field(
        "RAINSNOW_DATA", alias="RainSnowFraction"
    )
    hydrologic_processes: Sequence[Process] = Field(
        [
            p.SoilEvaporation(
                algo="SOILEVAP_LINEAR", source="SOIL[0]", to="ATMOSPHERE"
            ),
            p.SnowBalance(algo="SNOBAL_SIMPLE_MELT", source="SNOW", to="PONDED_WATER"),
            p.Precipitation(algo="RAVEN_DEFAULT", source="ATMOS_PRECIP", to="MULTIPLE"),
            p.Infiltration(algo="INF_HBV", source="PONDED_WATER", to="MULTIPLE"),
            p.Baseflow(algo="BASE_LINEAR", source="SOIL[0]", to="SURFACE_WATER"),
            p.Percolation(algo="PERC_LINEAR", source="SOIL[0]", to="SOIL[1]"),
            p.Baseflow(algo="BASE_LINEAR", source="SOIL[1]", to="SURFACE_WATER"),
        ]
    )
    sub_basin_properties: rc.SubBasinProperties = Field(
        {
            "parameters": ["GAMMA_SCALE", "GAMMA_SHAPE"],
            "records": [
                {"sb_id": "1", "values": (1 / P.X10, P.X09)},
            ],
        },
        alias="SubBasinProperties",
    )
    _nc_attrs = field_validator("netcdf_attribute")(nc_attrs)
