from typing import Dict, Literal, Sequence, Type, Union

from pydantic import Field
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

import ravenpy.new_config.processes as p
from ravenpy.new_config import commands as rc
from ravenpy.new_config import options as o
from ravenpy.new_config.base import Params, Sym, SymConfig
from ravenpy.new_config.commands import (
    HRU,
    PL,
    LandUseClasses,
    LandUseParameterList,
    Process,
    RainSnowTransition,
    SoilClasses,
    SoilModel,
    SoilParameterList,
    SoilProfile,
    VegetationClasses,
)
from ravenpy.new_config.rvs import Config


class LandHRU(HRU):
    land_use_class: str = "Landuse_Land_HRU"
    veg_class: str = "Veg_Land_HRU"
    soil_profile: str = "Soil_Land_HRU"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    hru_type: Literal["land"] = "land"


class LakeHRU(HRU):
    land_use_class: str = "Landuse_Lake_HRU"
    veg_class: str = "Veg_Lake_HRU"
    soil_profile: str = "Lake_Soil_Lake_HRU"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    hru_type: Literal["lake"] = "lake"


class HRUs(rc.Command):
    """HRUs command for GR4J.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    __root__: Sequence[Union[LandHRU, LakeHRU]]


class BasicRoute(Config):
    """Raven configuration performing routing only."""

    catchment_route: o.CatchmentRoute = Field(
        "ROUTE_DUMP",
        alias="CatchmentRoute",
        description="Catchment routing method, "
        "used to convey water from the catchment tributaries and rivulets to the subbasin outlets.",
    )
    routing: o.Routing = Field(
        "ROUTE_DIFFUSIVE_WAVE",
        alias="Routing",
        description="Channel routing method which is used to "
        "transport water from upstream to downstream within the main subbasin channels.",
    )
    precip_icept_frac: o.PrecipIceptFract = Field(
        "PRECIP_ICEPT_NONE",
        alias="PrecipIcepFract",
        description="Estimation of the precipitation interception fraction. In this routing "
        "model, stream input(s) are pretending to be precipitation going into "
        "Raven.",
    )
    potential_melt_method: o.PotentialMeltMethod = Field(
        "POTMELT_NONE",
        alias="PotentialMeltMethod",
        description="Estimation of the potential snow melt. In this routing model, snow melt processes are not "
        "relevant, thus using DEFAULT POTMELT_NONE method.",
    )
    soil_model: SoilModel = Field(
        1, alias="SoilModel", description="Single soil layer structure"
    )
    hydrologic_processes: Sequence[Process] = Field(
        [
            p.Precipitation(
                algo="PRECIP_RAVEN", source="ATMOS_PRECIP", to="PONDED_WATER"
            ),
            p.Flush(source="PONDED_WATER", to="SURFACE_WATER"),
        ],
        alias="HydrologicProcesses",
    )