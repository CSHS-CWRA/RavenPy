from collections.abc import Sequence
from typing import Literal, Union

from pydantic import Field, field_validator

import ravenpy.config.processes as p
from ravenpy.config import commands as rc
from ravenpy.config import options as o
from ravenpy.config.commands import HRU, Process, SoilModel
from ravenpy.config.defaults import nc_attrs
from ravenpy.config.rvs import Config


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
    soil_profile: str = "Soil_Lake_HRU"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    hru_type: Literal["lake"] = "lake"


class HRUs(rc.HRUs):
    """HRUs command for GR4J.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    root: Sequence[Union[LandHRU, LakeHRU]]


class BasicRoute(Config):
    """Raven configuration performing routing only."""

    hrus: HRUs = Field([LandHRU()], alias="HRUs")
    netcdf_attribute: dict[str, str] = {"model_id": "BasicRoute"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
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
    _nc_attrs = field_validator("netcdf_attribute")(nc_attrs)
