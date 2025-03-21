from collections.abc import Sequence
from typing import Literal, Union

from pydantic import Field, field_validator
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

import ravenpy.config.processes as p
from ravenpy.config import commands as rc
from ravenpy.config import options as o
from ravenpy.config.base import Params, Sym, SymConfig
from ravenpy.config.commands import (
    HRU,
    PL,
    LandUseClasses,
    LandUseParameterList,
    Process,
    RainSnowTransition,
    SoilClasses,
    SoilModel,
    SoilParameterList,
    SoilProfiles,
    VegetationClasses,
)
from ravenpy.config.defaults import nc_attrs
from ravenpy.config.rvs import Config


@dataclass(config=SymConfig)
class P(Params):
    GR4J_X1: Sym = Variable("GR4J_X1")
    GR4J_X2: Sym = Variable("GR4J_X2")
    GR4J_X3: Sym = Variable("GR4J_X3")
    GR4J_X4: Sym = Variable("GR4J_X4")
    CEMANEIGE_X1: Sym = Variable("CEMANEIGE_X1")
    CEMANEIGE_X2: Sym = Variable("CEMANEIGE_X2")


# Note that the `hru_type` field is not part of the Raven configuration specification, but used here to automatically
# detect whether an HRU is land or lake. Just pass `hru_type` within the dictionary of attributes used to instantiate
# HRUs.


class LandHRU(HRU):
    land_use_class: str = "LU_ALL"
    veg_class: str = "VEG_ALL"
    soil_profile: str = "DEFAULT_P"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    hru_type: Literal["land"] = "land"


class LakeHRU(HRU):
    land_use_class: str = "LU_WATER"
    veg_class: str = "VEG_WATER"
    soil_profile: str = "LAKE"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    hru_type: Literal["lake"] = "lake"


class HRUs(rc.HRUs):
    """
    HRUs command for GR4J.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    root: Sequence[Union[LandHRU, LakeHRU]]


class GR4JCN(Config):
    """
    GR4J + Cemaneige global hydrological model.

    References
    ----------
    Perrin, C., C. Michel and V. Andréassian (2003). Improvement of a parsimonious model for streamflow simulation.
    Journal of Hydrology, 279(1-4), 275-289. doi: 10.1016/S0022-1694(03)00225-7.

    Valéry, Audrey, Vazken Andréassian, and Charles Perrin. 2014. “’As Simple as Possible but Not Simpler’: What Is
    Useful in a Temperature-Based Snow-Accounting Routine? Part 2 - Sensitivity Analysis of the Cemaneige Snow
    Accounting Routine on 380 Catchments.” Journal of Hydrology, no. 517(0): 1176–87,
    doi: 10.1016/j.jhydrol.2014.04.058.
    """

    params: P = P()
    hrus: HRUs = Field([LandHRU()], alias="HRUs")
    netcdf_attribute: dict[str, str] = {"model_id": "GR4JCN"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    uniform_initial_conditions: Union[dict[str, Sym], None] = Field(
        {"SOIL[0]": P.GR4J_X1 * 1000 / 2, "SOIL[1]": 15},
        alias="UniformInitialConditions",
    )
    evaporation: o.Evaporation = Field("PET_OUDIN", alias="Evaporation")
    rain_snow_fraction: o.RainSnowFraction = Field(
        "RAINSNOW_DINGMAN", alias="RainSnowFraction"
    )
    soil_model: Union[int, SoilModel] = Field(4, alias="SoilModel")
    routing: o.Routing = Field("ROUTE_NONE", alias="Routing")
    catchment_route: o.CatchmentRoute = Field("ROUTE_DUMP", alias="CatchmentRoute")
    potential_melt: o.PotentialMeltMethod = Field(
        "POTMELT_DEGREE_DAY", alias="PotentialMeltMethod"
    )
    oro_temp_correct: o.OroTempCorrect = Field(
        o.OroTempCorrect.SIMPLELAPSE, alias="OroTempCorrect"
    )
    oro_precip_correct: o.OroPrecipCorrect = Field(
        o.OroPrecipCorrect.SIMPLELAPSE, alias="OroPrecipCorrect"
    )

    hydrologic_processes: Sequence[Union[Process, p.Conditional]] = Field(
        [
            p.Precipitation(algo="PRECIP_RAVEN", source="ATMOS_PRECIP", to="MULTIPLE"),
            p.SnowTempEvolve(algo="SNOTEMP_NEWTONS", source="SNOW_TEMP"),
            p.SnowBalance(algo="SNOBAL_CEMA_NIEGE", source="SNOW", to="PONDED_WATER"),
            p.OpenWaterEvaporation(
                algo="OPEN_WATER_EVAP", source="PONDED_WATER", to="ATMOSPHERE"
            ),
            p.Infiltration(algo="INF_GR4J", source="PONDED_WATER", to="MULTIPLE"),
            p.SoilEvaporation(algo="SOILEVAP_GR4J", source="SOIL[0]", to="ATMOSPHERE"),
            p.Percolation(algo="PERC_GR4J", source="SOIL[0]", to="SOIL[2]"),
            p.Flush(
                algo="RAVEN_DEFAULT",
                source="SURFACE_WATER",
                to="SOIL[2]",
            ),
            p.Split(
                algo="RAVEN_DEFAULT",
                source="SOIL[2]",
                to=("CONVOLUTION[0]", "CONVOLUTION[1]"),
                p=0.9,
            ),
            p.Convolve(algo="CONVOL_GR4J_1", source="CONVOLUTION[0]", to="SOIL[1]"),
            p.Convolve(algo="CONVOL_GR4J_2", source="CONVOLUTION[1]", to="SOIL[2]"),
            p.Percolation(algo="PERC_GR4JEXCH", source="SOIL[1]", to="SOIL[3]"),
            p.Percolation(algo="PERC_GR4JEXCH2", source="SOIL[2]", to="SOIL[3]"),
            p.Flush(algo="RAVEN_DEFAULT", source="SOIL[2]", to="SURFACE_WATER"),
            p.Baseflow(algo="BASE_GR4J", source="SOIL[1]", to="SURFACE_WATER"),
        ],
        alias="HydrologicProcesses",
    )

    rain_snow_transition: RainSnowTransition = Field(
        {"temp": 0, "delta": 1}, alias="RainSnowTransition"
    )
    global_parameter: dict[str, Sym] = Field(
        {
            "PRECIP_LAPSE": 0.0004,
            "ADIABATIC_LAPSE": 0.0065,
            "AVG_ANNUAL_SNOW": P.CEMANEIGE_X1,
            "AIRSNOW_COEFF": 1 - P.CEMANEIGE_X2,
        },
        alias="GlobalParameter",
    )

    soil_classes: SoilClasses = Field(
        [
            {"name": "SOIL_PROD"},
            {"name": "SOIL_ROUT"},
            {"name": "SOIL_TEMP"},
            {"name": "SOIL_GW"},
            {"name": "AQUIFER"},
        ],
        alias="SoilClasses",
    )

    soil_parameter_list: SoilParameterList = Field(
        {
            "parameters": ("POROSITY", "GR4J_X3", "GR4J_X2"),
            "pl": [PL(name="[DEFAULT]", values=(1, P.GR4J_X3, P.GR4J_X2))],
        },
        alias="SoilParameterList",
    )

    soil_profiles: SoilProfiles = Field(
        [
            dict(
                name="DEFAULT_P",
                soil_classes=["SOIL_PROD", "SOIL_ROUT", "SOIL_TEMP", "SOIL_GW"],
                thicknesses=[P.GR4J_X1, 0.3, 1, 1],
            ),
            dict(name="LAKE"),
        ],
        alias="SoilProfiles",
    )

    vegetation_classes: VegetationClasses = Field(
        [{"name": "VEG_ALL"}, {"name": "VEG_WATER"}],
        alias="VegetationClasses",
    )

    land_use_classes: LandUseClasses = Field(
        [{"name": "LU_ALL"}, {"name": "LU_WATER"}],
        alias="LandUseClasses",
    )

    land_use_parameter_list: LandUseParameterList = Field(
        {
            "parameters": ("GR4J_X4", "MELT_FACTOR"),
            "pl": [PL(name="[DEFAULT]", values=(P.GR4J_X4, 7.73))],
        },
        alias="LandUseParameterList",
    )
    _nc_attrs = field_validator("netcdf_attribute")(nc_attrs)
