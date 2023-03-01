from typing import Dict, Sequence, Type, Union

from pydantic import Field
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

import ravenpy.new_config.processes as p
from ravenpy.config import options as o
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


@dataclass(config=SymConfig)
class P(Params):
    GR4J_X1: Sym = Variable("GR4J_X1")
    GR4J_X2: Sym = Variable("GR4J_X2")
    GR4J_X3: Sym = Variable("GR4J_X3")
    GR4J_X4: Sym = Variable("GR4J_X4")
    CEMANEIGE_X1: Sym = Variable("CEMANEIGE_X1")
    CEMANEIGE_X2: Sym = Variable("CEMANEIGE_X2")


class GR4JCN(Config):
    params: P
    soil_model: Union[int, SoilModel] = Field(4, alias="SoilModel")
    catchment_route: o.CatchmentRoute = Field("ROUTE_DUMP", alias="CatchmentRoute")
    potential_melt: o.PotentialMeltMethod = Field(
        "POTMELT_DEGREE_DAY", alias="PotentialMeltMethod"
    )
    alias: Dict[str, str] = Field(
        {
            "PRODUCT_STORE": "SOIL[0]",
            "ROUTING_STORE": "SOIL[1]",
            "TEMP_STORE": "SOIL[2]",
            "GW_STORE": "SOIL[3]",
        },
        alias="Alias",
    )
    hydrologic_processes: Sequence[Process] = Field(
        [
            p.Precipitation(algo="PRECIP_RAVEN", source=["ATMOS_PRECIP", "MULTIPLE"]),
            p.SnowTempEvolve(algo="SNOTEMP_NEWTONS", source=["SNOW_TEMP"]),
            p.SnowBalance(algo="SNOWBAL_CEMA_NEIGE", source=["SNOW", "PONDED_WATER"]),
            p.OpenWaterEvaporation(
                algo="OPEN_WATER_EVAP", source=["PONDED_WATER", "ATMOSPHERE"]
            ),
            p.Infiltration(algo="INF_GR4J", source=["PONDED_WATER", "MULTIPLE"]),
            p.SoilEvaporation(
                algo="SOILEVAP_GR4J", source=["PRODUCT_STORE", "ATMOSPHERE"]
            ),
            p.Percolation(algo="PERC_GR4J", source=["PRODUCT_STORE", "TEMP_STORE"]),
            p.Flush(
                algo="RAVEN_DEFAULT",
                source=["RAVEN_DEFAULT", "SURFACE_WATER", "TEMP_STORE"],
            ),
            p.Split(
                algo="RAVEN_DEFAULT",
                source=["TEMP_STORE", "CONVOLUTION[0]", "CONVOLUTION[1]", "0.9"],
            ),
            p.Convolve(
                algo="CONVOL_GR4J_1", source=["CONVOLUTION[0]", "ROUTING_STORE"]
            ),
            p.Convolve(algo="CONVOL_GR4J_2", source=["CONVOLUTION[1]", "TEMP_STORE"]),
            p.Percolation(algo="PERC_GR4JEXCH", source=["ROUTING_STORE", "GW_STORE"]),
            p.Percolation(algo="PERC_GR4JEXCH2", source=["TEMP_STORE", "GW_STORE"]),
            p.Flush(algo="RAVEN_DEFAULT", source=["TEMP_STORE", "SURFACE_WATER"]),
            p.BaseFlow(algo="BASE_GR4J", source=["ROUTING_STORE", "SURFACE_WATER"]),
        ],
        alias="HydrologicProcesses",
    )

    air_snow_coeff: Sym = Field(1 - P.CEMANEIGE_X2, alias="AirSnowCoeff")
    avg_annual_snow: Sym = Field(P.CEMANEIGE_X1, alias="AvgAnnualSnow")
    rain_snow_transition: RainSnowTransition = Field(
        {"temp": 0, "delta": 1}, alias="RainSnowTransition"
    )
    global_parameter: Dict[str, str] = Field(
        {"PrecipitationLapseRate": 0.0004, "AdiabaticLapseRate": 0.0065},
        alias="GlobalParameter",
    )

    soil_classes: SoilClasses = Field(
        SoilClasses(
            names=["SOIL_PROD", "SOIL_ROUT", "SOIL_TEMP", "SOIL_GW", "AQUIFER"]
        ),
        alias="SoilClasses",
    )

    soil_parameter_list: SoilParameterList = Field(
        SoilParameterList(
            names=("POROSITY", "GR4J_X3", "GR4J_X2"),
            records=[PL(name="[DEFAULT]", vals=(1, P.GR4J_X3, P.GR4J_X2))],
        ),
        alias="SoilParameterList",
    )

    soil_profiles: Sequence[SoilProfile] = Field(
        [
            SoilProfile(
                name="DEFAULT_P",
                soil_classes=["SOIL_PROD", "SOIL_ROUT", "SOIL_TEMP", "SOIL_GW"],
                thicknesses=[P.GR4J_X1, 0.3, 1, 1],
            ),
            SoilProfile(name="LAKE"),
        ],
        alias="SoilProfiles",
    )

    vegetation_classes: VegetationClasses = Field(
        VegetationClasses.parse_obj([{"name": "VEG_ALL"}, {"name": "VEG_WATER"}]),
        alias="VegetationClasses",
    )

    land_use_classes: LandUseClasses = Field(
        LandUseClasses.parse_obj([{"name": "LU_ALL"}, {"name": "LU_WATER"}]),
        alias="LandUseClasses",
    )

    land_use_parameter_list: LandUseParameterList = Field(
        LandUseParameterList(
            names=("GR4J_X4", "MELT_FACTOR"),
            records=[PL(name="[DEFAULT]", vals=(P.GR4J_X4, 7.73))],
        ),
        alias="LandUseParameterList",
    )

    @dataclass
    class LandHRU(HRU):
        land_use_class: str = "LU_ALL"
        veg_class: str = "VEG_ALL"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        hru_type: str = "land"

    @dataclass
    class LakeHRU(HRU):
        land_use_class: str = "LU_WATER"
        veg_class: str = "VEG_WATER"
        soil_profile: str = "LAKE"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        hru_type: str = "lake"
