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
    GAMMA_SHAPE: Sym = Variable("GAMMA_SHAPE")
    GAMMA_SCALE: Sym = Variable("GAMMA_SCALE")
    GAMMA_SHAPE2: Sym = Variable("GAMMA_SHAPE2")
    GAMMA_SCALE2: Sym = Variable("GAMMA_SCALE2")
    MIN_MELT_FACTOR: Sym = Variable("MIN_MELT_FACTOR")
    MAX_MELT_FACTOR: Sym = Variable("MAX_MELT_FACTOR")
    DD_MELT_TEMP: Sym = Variable("DD_MELT_TEMP")
    DD_AGGRADATION: Sym = Variable("DD_AGGRADATION")
    SNOW_SWI_MIN: Sym = Variable("SNOW_SWI_MIN")
    SNOW_SWI_MAX: Sym = Variable("SNOW_SWI_MAX")
    SWI_REDUCT_COEFF: Sym = Variable("SWI_REDUCT_COEFF")
    DD_REFREEZE_TEMP: Sym = Variable("DD_REFREEZE_TEMP")
    REFREEZE_FACTOR: Sym = Variable("REFREEZE_FACTOR")
    REFREEZE_EXP: Sym = Variable("REFREEZE_EXP")
    PET_CORRECTION: Sym = Variable("PET_CORRECTION")
    HMETS_RUNOFF_COEFF: Sym = Variable("HMETS_RUNOFF_COEFF")
    PERC_COEFF: Sym = Variable("PERC_COEFF")
    BASEFLOW_COEFF_1: Sym = Variable("BASEFLOW_COEFF_1")
    BASEFLOW_COEFF_2: Sym = Variable("BASEFLOW_COEFF_2")
    TOPSOIL: Sym = Variable("TOPSOIL")
    PHREATIC: Sym = Variable("PHREATIC")


class ForestHRU(HRU):
    land_use_class: str = "FOREST"
    veg_class: str = "FOREST"
    soil_profile: str = "DEFAULT_P"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    hru_type: Literal["land"] = "land"


class HRUs(rc.HRUs):
    """
    HRUs command for HMETS.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    root: Sequence[ForestHRU]


class HMETS(Config):
    """
    Hydrology Model - École de technologie supérieure (HMETS).

    References
    ----------
    Martel, J.-L., Demeester, K., Brissette, F., Arsenault, R., Poulin, A. 2017.
    HMETS: A simple and efficient hydrology model for teaching hydrological modelling,
    flow forecasting and climate change impacts. Int. J. Eng. Educ., 33, 1307–1316.
    """

    params: P = P()
    hrus: HRUs = Field([ForestHRU()], alias="HRUs")
    netcdf_attribute: dict[str, str] = {"model_id": "HMETS"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    uniform_initial_conditions: dict[str, Sym] = Field(
        {"SOIL[0]": P.TOPSOIL / 2, "SOIL[1]": P.PHREATIC / 2},
        alias="UniformInitialConditions",
    )
    potential_melt_method: o.PotentialMeltMethod = Field(
        "POTMELT_HMETS", alias="PotentialMeltMethod"
    )
    rain_snow_fraction: o.RainSnowFraction = Field(
        "RAINSNOW_DATA", alias="RainSnowFraction"
    )
    evaporation: o.Evaporation = Field("PET_OUDIN", alias="Evaporation")
    catchment_route: o.CatchmentRoute = Field("ROUTE_DUMP", alias="CatchmentRoute")
    routing: o.Routing = Field("ROUTE_NONE", alias="Routing")
    soil_model: rc.SoilModel = Field(2, alias="SoilModel")

    hydrologic_processes: Sequence[Process] = Field(
        [
            p.SnowBalance(algo="SNOBAL_HMETS", source="MULTIPLE", to="MULTIPLE"),
            p.Precipitation(algo="RAVEN_DEFAULT", source="ATMOS_PRECIP", to="MULTIPLE"),
            p.Infiltration(algo="INF_HMETS", source="PONDED_WATER", to="MULTIPLE"),
            p.Overflow(algo="OVERFLOW_RAVEN", source="SOIL[0]", to="CONVOLUTION[1]"),
            p.Baseflow(algo="BASE_LINEAR", source="SOIL[0]", to="SURFACE_WATER"),
            p.Percolation(algo="PERC_LINEAR", source="SOIL[0]", to="SOIL[1]"),
            p.Overflow(algo="OVERFLOW_RAVEN", source="SOIL[1]", to="CONVOLUTION[1]"),
            p.SoilEvaporation(algo="SOILEVAP_ALL", source="SOIL[0]", to="ATMOSPHERE"),
            p.Convolve(
                algo="CONVOL_GAMMA", source="CONVOLUTION[0]", to="SURFACE_WATER"
            ),
            p.Convolve(
                algo="CONVOL_GAMMA_2", source="CONVOLUTION[1]", to="SURFACE_WATER"
            ),
            p.Baseflow(algo="BASE_LINEAR", source="SOIL[1]", to="SURFACE_WATER"),
        ]
    )

    soil_classes: rc.SoilClasses = Field(
        [{"name": "TOPSOIL"}, {"name": "PHREATIC"}], alias="SoilClasses"
    )
    land_use_classes: rc.LandUseClasses = Field(
        [{"name": "FOREST", "forest_coverage": 1}], alias="LandUseClasses"
    )
    vegetation_classes: rc.VegetationClasses = Field(
        [{"name": "FOREST", "max_ht": 4, "max_lai": 5, "max_leaf_cond": 5}],
        alias="VegetationClasses",
    )
    soil_profiles: rc.SoilProfiles = Field(
        [
            {"name": "LAKE"},
            {"name": "ROCK"},
            {
                "name": "DEFAULT_P",
                "soil_classes": ("TOPSOIL", "PHREATIC"),
                "thicknesses": (P.TOPSOIL / 1000.0, P.PHREATIC / 1000.0),
            },
        ]
    )
    global_parameter: dict[str, Sym] = Field(
        {
            "SNOW_SWI_MIN": P.SNOW_SWI_MIN,
            "SNOW_SWI_MAX": P.SNOW_SWI_MAX,
            "SWI_REDUCT_COEFF": P.SWI_REDUCT_COEFF,
            "SNOW_SWI": 0.05,
        },
        alias="GlobalParameter",
    )
    soil_parameter_list: rc.SoilParameterList = Field(
        {
            "parameters": (
                "POROSITY",
                "PERC_COEFF",
                "PET_CORRECTION",
                "BASEFLOW_COEFF",
            ),
            "pl": [
                PL(
                    name="TOPSOIL",
                    values=(1, P.PERC_COEFF, P.PET_CORRECTION, P.BASEFLOW_COEFF_1),
                ),
                PL(name="PHREATIC", values=(1.0, 0, 0, P.BASEFLOW_COEFF_2)),
            ],
        },
        alias="SoilParameterList",
    )
    land_use_parameter_list: rc.LandUseParameterList = Field(
        {
            "parameters": [
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
            ],
            "pl": [
                PL(
                    name="[DEFAULT]",
                    values=(
                        P.MIN_MELT_FACTOR,
                        P.MAX_MELT_FACTOR,
                        P.DD_MELT_TEMP,
                        P.DD_AGGRADATION,
                        P.REFREEZE_FACTOR,
                        P.REFREEZE_EXP,
                        P.DD_REFREEZE_TEMP,
                        P.HMETS_RUNOFF_COEFF,
                        P.GAMMA_SHAPE,
                        P.GAMMA_SCALE,
                        P.GAMMA_SHAPE2,
                        P.GAMMA_SCALE2,
                    ),
                )
            ],
        },
        alias="LandUseParameterList",
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        {
            "parameters": ["RAIN_ICEPT_PCT", "SNOW_ICEPT_PCT"],
            "pl": [PL(name="[DEFAULT]", values=(0, 0))],
        },
        alias="VegetationParameterList",
    )
    _nc_attrs = field_validator("netcdf_attribute")(nc_attrs)
