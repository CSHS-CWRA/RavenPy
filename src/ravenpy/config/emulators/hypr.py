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
    land_use_class: str = "OPEN_1"
    veg_class: str = "FOREST"
    soil_profile: str = "DEFAULT_P"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"


class HRUs(rc.HRUs):
    """HRUs command for GR4J.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    root: Sequence[LandHRU]


class HYPR(Config):
    """
    The HYdrological model for Prairie Region (HYPR) is based on the
    conceptual Hydrologiska Byrans Vattenbalansavdelning (HBV)-light
    model. HBV is modified to work in the prairies by incorporating a
    conceptual lateral-flow component to represent the pothole storage
    complexities. HYPR can be used for prairie streamflow simulation.

    References
    ----------
    Mohamed I. Ahmed, Amin Elshorbagy, and Alain Pietroniro (2020).
    Toward Simple Modeling Practices in the Complex Canadian Prairie
    Watersheds. Journal of Hydrologic Engineering, 25(6), 04020024.
    doi: 10.1061/(ASCE)HE.1943-5584.0001922.
    """

    # TRANSLATION OF HYPR PARAMETERS
    # -------------------------------------------------------
    # HYPR --> RAVEN
    # ................
    # TT              --> DD_MELT_TEMP                                         x01
    # TT              --> RAINSNOW_TEMP                                        x02
    # CWH             --> SNOW_SWI                                             x03
    # LP              --> FIELD_CAPACITY                                       x04
    # K0              --> BASEFLOW_COEFF[2]                                    x05
    # K1              --> BASEFLOW_COEFF[1]                                    x06
    # B               --> PDMROF_B                                             x07
    # MAXPA           --> MAX_DEP_AREA_FRAC                                    x08
    # PWR             --> PONDED_EXP = 2/PWR                                   x09
    # UZL             --> Threshold storage = STORAGE_THRESHOLD               x10
    # FC              --> max storage capacity  = THICKNESS[0]*POROSITY[0]     x11
    # CRFR            --> REFREEZE_FACTOR                                      x12
    # MRF             --> HBV_MELT_FOR_CORR                                    x13
    # CMAX            --> DEP_MAX = CMAX/(B+1)                                 x14
    # MAXBAS          --> TIME_CONC                                            x15
    # BETA            --> HBV_BETA                                             x16
    # C0=KMIN         --> MELT_FACTOR                                          x17
    # C0=KMIN         --> MIN_MELT_FACTOR                                      x18
    # TCALT           --> AdiabaticLapseRate (from HBV)                        x19
    # PCALT           --> PRECIP_LAPSE       (from HBV)                        x20
    # SCF             --> :SnowCorrection                                      x21

    params: P = P()
    hrus: HRUs = Field(
        [
            LandHRU(),
        ],
        alias="HRUs",
    )
    netcdf_attribute: dict[str, str] = {"model_id": "HYPR"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    routing: o.Routing = Field("ROUTE_NONE", alias="Routing")
    catchment_route: o.CatchmentRoute = Field(
        "ROUTE_TRI_CONVOLUTION", alias="CatchmentRouting"
    )
    evaporation: o.Evaporation = Field(o.Evaporation.FROMMONTHLY, alias="Evaporation")
    ow_evaporation: o.Evaporation = Field(
        o.Evaporation.FROMMONTHLY, alias="OW_Evaporation"
    )
    rain_snow_fraction: o.RainSnowFraction = Field(
        o.RainSnowFraction.HBV, alias="RainSnowFraction"
    )
    sw_radiation_method: o.SWRadiationMethod = Field(
        o.SWRadiationMethod.DEFAULT, alias="SWRadiationMethod"
    )
    sw_cloud_correct: o.SWCloudCorrect = Field(
        o.SWCloudCorrect.NONE, alias="SWCloudCorrect"
    )
    sw_canopy_correct: o.SWCanopyCorrect = Field(
        o.SWCanopyCorrect.NONE, alias="SWCanopyCorrect"
    )
    lw_radiation_method: o.LWRadiationMethod = Field(
        o.LWRadiationMethod.DEFAULT, alias="LWRadiationMethod"
    )
    potential_melt_method: o.PotentialMeltMethod = Field(
        o.PotentialMeltMethod.HBV, alias="PotentialMeltMethod"
    )
    cloud_cover_method: o.CloudCoverMethod = Field(
        o.CloudCoverMethod.NONE, alias="CloudCoverMethod"
    )
    precip_icept_frac: o.PrecipIceptFract = Field(
        o.PrecipIceptFract.USER, alias="PrecipIceptFrac"
    )
    soil_model: rc.SoilModel = Field(3, alias="SoilModel")

    hydrologic_processes: Sequence[Union[rc.Process, p.Conditional]] = Field(
        [
            p.SnowRefreeze(algo="FREEZE_DEGREE_DAY", source="SNOW_LIQ", to="SNOW"),
            p.Precipitation(algo="PRECIP_RAVEN", source="ATMOS_PRECIP", to="MULTIPLE"),
            p.CanopyEvaporation(algo="CANEVP_ALL", source="CANOPY", to="ATMOSPHERE"),
            p.CanopySublimation(
                algo="CANEVP_ALL", source="CANOPY_SNOW", to="ATMOSPHERE"
            ),
            p.SnowBalance(algo="SNOBAL_SIMPLE_MELT", source="SNOW", to="PONDED_WATER"),
            p.Infiltration(algo="INF_HBV", source="PONDED_WATER", to="MULTIPLE"),
            p.Flush(algo="RAVEN_DEFAULT", source="SURFACE_WATER", to="PONDED_WATER"),
            p.Abstraction(algo="ABST_PDMROF", source="PONDED_WATER", to="DEPRESSION"),
            p.Flush(algo="RAVEN_DEFAULT", source="SURFACE_WATER", to="SOIL[1]"),
            p.SoilEvaporation(algo="SOILEVAP_HYPR", source="MULTIPLE", to="ATMOSPHERE"),
            p.Baseflow(algo="BASE_LINEAR", source="SOIL[1]", to="SURFACE_WATER"),
            p.Baseflow(algo="BASE_THRESH_STOR", source="SOIL[1]", to="SURFACE_WATER"),
        ],
        alias="HydrologicProcesses",
    )
    lake_storage: o.StateVariables = Field("SOIL[2]", alias="LakeStorage")

    soil_classes: SoilClasses = Field(
        [{"name": "TOPSOIL"}, {"name": "SLOW_RES"}, {"name": "FAST_RES"}],
        alias="SoilClasses",
    )
    vegetation_classes: VegetationClasses = Field(
        [{"name": "FOREST", "max_ht": 0, "max_lai": 0, "max_leaf_cond": 1e99}],
        alias="VegetationClasses",
    )
    land_use_classes: LandUseClasses = Field(
        [rc.LU(name="OPEN_1", impermeable_frac=0.0, forest_coverage=0.0)],
        alias="LandUseClasses",
    )

    soil_profiles: SoilProfiles = Field(
        [
            {
                "name": "DEFAULT_P",
                "soil_classes": ("TOPSOIL", "FAST_RES", "SLOW_RES"),
                "thicknesses": (P.X11, 1e99, 1e99),
            },
        ],
        alias="SoilProfiles",
    )

    global_parameter: dict = Field(
        {
            "RAINSNOW_TEMP": P.X02,
            "RAINSNOW_DELTA": 0,
            "SNOW_SWI": P.X03,
            "AdiabaticLapseRate": P.X19,
            "PrecipLapse": P.X20,
        }
    )
    soil_parameter_list: SoilParameterList = Field(
        {
            "parameters": (
                "POROSITY",
                "FIELD_CAPACITY",
                "SAT_WILT",
                "HBV_BETA",
                "MAX_CAP_RISE_RATE",
                "MAX_PERC_RATE",
                "BASEFLOW_COEFF",
                "BASEFLOW_N",
                "BASEFLOW_COEFF2",
                "STORAGE_THRESHOLD",
            ),
            "pl": [
                PL(
                    name="[DEFAULT]",
                    values=(1.0, P.X04, 0.0, P.X16, 0, 0, 0, 0, 0, 0),
                ),
                PL(
                    name="FAST_RES",
                    values=(
                        "_DEFAULT",
                        "_DEFAULT",
                        0.0,
                        "_DEFAULT",
                        "_DEFAULT",
                        0,
                        P.X06,
                        1.0,
                        P.X05,
                        P.X10,
                    ),
                ),
                PL(
                    name="SLOW_RES",
                    values=(
                        "_DEFAULT",
                        "_DEFAULT",
                        0.0,
                        "_DEFAULT",
                        "_DEFAULT",
                        "_DEFAULT",
                        0.01,
                        1.0,
                        0.05,
                        0,
                    ),
                ),
            ],
        },
        alias="SoilParameterList",
    )
    land_use_parameter_list: LandUseParameterList = Field(
        {
            "parameters": (
                "MELT_FACTOR",
                "MIN_MELT_FACTOR",
                "HBV_MELT_FOR_CORR",
                "REFREEZE_FACTOR",
                "HBV_MELT_ASP_CORR",
                "DD_MELT_TEMP",
                "FOREST_COVERAGE",
                "PONDED_EXP",
                "PDMROF_B",
                "DEP_MAX",
                "MAX_DEP_AREA_FRAC",
            ),
            "pl": [
                PL(
                    name="[DEFAULT]",
                    values=(
                        P.X17,
                        P.X18,
                        P.X13,
                        P.X12,
                        0,
                        P.X01,
                        0,
                        P.X09,
                        P.X07,
                        P.X14,
                        P.X08,
                    ),
                )
            ],
        },
        alias="LandUseParameterList",
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        {
            "parameters": (
                "SAI_HT_RATIO",
                "MAX_CAPACITY",
                "MAX_SNOW_CAPACITY",
                "TFRAIN",
                "TFSNOW",
            ),
            "pl": [PL(name="FOREST", values=(0, 10000, 10000, 1, 1))],
        },
        alias="VegetationParameterList",
    )
    sub_basin_properties: rc.SubBasinProperties = Field(
        {"parameters": ["TIME_CONC"], "records": [{"sb_id": 1, "values": (P.X15,)}]},
        alias="SubBasinProperties",
    )

    hru_state_variable_table: rc.HRUStateVariableTable = Field(
        [
            rc.HRUState(
                hru_id=1,
                data={
                    "SOIL[0]": 10,
                    "SOIL[1]": 10,
                    "DEPRESSION": 20,
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
                gauge.snow_correction = self.params.X21
