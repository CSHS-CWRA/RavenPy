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
    """Model parameters.

    Default parameters values are mainly taken from the Raven RavenParameters.dat.
    Min and max values are taken either from the Raven RavenParameters.dat or the manual (v3.7),
    when both sources differed or when values were missing (value of -9999) in the RavenParameters.dat file.
    """

    X01: Sym = Variable(
        "X01"
    )  # RAINSNOW_TEMP      | rain/snow halfway transition temperature [◦C]                  | default=-0.15, min=-3.0, max=3.0
    X02: Sym = Variable(
        "X02"
    )  # MELT_FACTOR        | maximum snow melt factor used in degree day models [mm/d/◦C]   | default=5.04, min=-9999, max=-9999
    X03: Sym = Variable(
        "X03"
    )  # REFREEZE_FACTOR    | maximum refreeze factor used in degree day models [mm/d/◦C]    | default=5.04, min=0.0, max=10.0
    X04: Sym = Variable(
        "X04"
    )  # SNOW_SWI           | water saturation fraction of snow [-]                          | default=0.05, min=0.04, max=0.07
    X05: Sym = Variable(
        "X05"
    )  # POROSITY           | effective porosity of the soil [-]                             | default=0.4, min=0.0, max=1.0
    X06: Sym = Variable(
        "X06"
    )  # FIELD_CAPACITY     | field capacity saturation of the soil [-]                      | default=0.0, min=0.0, max=1.0
    X07: Sym = Variable(
        "X07"
    )  # HBV_BETA           | HBV infiltration exponent [-]                                  | default=1.0, min=0.0, max=7.0
    X08: Sym = Variable(
        "X08"
    )  # MAX_PERC_RATE      | VIC/ARNO/GAWSER percolation rate [mm/d]                        | default=-9999, min=0.01, max=1000
    X09: Sym = Variable(
        "X09"
    )  # BASEFLOW_COEFF     | linear baseflow storage/routing coeff (FAST_RES) [1/d]         | default=0.1, min=0.0, max=9999
    X10: Sym = Variable(
        "X10"
    )  # BASEFLOW_COEFF     | linear baseflow storage/routing coeff (SLOW_RES) [1/d]         | default=0.1, min=0.0, max=9999
    X11: Sym = Variable(
        "X11"
    )  # TIME_CONC          | the time of concentration of the unit hydrograph [d]           | default=1.0, min=0.0, max=200.0
    X12: Sym = Variable(
        "X12"
    )  # PRECIP_LAPSE       | precipitation lapse rate for orographic correction [mm/d/km]   | default=-9999, min=0.0, max=100.0
    X13: Sym = Variable(
        "X13"
    )  # ADIABATIC_LAPSE    | adiabatic temperature lapse rate [◦C/km]                       | default=6.4, min=5.0, max=8.0
    X14: Sym = Variable(
        "X14"
    )  # SAT_WILT           | hydroscopic minimum saturation [-]                             | default=0.0, min=0.0, max=1.0
    X15: Sym = Variable(
        "X15"
    )  # BASEFLOW_N         | VIC/ARNO baseflow exponent [-]                                 | default=5.0, min=1.0, max=10.0
    X16: Sym = Variable(
        "X16"
    )  # MAX_CAP_RISE_RATE  | HBV max capillary rise rate [mm/d]                             | default=0.0, min=-9999, max=-9999
    X17: Sym = Variable(
        "X17"
    )  # TOPSOIL THICKNESS  | soil layer thickness (TOPSOIL) [m]                             | default=-9999, min=0.0, max=-9999
    X18: Sym = Variable(
        "X18"
    )  # HBV_MELT_FOR_CORR  | HBV snowmelt forest correction (MRF in HBV-EC) [-]             | default=1.0, min=0.1, max=1.0
    X19: Sym = Variable(
        "X19"
    )  # GLAC_STORAGE_COEFF | maximum linear storage coefficient for glacial melt [-]        | default=0.1, min=0.0, max=1.0
    X20: Sym = Variable(
        "X20"
    )  # RAIN_CORR          | rain correction factor for subbasin (multiplier) [-]           | default=1.0, min=0.5, max=2.0
    X21: Sym = Variable(
        "X21"
    )  # SNOW_CORR          | snow correction factor for subbasin (multiplier) [-]           | default=1.0, min=0.5, max=2.0


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


class HBVEC(Config):
    """
    Hydrologiska Byråns Vattenbalansavdelning – Environment Canada (HBV-EC).

    References
    ----------
    Lindström, G., et al., 1997. Development and test of the distributed HBV-96 hydrological model.
    Journal of Hydrology, 201 (1–4), 272–288. doi:10.1016/S0022-1694(97)00041-3

    Hamilton, A.S., Hutchinson, D.G., and Moore, R.D., 2000. Estimating winter streamflow using
    conceptual streamflow model. Journal of Cold Regions Engineering, 14 (4), 158–175.
    doi:10.1061/(ASCE)0887-381X(2000)14:4(158)

    Canadian Hydraulics Centre, 2010. Green kenue reference manual.
    Ottawa, Ontario: National Research Council.
    """

    params: P = P()
    hrus: HRUs = Field([LandHRU()], alias="HRUs")
    netcdf_attribute: dict[str, str] = {"model_id": "HBVEC"}
    sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
    write_netcdf_format: bool = Field(True, alias="WriteNetcdfFormat")
    time_step: Union[float, str] = Field(1.0, alias="TimeStep")
    calendar: o.Calendar = Field("PROLEPTIC_GREGORIAN", alias="Calendar")
    interpolation: o.Interpolation = Field(
        "INTERP_NEAREST_NEIGHBOR", alias="Interpolation"
    )
    routing: o.Routing = Field("ROUTE_NONE", alias="Routing")
    catchment_route: o.CatchmentRoute = Field("TRIANGULAR_UH", alias="CatchmentRoute")
    evaporation: o.Evaporation = Field("PET_FROMMONTHLY", alias="Evaporation")
    ow_evaporation: o.Evaporation = Field("PET_FROMMONTHLY", alias="OW_Evaporation")
    sw_radiation_method: o.SWRadiationMethod = Field(
        "SW_RAD_DEFAULT", alias="SWRadiationMethod"
    )
    sw_cloud_correct: o.SWCloudCorrect = Field(
        "SW_CLOUD_CORR_NONE", alias="SWCloudCorrect"
    )
    sw_canopy_correct: o.SWCanopyCorrect = Field(
        "SW_CANOPY_CORR_NONE", alias="SWCanopyCorrect"
    )
    lw_radiation_method: o.LWRadiationMethod = Field(
        "LW_RAD_DEFAULT", alias="LWRadiationMethod"
    )
    rain_snow_fraction: o.RainSnowFraction = Field(
        "RAINSNOW_HBV", alias="RainSnowFraction"
    )
    potential_melt_method: o.PotentialMeltMethod = Field(
        "POTMELT_HBV", alias="PotentialMeltMethod"
    )
    oro_temp_correct: o.OroTempCorrect = Field("OROCORR_HBV", alias="OroTempCorrect")
    oro_precip_correct: o.OroPrecipCorrect = Field(
        "OROCORR_HBV", alias="OroPrecipCorrect"
    )
    oro_pet_correct: o.OroPETCorrect = Field("OROCORR_HBV", alias="OroPETCorrect")
    cloud_cover_method: o.CloudCoverMethod = Field(
        "CLOUDCOV_NONE", alias="CloudCoverMethod"
    )
    precip_icept_frac: o.PrecipIceptFract = Field(
        "PRECIP_ICEPT_USER", alias="PrecipIcepFract"
    )
    monthly_interpolation_method: o.MonthlyInterpolationMethod = Field(
        "MONTHINT_LINEAR_21", alias="MonthlyInterpolationMethod"
    )
    soil_model: rc.SoilModel = Field(3, alias="SoilModel")
    lake_storage: o.StateVariables = Field("SOIL[2]", alias="LakeStorage")
    hydrologic_processes: Sequence[Union[Process, p.Conditional]] = Field(
        [
            p.SnowRefreeze(algo="FREEZE_DEGREE_DAY", source="SNOW_LIQ", to="SNOW"),
            p.Precipitation(algo="PRECIP_RAVEN", source="ATMOS_PRECIP", to="MULTIPLE"),
            p.CanopyEvaporation(algo="CANEVP_ALL", source="CANOPY", to="ATMOSPHERE"),
            p.CanopySublimation(
                algo="CANEVP_ALL", source="CANOPY_SNOW", to="ATMOSPHERE"
            ),
            p.SnowBalance(algo="SNOBAL_SIMPLE_MELT", source="SNOW", to="SNOW_LIQ"),
            p.Overflow(algo="RAVEN_DEFAULT", source="SNOW_LIQ", to="PONDED_WATER"),
            p.Flush(algo="RAVEN_DEFAULT", source="PONDED_WATER", to="GLACIER"),
            p.Conditional(kind="HRU_TYPE", op="IS", value="GLACIER"),
            p.GlacierMelt(algo="GMELT_HBV", source="GLACIER_ICE", to="GLACIER"),
            p.GlacierRelease(
                algo="GRELEASE_HBV_EC", source="GLACIER", to="SURFACE_WATER"
            ),
            p.Infiltration(algo="INF_HBV", source="PONDED_WATER", to="MULTIPLE"),
            p.Flush(algo="RAVEN_DEFAULT", source="SURFACE_WATER", to="SOIL[1]"),
            p.Conditional(kind="HRU_TYPE", op="IS_NOT", value="GLACIER"),
            p.SoilEvaporation(algo="SOILEVAP_HBV", source="SOIL[0]", to="ATMOSPHERE"),
            p.CapillaryRise(algo="CRISE_HBV", source="SOIL[1]", to="SOIL[0]"),
            p.LakeEvaporation(
                algo="LAKE_EVAP_BASIC", source="SOIL[2]", to="ATMOSPHERE"
            ),
            p.Percolation(algo="PERC_CONSTANT", source="SOIL[1]", to="SOIL[2]"),
            p.Baseflow(algo="BASE_POWER_LAW", source="SOIL[1]", to="SURFACE_WATER"),
            p.Baseflow(algo="BASE_LINEAR", source="SOIL[2]", to="SURFACE_WATER"),
        ]
    )

    sub_basin_properties: rc.SubBasinProperties = Field(
        {
            "parameters": ["TIME_CONC", "TIME_TO_PEAK"],
            "records": [{"sb_id": 1, "values": (P.X11, P.X11 / 2)}],
        },
        alias="SubBasinProperties",
    )

    rain_snow_transition: rc.RainSnowTransition = Field(
        {"temp": P.X01, "delta": 2}, alias="RainSnowTransition"
    )
    global_parameter: dict = Field(
        {
            "ADIABATIC_LAPSE": P.X13,
            "SNOW_SWI": P.X04,
            "PRECIP_LAPSE": P.X12,
        },
        alias="GlobalParameters",
    )
    soil_classes: rc.SoilClasses = Field(
        [
            {"name": "TOPSOIL", "mineral": (1, 0, 0), "organic": 0},
            {"name": "FAST_RES", "mineral": (1, 0, 0), "organic": 0},
            {"name": "SLOW_RES", "mineral": (1, 0, 0), "organic": 0},
        ],
        alias="SoilClasses",
    )
    soil_parameter_list: rc.SoilParameterList = Field(
        {
            "parameters": [
                "POROSITY",
                "FIELD_CAPACITY",
                "SAT_WILT",
                "HBV_BETA",
                "MAX_CAP_RISE_RATE",
                "MAX_PERC_RATE",
                "BASEFLOW_COEFF",
                "BASEFLOW_N",
            ],
            "pl": [
                PL(
                    name="[DEFAULT]",
                    values=(P.X05, P.X06, P.X14, P.X07, P.X16, 0.0, 0.1, 1.0),
                ),
                PL(
                    name="TOPSOIL",
                    values=(P.X05, P.X06, P.X14, P.X07, P.X16, 0.0, 0.1, 1.0),
                ),
                PL(
                    name="FAST_RES",
                    values=(
                        "_DEFAULT",
                        "_DEFAULT",
                        0.0,
                        "_DEFAULT",
                        "_DEFAULT",
                        P.X08,
                        P.X09,
                        P.X15 + 1,
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
                        P.X10,
                        1.0,
                    ),
                ),
            ],
        },
        alias="SoilParameterList",
    )
    soil_profiles: rc.SoilProfiles = Field(
        [
            {
                "name": "DEFAULT_P",
                "soil_classes": ["TOPSOIL", "FAST_RES", "SLOW_RES"],
                "thicknesses": (P.X17, 100, 100),
            }
        ],
        alias="SoilProfiles",
    )
    vegetation_classes: rc.VegetationClasses = Field(
        [{"name": "VEG_ALL", "max_ht": 25, "max_lai": 6, "max_leaf_cond": 5.3}],
        alias="VegetationClasses",
    )
    vegetation_parameter_list: rc.VegetationParameterList = Field(
        {
            "parameters": ["MAX_CAPACITY", "MAX_SNOW_CAPACITY", "TFRAIN", "TFSNOW"],
            "pl": [PL(name="VEG_ALL", values=(10000, 10000, 0.88, 0.88))],
        },
        alias="VegetationParameterList",
    )
    land_use_classes: rc.LandUseClasses = Field(
        [{"name": "LU_ALL", "impermeable_frac": 0, "forest_coverage": 1}],
        alias="LandUseClasses",
    )
    land_use_parameter_list: rc.LandUseParameterList = Field(
        {
            "parameters": [
                "MELT_FACTOR",
                "MIN_MELT_FACTOR",
                "HBV_MELT_FOR_CORR",
                "REFREEZE_FACTOR",
                "HBV_MELT_ASP_CORR",
                "HBV_MELT_GLACIER_CORR",
                "HBV_GLACIER_KMIN",
                "GLAC_STORAGE_COEFF",
                "HBV_GLACIER_AG",
            ],
            "pl": [
                PL(
                    name="[DEFAULT]",
                    values=(P.X02, 2.2, P.X18, P.X03, 0.48, 1.64, 0.05, P.X19, 0.05),
                )
            ],
        },
        alias="LandUseParameterList",
    )

    hru_state_variable_table: rc.HRUStateVariableTable = Field(
        [{"hru_id": 1, "data": {"SOIL[2]": 0.50657}}]
    )
    _nc_attrs = field_validator("netcdf_attribute")(nc_attrs)

    def __init__(self, **data):
        super().__init__(**data)

        if self.gauge:
            for gauge in self.gauge:
                gauge.rain_correction = self.params.X20
                gauge.snow_correction = self.params.X21
