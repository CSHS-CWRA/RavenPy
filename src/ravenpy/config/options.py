from enum import Enum
from typing import Literal

Forcings = Literal[
    "PRECIP",
    "PRECIP_DAILY_AVE",
    "PRECIP_5DAY",
    "SNOW_FRAC",
    "SNOWFALL",
    "RAINFALL",
    "RECHARGE",
    "TEMP_AVE",
    "TEMP_DAILY_AVE",
    "TEMP_MIN",
    "TEMP_DAILY_MIN",
    "TEMP_MAX",
    "TEMP_DAILY_MAX",
    "TEMP_MONTH_MAX",
    "TEMP_MONTH_MIN",
    "TEMP_MONTH_AVE",
    "TEMP_AVE_UNC",
    "TEMP_MAX_UNC",
    "TEMP_MIN_UNC",
    "AIR_DENS",
    "AIR_PRES",
    "REL_HUMIDITY",
    "ET_RADIA",
    "SHORTWAVE",
    "SW_RADIA",
    "SW_RADIA_NET",
    "LW_RADIA_NET",
    "LW_INCOMING",
    "CLOUD_COVER",
    "DAY_LENGTH",
    "DAY_ANGLE",
    "WIND_VEL",
    "PET",
    "OW_PET",
    "PET_MONTH_AVE",
    "POTENTIAL_MELT",
    "SUBDAILY_CORR",
]

# Allowed soil parameters
SoilParameters = Literal[
    "SAND_CON",
    "CLAY_CON",
    "SILT_CON",
    "ORG_CON",
    "POROSITY",
    "STONE_FRAC",
    "SAT_WILT",
    "FIELD_CAPACITY",
    "BULK_DENSITY",
    "HYDRAUL_COND",
    "CLAPP_B",
    "CLAPP N,CLAPP M",
    "SAT_RES",
    "AIR_ENTRY_PRESSURE",
    "WILTING_PRESSURE",
    "HEAT_CAPACITY",
    "THERMAL_COND",
    "WETTING_FRONT_PSI",
    "EVAP_RES_FC",
    "SHUTTLEWORTH_B",
    "ALBEDO_WET",
    "ALBEDO_DRY",
    "VIZ_ZMIN",
    "VIC_ZMAX",
    "VIC ALPHA",
    "VIC_EVAP_GAMMA",
    "MAX_PERC_RATE",
    "PERC_N",
    "PERC_COEFF",
    "SAC_PERC_ALPHA",
    "SAC_PERC_EXPON",
    "SAC_PERC_PFREE",
    "UNAVAIL_FRAC",
    "HBV_BETA",
    "MAX_BASEFLOW_RATE",
    "BASEFLOW_N",
    "BASEFLOW_COEFF",
    "BASEFLOW_COEFF2",
    "BASEFLOW_THRESH",
    "BF_LOSS_FRACTION",
    "STORAGE_THRESHOLD",
    "MAX_CAP_RISE_RATE",
    "MAX_INTERFLOW_RATE",
    "INTERFLOW_COEF",
    "UBC_EVAL_SOIL_DEF",
    "UBC_INFIL_SOIL_DEF",
    "GR4J_X2",
    "GR4J_X3",
    "B_EXP",
    "PET_CORRECTION",  # Not in Raven Docs Table A.3
]

StateVariables = Literal[
    "SURFACE_WATER",
    "ATMOSPHERE",
    "ATMOS_PRECIP",
    "PONDED_WATER",
    "SOIL",
    "SOIL[0]",
    "SOIL[1]",
    "SOIL[2]",
    "GROUNDWATER",
    "CANOPY",
    "CANOPY_SNOW",
    "TRUNK",
    "ROOT",
    "DEPRESSION",
    "WETLAND",
    "LAKE_STORAGE",
    "SNOW",
    "SNOW_LIQ",
    "GLACIER",
    "GLACIER_ICE",
    "CONVOLUTION",
    "CONV_STOR",
    "SURFACE_WATER_TEMP",
    "SNOW_TEMP",
    "COLD_CONTENT",
    "GLACIER_CC",
    "SOIL_TEMP",
    "CANOPY_TEMP",
    "SNOW_DEPTH",
    "PERMAFROST_DEPTH",
    "SNOW_COVER",
    "SNOW_AGE",
    "SNOW_ALBEDO",
    "CROP_HEAT_UNITS",
    "CUM_INFIL",
    "CUM_SNOWMELT",
    "CONSTITUENT",
    "CONSTITUENT_SRC",
    "CONSTITUENT_SW",
    "CONSTITUENT_SINK",
    "MULTIPLE",
]

LandUseParameters = Literal[
    "FOREST_COVERAGE",
    "IMPERMEABLE_FRAC",
    "ROUGHNESS",
    "FOREST_SPARSENESS",
    "DEP_MAX",
    "MAX_DEP_AREA_FRAC",
    "DD_MELT_TEMP",
    "MELT_FACTOR",
    "DD_REFREEZE_TEMP",
    "MIN_MELT_FACTOR",
    "MAX_MELT_FACTOR",
    "REFREEZE_FACTOR",
    "REFREEZE_EXP",
    "DD_AGGRADATION",
    "SNOW_PATCH_LIMIT",
    "HBV_MELT_FOR_CORR",
    "HBV_MELT_ASP_CORR",
    "GLAC_STORAGE_COEFF",
    "HBV_MELT_GLACIER_CORR",
    "HBV_GLACIER_KMIN",
    "HBV_GLACIER_AG",
    "CC_DECAY_COEFF",
    "SCS_CN",
    "SCS_IA_FRACTION",
    "PARTITION_COEFF",
    "MAX_SAT_AREA_FRAC",
    "B_EXP",
    "ABST_PERCENT",
    "DEP_MAX_FLOW",
    "DEP_N",
    "DEP_SEEP_K",
    "DEP_K",
    "DEP_THRESHOLD",
    "PDMROF_B",
    "PONDED_EXP",
    "OW_PET_CORR",
    "LAKE_PET_CORR",
    "LAKE_REL_COEFF",
    "FOREST_PET_CORR",
    "GAMMA_SCALE",
    "GAMMA_SHAPE",
    "GAMMA_SCALE2",
    "GAMMA_SHAPE2",
    "HMETS_RUNOFF_COEFF",
    "AET_COEFF",
    "GR4J_X4",
    "UBC_ICEPT_FACTOR",
    "STREAM_FRACTION",
    "BF_LOSS_FRACTION",
]

VegetationParameters = Literal[
    "MAX_HEIGHT",
    "MAX_LEAF_COND",
    "MAX_LAI",
    "SVF_EXTINCTION",
    "RAIN_ICEPT_PCT",
    "SNOW_ICEPT_PCT",
    "RAIN_ICEPT_FACT",
    "SNOW_ICEPT_FACT",
    "SAI_HT_RATIO",
    "TRUNK_FRACTION",
    "STEMFLOW_FRAC",
    "ALBEDO",
    "ALBEDO_WET",
    "MAX_CAPACITY",
    "MAX_SNOW_CAPACITY",
    "ROOT_EXTINCT",
    "MAX_ROOT_LENGTH",
    "MIN_RESISTIVITY",
    "XYLEM_FRAC",
    "ROOTRADIUS",
    "PSI_CRITICAL",
    "DRIP_PROPORTION",
    "MAX_INTERCEPT_RATE",
    "CHU_MATURITY",
    "VEG_DIAM",
    "VEG_MBETA",
    "VEG_DENS" "PET_VEG_CORR",
    "TFRAIN",
    "TFSNOW",
    "RELATIVE_HT",
    "RELATIVE_LAI",
    "CAP_LAI_RATIO",
    "SNOCAP_LAI_RATIO",
]


class AirPressureMethod(Enum):
    BASIC = "AIRPRESS_BASIC"  # Default
    CONST = "AIRPRESS_CONST"
    DATA = "AIRPRESS_DATA"
    UBC = "AIRPRESS_UBC"


class Calendar(Enum):
    PROLEPTIC_GREGORIAN = "PROLEPTIC_GREGORIAN"
    JULIAN = "JULIAN"
    GREGORIAN = "GREGORIAN"
    STANDARD = "STANDARD"
    NOLEAP = "NOLEAP"
    _365_DAY = "365_DAY"
    ALL_LEAP = "ALL_LEAP"
    _366_DAY = "366_DAY"


class CatchmentRoute(Enum):
    """:CatchmentRoute"""

    DUMP = "ROUTE_DUMP"
    GAMMA = "ROUTE_GAMMA_CONVOLUTION"
    TRI = "ROUTE_TRI_CONVOLUTION"
    TRI_CONVOLUTION = "TRI_CONVOLUTION"
    TRIANGULAR_UH = "TRIANGULAR_UH"
    RESERVOIR = "ROUTE_RESERVOIR_SERIES"
    EXP = "ROUTE_EXPONENTIAL"
    DELAYED_FIRST_ORDER = "ROUTE_DELAYED_FIRST_ORDER"
    EXPONENTIAL_UH = "EXPONENTIAL_UH"


class CloudCoverMethod(Enum):
    NONE = "CLOUDCOV_NONE"  # default
    DATA = "CLOUDCOV_DATA"  # gauge or gridded time series used
    UBC = "CLOUDCOV_UBC"


class EvaluationMetrics(Enum):
    NASH_SUTCLIFFE = "NASH_SUTCLIFFE"
    LOG_NASH = "LOG_NASH"
    RMSE = "RMSE"
    PCT_BIAS = "PCT_BIAS"
    ABSERR = "ABSERR"
    ABSMAX = "ABSMAX"
    PDIFF = "PDIFF"
    TMVOL = "TMVOL"
    RCOEFF = "RCOEFF"
    NSC = "NSC"
    KLING_GUPTA = "KLING_GUPTA"
    DIAG_SPEARMAN = "DIAG_SPEARMAN"


evaluation_metrics_multiplier = dict(
    NASH_SUTCLIFFE=1,
    LOG_NASH=1,
    RMSE=-1,
    PCT_BIAS="Not Supported",
    ABSERR=-1,
    ABSMAX=-1,
    PDIFF=-1,
    TMVOL=-1,
    RCOEFF=-1,
    NSC=1,
    KLING_GUPTA=1,
    DIAG_SPEARMAN=1,
)


class Evaporation(Enum):
    CONSTANT = "PET_CONSTANT"
    PENMAN_MONTEITH = "PET_PENMAN_MONTEITH"
    PENMAN_COMBINATION = "PET_PENMAN_COMBINATION"
    PRIESTLEY_TAYLOR = "PET_PRIESTLEY_TAYLOR"
    HARGREAVES = "PET_HARGREAVES"
    HARGREAVES_1985 = "PET_HARGREAVES_1985"
    FROMMONTHLY = "PET_FROMMONTHLY"
    DATA = "PET_DATA"
    HAMON_1961 = "PET_HAMON_1961"
    TURC_1961 = "PET_TURC_1961"
    MAKKINK_1957 = "PET_MAKKINK_1957"
    MONTHLY_FACTOR = "PET_MONTHLY_FACTOR"
    MOHYSE = "PET_MOHYSE"
    OUDIN = "PET_OUDIN"
    VAP_DEFICIT = "PET_VAPDEFICIT"


class LWIncomingMethod(Enum):
    DATA = "LW_INC_DATA"
    DEFAULT = "LW_INC_DEFAULT"
    SICART = "LW_INC_SICART"
    SKYVIEW = "LW_INC_SKYVIEW"
    DINGMAN = "LW_INC_DINGMAN"


class Interpolation(Enum):
    FROM_FILE = "INTERP_FROM_FILE"
    AVERAGE_ALL = "INTERP_AVERAGE_ALL"
    NEAREST_NEIGHBOR = "INTERP_NEAREST_NEIGHBOR"
    INVERSE_DISTANCE = "INTERP_INVERSE_DISTANCE"


class LWRadiationMethod(Enum):
    DATA = "LW_RAD_DATA"
    DEFAULT = "LW_RAD_DEFAULT"
    UBCWM = "LW_RAD_UBC"


class MonthlyInterpolationMethod(Enum):
    UNIFORM = "MONTHINT_UNIFORM"
    LINEAR_MID = "MONTHINT_LINEAR_MID"
    LINEAR_FOM = "MONTHINT_LINEAR_FOM"
    LINEAR_21 = "MONTHINT_LINEAR_21"


class OroPETCorrect(Enum):
    NONE = "OROCORR_NONE"
    SIMPLELAPSE = "OROCORR_SIMPLELAPSE"
    HBV = "OROCORR_HBV"


class OroPrecipCorrect(Enum):
    NONE = "OROCORR_NONE"
    UBC = "OROCORR_UBC"
    HBV = "OROCORR_HBV"
    SIMPLELAPSE = "OROCORR_SIMPLELAPSE"


class OroTempCorrect(Enum):
    NONE = "OROCORR_NONE"
    HBV = "OROCORR_HBV"
    UBC = "OROCORR_UBC"
    UBC2 = "OROCORR_UBC_2"
    SIMPLELAPSE = "OROCORR_SIMPLELAPSE"


class PotentialMeltMethod(Enum):
    """:PotentialMelt algorithms"""

    DEGREE_DAY = "POTMELT_DEGREE_DAY"
    NONE = "POTMELT_NONE"
    RESTRICTED = "POTMELT_RESTRICTED"
    DATA = "POTMELT_DATA"
    EB = "POTMELT_EB"
    USACE = "POTMELT_USACE"
    HMETS = "POTMELT_HMETS"
    HBV = "POTMELT_HBV"
    UBC = "POTMELT_UBC"


class Precipitation(Enum):
    DEFAULT = "PRECIP_RAVEN"


class PrecipIceptFract(Enum):
    USER = "PRECIP_ICEPT_USER"  # default
    LAI = "PRECIP_ICEPT_LAI"
    EXPLAI = "PRECIP_ICEPT_EXPLAI"
    NONE = "PRECIP_ICEPT_NONE"
    HEDSTROM = "PRECIP_ICEPT_HEDSTROM"


class RainSnowFraction(Enum):
    DATA = "RAINSNOW_DATA"
    DINGMAN = "RAINSNOW_DINGMAN"
    UBC = "RAINSNOW_UBC"
    HBV = "RAINSNOW_HBV"
    HARDER = "RAINSNOW_HARDER"
    HSPF = "RAINSNOW_HSPF"
    WANG = "RAINSNOW_WANG"
    SNTHERM89 = "RAINSNOW_SNTHERM89"


class RelativeHumidityMethod(Enum):
    CONSTANT = "RELHUM_CONSTANT"
    DATA = "RELHUM_DATA"
    MINDEWPT = "RELHUM_MINDEWPT"
    CORR = "RELHUM_CORR"
    WINDVEL = "WINDVEL_CORR"


class Routing(Enum):
    DIFFUSIVE_WAVE = "ROUTE_DIFFUSIVE_WAVE"
    HYDROLOGIC = "ROUTE_HYDROLOGIC"
    NONE = "ROUTE_NONE"
    STORAGE_COEFF = "ROUTE_STORAGE_COEFF"
    PLUG_FLOW = "ROUTE_PLUG_FLOW"
    MUSKINGUM = "MUSKINGUM"


class SoilModel(Enum):
    ONE_LAYER = "SOIL_ONE_LAYER"
    TWO_LAYER = "SOIL_TWO_LAYER"
    MULTILAYER = "SOIL_MULTILAYER"


SubBasinProperties = Literal[
    "TIME_TO_PEAK",
    "TIME_CONC",
    "TIME_LAG",
    "NUM_RESERVOIRS",
    "RES_CONSTANT",
    "GAMMA_SHAPE",
    "GAMMA_SCALE",
    "Q_REFERENCE",
    "MANNINGS_N",
    "SLOPE",
    "DIFFUSIVITY",
    "CELERITY",
    "RAIN_CORR",
    "SNOW_CORR",
]


class SubdailyMethod(Enum):
    NONE = "SUBDAILY_NONE"
    SIMPLE = "SUBDAILY_SIMPLE"
    UBC = "SUBDAILY_UBC"


class SWCanopyCorrect(Enum):
    NONE = "SW_CANOPY_CORR_NONE"  # Default
    STATIC = "SW_CANOPY_CORR_STATIC"
    DYNAMIC = "SW_CANOPY_CORR_DYNAMIC"
    UBC = "SW_CANOPY_CORR_UBC"


class SWCloudCorrect(Enum):
    NONE = "SW_CLOUD_CORR_NONE"  # Default
    DINGMAN = "SW_CLOUD_CORR_DINGMAN"
    UBC = "SW_CLOUD_CORR_UBCWM"
    ANNANDALE = "SW_CLOUD_CORR_ANNANDALE"


class SWRadiationMethod(Enum):
    DATA = "SW_RAD_DATA"
    DEFAULT = "SW_RAD_DEFAULT"
    UBCWM = "SW_RAD_UBCWM"


class WindspeedMethod(Enum):
    CONSTANT = "WINDVEL_CONSTANT"
    DATA = "WINDVEL_DATA"
    UBC = "WINDVEL_UBC"


class EnKFMode(Enum):
    SPINUP = "ENKF_SPINUP"
    CLOSED_LOOP = "ENKF_CLOSED_LOOP"
    FORECAST = "ENKF_FORECAST"
    OPEN_LOOP = "ENKF_OPEN_LOOP"
    OPEN_FORECAST = "ENKF_OPEN_FORECAST"
