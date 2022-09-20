from typing import Literal

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
    "SAC_PERC_ALPHA",
    "SAC_PERC_EXPON",
    "HBV_BETA",
    "MAX_BASEFLOW_RATE",
    "BASEFLOW_N",
    "BASEFLOW_COEFF",
    "BASEFLOW_COEF2",
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
    "HMETS_RUNOFF_COEFF",
    "AET_COEFF",
    "GR4J_X4",
    "UBC_ICEPT_FACTOR",
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
    "PET_VEG_CORR",
]
