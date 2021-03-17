#########################################################################
:FileType          rvp ASCII Raven 3.0.4
:WrittenBy         Robert Chlumsky, James Craig & Juliane Mai
:CreationDate      Feb 2021
#
# Emulation of Canadian Shield simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#

# special parameter for Canadian Shield model in lumped model
## distributes the total area between bedrock and organic soil dominated HRUs

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this some parameters wouldn't be detectable)
#    para_sum_x05_x06 = sum_x05_x06 = par_x05 + par_x06 = para_x05 + para_x06
#    para_sum_x16_x17 = sum_x16_x17 = par_x16 + par_x17 = para_x16 + para_x17
#    para_pow_x08     = par_pow_x08 = 10^(par_x08) = 10^para_x08
#    para_pow_x09     = par_pow_x09 = 10^(par_x09) = 10^para_x09


#-----------------------------------------------------------------
# Soil Classes
#-----------------------------------------------------------------
:SoilClasses
  :Attributes,
  :Units,
  TOPSOIL,
  VADOSE,
  FRACBEDROCK,
:EndSoilClasses

#-----------------------------------------------------------------
# Soil Profiles
#-----------------------------------------------------------------
:SoilProfiles
  LAKE, 0
  ROCK, 0
  SOILP_ORG,     3, TOPSOIL, par_x01,  VADOSE, par_x02,  FRACBEDROCK, par_x03,
  SOILP_BEDROCK, 3, TOPSOIL, par_x01,  VADOSE, 0.00,     FRACBEDROCK, par_x03,
# SOILP_ORG,     3, TOPSOIL, para_x01, VADOSE, para_x02, FRACBEDROCK, para_x03,
# SOILP_BEDROCK, 3, TOPSOIL, para_x01, VADOSE, 0.00,     FRACBEDROCK, para_x03,
:EndSoilProfiles

#-----------------------------------------------------------------
# Land Use Classes
#-----------------------------------------------------------------
# :LandUseClasses,
#   :Attributes,        IMPERM,    FOREST_COV,
#        :Units,          frac,          frac,
#   FOREST,                0.0,       0.02345 ## forest_frac,
# :EndLandUseClasses
{land_use_classes_cmd}

#-----------------------------------------------------------------
# Vegetation Classes
#-----------------------------------------------------------------
:VegetationClasses,
        :Attributes,            MAX_HT, MAX_LAI, MAX_LEAF_COND,
        :Units,                      m,    none,      mm_per_s,
        FOREST,                    5.0,     5.0,           5.0,
:EndVegetationClasses

#-----------------------------------------------------------------
# Global Parameters
#-----------------------------------------------------------------
:GlobalParameter          SNOW_SWI par_x15            # para_x15
:GlobalParameter      SNOW_SWI_MIN par_x16            # para_x16
:GlobalParameter      SNOW_SWI_MAX par_sum_x16_x17    # para_x16+para_x17
:GlobalParameter  SWI_REDUCT_COEFF par_x18            # para_x18
:GlobalParameter     RAINSNOW_TEMP par_x19            # para_x19
:GlobalParameter    RAINSNOW_DELTA par_x20            # para_x20
:GlobalParameter   MAX_SWE_SURFACE par_x21            # para_x21
:GlobalParameter    TOC_MULTIPLIER par_x22            # para_x22

#-----------------------------------------------------------------
# Soil Parameters
#-----------------------------------------------------------------
:SoilParameterList
  :Parameters,        POROSITY,  HBV_BETA,  BASEFLOW_COEFF,    BASEFLOW_N,  MAX_INTERFLOW_RATE,   FIELD_CAPACITY,   SAT_WILT,   MAX_PERC_RATE,  PET_CORRECTION,
       :Units,               -,         -,               -,             -,                   -,                -,          -,               -,               -,
      TOPSOIL,             1.0,   par_x07,             0.0,           0.0,             par_x12,  par_sum_x05_x06,    par_x05,         par_x13,         par_x04,
       VADOSE,             1.0,       0.0,     par_pow_x08,       par_x10,                 0.0,              0.0,        0.0,         par_x14,             0.0,
  FRACBEDROCK,             1.0,       0.0,     par_pow_x09,       par_x11,                 0.0,              0.0,        0.0,             0.0,             0.0,
#     TOPSOIL,             1.0,  para_x07,             0.0,           0.0,            para_x12, para_sum_x05_x06,   para_x05,        para_x13,        para_x04,
#      VADOSE,             1.0,       0.0,    para_pow_x08,      para_x10,                 0.0,              0.0,        0.0,        para_x14,             0.0,
# FRACBEDROCK,             1.0,       0.0,    para_pow_x09,      para_x11,                 0.0,              0.0,        0.0,             0.0,             0.0,
:EndSoilParameterList

# note: TOPSOIL FIELD_CAPACITY calculated as par_x05 + par_x06
# note: TOPSOIL BASEFLOW_COEFF calculated as  10^(par_x08)
# note: PHREATIC BASEFLOW_COEFF calculated as  10^(par_x09)

#-----------------------------------------------------------------
# Land Use Parameters
#-----------------------------------------------------------------
:LandUseParameterList
  :Parameters, FOREST_SPARSENESS,     MELT_FACTOR,    DD_MELT_TEMP, REFREEZE_FACTOR,         DEP_MAX,     OW_PET_CORR,
       :Units,                 -,               -,               -,               -,               -,               -,
       FOREST,               0.0,         par_x25,         par_x24,         par_x23,         par_x26,         par_x27,
#      FOREST,               0.0,        para_x25,        para_x24,        para_x23,        para_x26,        para_x27,
:EndLandUseParameterList

#-----------------------------------------------------------------
# Vegetation Parameters
#-----------------------------------------------------------------
:VegetationParameterList
  :Parameters,  SVF_EXTINCTION,   SAI_HT_RATIO, RAIN_ICEPT_FACT, SNOW_ICEPT_FACT,    MAX_CAPACITY, MAX_SNOW_CAPACITY,
       :Units,               -,              -,               -,               -,               -,                 -,
    FOREST,                0.5,            1.0,         par_x28,         par_x29,         par_x30,           par_x31,
#   FOREST,                0.5,            1.0,        para_x28,        para_x29,        para_x30,          para_x31,
:EndVegetationParameterList


# Can leave all as 1.0 to keep generic across any watershed. In Ontario would hardcode to reduce ratio over winter months.
:SeasonalRelativeLAI
        #                 J       F       M       A       M       J       J       A       S       O       N       D
        [DEFAULT]       1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0
:EndSeasonalRelativeLAI
:SeasonalRelativeHeight
        #                 J       F       M       A       M       J       J       A       S       O       N       D
        [DEFAULT]       1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0
:EndSeasonalRelativeHeight
