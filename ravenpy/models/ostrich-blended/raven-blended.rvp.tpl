#########################################################################
:FileType          rvp ASCII Raven 3.0.4
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Feb 2021
#
# Emulation of Blended model simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x25" and "par_x14" and "par_x10" wouldn't be detectable)
#    para_sum_x24_x25 = sum_x24_x25 =  para_x24 + para_x25 = par_x24 + par_x25
#    para_sum_x13_x14 = sum_x13_x14 =  para_x13 + para_x14 = par_x13 + par_x14
#    para_sum_x09_x10 = sum_x09_x10 =  para_x09 + para_x10 = par_x09 + par_x10
#    para_pow_x04     = pow_x04     =  10^(para_x04)       = 10^par_x04
#    para_pow_x11     = pow_x11     =  10^(para_x11)       = 10^par_x11

# RVT should contain
#  :RainCorrection par_x33
#  :SnowCorrection par_x34

#-----------------------------------------------------------------
# Soil Classes
#-----------------------------------------------------------------
:SoilClasses
  :Attributes,
  :Units,
  TOPSOIL,
  PHREATIC,
  DEEP_GW
:EndSoilClasses

#-----------------------------------------------------------------
# Land Use Classes
#-----------------------------------------------------------------
{land_use_classes_cmd}

#-----------------------------------------------------------------
# Vegetation Classes
#-----------------------------------------------------------------
:VegetationClasses,
  :Attributes,        MAX_HT,       MAX_LAI, MAX_LEAF_COND,
       :Units,             m,          none,      mm_per_s,
       FOREST,             4,             5,             5,
:EndVegetationClasses

#-----------------------------------------------------------------
# Soil Profiles
#-----------------------------------------------------------------
:SoilProfiles
         LAKE, 0
         ROCK, 0
  DEFAULT_P, 3, TOPSOIL,    par_x29, PHREATIC,    par_x30, DEEP_GW, 1e6
# DEFAULT_P, 3, TOPSOIL,      x(29), PHREATIC,      x(30), DEEP_GW, 1e6
:EndSoilProfiles

#-----------------------------------------------------------------
# Terrain Classes
#-----------------------------------------------------------------
:TerrainClasses
  :Attributes,        hillslope_len, drainage_dens,            lambda,
       :Units,                   ??,            ??,                ??
    DEFAULT_T,                  1.0,           1.0,           par_x07
#                                                     TOPMODEL_LAMBDA x(7)
:EndTerrainClasses

#-----------------------------------------------------------------
# Global Parameters
#-----------------------------------------------------------------
:GlobalParameter         SNOW_SWI_MIN par_x13            # x(13)
:GlobalParameter         SNOW_SWI_MAX sum_x13_x14   # x(13)+x(14)
:GlobalParameter     SWI_REDUCT_COEFF par_x15            # x(15)
:GlobalParameter             SNOW_SWI par_x19            # x(19)
:GlobalParameter        RAINSNOW_TEMP par_x31            # x(31)
:GlobalParameter       RAINSNOW_DELTA par_x32            # x(32)
#:GlobalParameter      TOC_MULTIPLIER 1.0                   #


#-----------------------------------------------------------------
# Soil Parameters
#-----------------------------------------------------------------
:SoilParameterList
  :Parameters,        POROSITY,      PERC_COEFF,  PET_CORRECTION,  BASEFLOW_COEFF,      B_EXP,   HBV_BETA, MAX_BASEFLOW_RATE, BASEFLOW_N,      FIELD_CAPACITY,   SAT_WILT,
       :Units,               -,             1/d,               -,             1/d
      TOPSOIL,             1.0,         par_x28,         par_x08,         pow_x04,    par_x02,    par_x03,           par_x05,    par_x06,         sum_x09_x10,    par_x09,
     PHREATIC,             1.0,         par_x35,             0.0,         pow_x11,        0.0,        0.0,               0.0,    par_x12,                 0.0,        0.0,
      DEEP_GW,             1.0,             0.0,             0.0,             0.0,        0.0,        0.0,               0.0,        0.0,                 0.0,        0.0,
 #    TOPSOIL,             1.0,           x(28),           x(08),           x(04),      x(02),      x(03),             x(05),      x(06),         x(09)+x(10),      x(09),
 #   PHREATIC,             1.0,           x(35),             0.0,           x(11),        0.0,        0.0,               0.0,      x(12),                 0.0,        0.0,
 #    DEEP_GW,             1.0,             0.0,             0.0,             0.0,        0.0,        0.0,               0.0,        0.0,                 0.0,        0.0,
:EndSoilParameterList

#-----------------------------------------------------------------
# Land Use Parameters
#-----------------------------------------------------------------
:LandUseParameterList
  :Parameters, MIN_MELT_FACTOR,     MAX_MELT_FACTOR,    DD_MELT_TEMP,  DD_AGGRADATION, REFREEZE_FACTOR, REFREEZE_EXP, DD_REFREEZE_TEMP, HMETS_RUNOFF_COEFF,
       :Units,          mm/d/C,              mm/d/C,               C,            1/mm,          mm/d/C,            -,                C,                  -,
    [DEFAULT],         par_x24,         sum_x24_x25,         par_x26,         par_x27,         par_x18,      par_x17,          par_x16,            par_x01,
#                        x(24),         x(24)+x(25),           x(26),           x(27),           x(18),        x(17),            x(16),              x(01),
:EndLandUseParameterList
:LandUseParameterList
  :Parameters,   GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2,    FOREST_SPARSENESS,
       :Units,             -,               -,               -,               -,                    -,
    [DEFAULT],       par_x20,         par_x21,         par_x22,         par_x23,                  0.0,
    #                  x(20),           x(21),           x(22),           x(23),                  0.0,
:EndLandUseParameterList

#-----------------------------------------------------------------
# Vegetation Parameters
#-----------------------------------------------------------------
:VegetationParameterList
  :Parameters,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,    SAI_HT_RATIO
       :Units,               -,               -,               -
    [DEFAULT],             0.0,             0.0,             0.0
:EndVegetationParameterList

:SeasonalRelativeLAI
  FOREST, 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
:EndSeasonalRelativeLAI
:SeasonalRelativeHeight
  FOREST, 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
:EndSeasonalRelativeHeight
