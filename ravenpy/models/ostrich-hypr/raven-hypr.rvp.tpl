#########################################################################
:FileType          rvp ASCII Raven 3.0.4
:WrittenBy         Mahkameh Taheri, Juliane Mai & James Craig
:CreationDate      Feb 2021
#
# Emulation of the HYPR model for Salmon River near Prince George
#------------------------------------------------------------------------


# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "para_x05" and "para_x06" wouldn't be detectable)
#    para_pow_x5     = pow_x05     =  10^(para_x05)       = 10^par_x05
#    para_pow_x6     = pow_x06     =  10^(para_x06)       = 10^par_x06


:SoilClasses
 :Attributes,
 :Units,
   TOPSOIL,
   SLOW_RES,
   FAST_RES,
:EndSoilClasses

:VegetationClasses
 :Attributes,   MAX_HT, MAX_LAI, MAX_LEAF_COND
 :Units,        m,      none,    mm_per_s
   FOREST,      0.0,     0.0,    1e99
:EndVegetationClasses


:LandUseClasses
 :Attributes, IMPERM, FOREST_COV
 :Units     ,   frac,       frac
      OPEN_1,    0.0,        0.0
:EndLandUseClasses


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
# UZL             --> Threshhold storage = STORAGE_THRESHOLD               x10
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
#

#
# ETF hard coded as 0.5 in Raven (here =0.0687) line 318 evaporation.cpp


#---------------------------------------------------------


#:SoilProfiles
# name, layers, [soilClass, thickness] x layers
#
:SoilProfiles
   DEFAULT_P, 3, TOPSOIL, par_x11, FAST_RES,1e99, SLOW_RES, 1e99
:EndSoilProfiles



# --- Parameter Specification ----------------------------

:GlobalParameter RAINSNOW_TEMP   par_x02
#                                para_x02 = TT
:GlobalParameter RAINSNOW_DELTA  0.0
:GlobalParameter SNOW_SWI        par_x03 # para_x03 = CWH
:AdiabaticLapseRate              par_x19 # para_x19 = TCALT
:GlobalParameter PRECIP_LAPSE    par_x20 # para_x20 = PCALT



:SoilParameterList
  :Parameters, POROSITY,   FIELD_CAPACITY, SAT_WILT,         HBV_BETA, MAX_CAP_RISE_RATE, MAX_PERC_RATE, BASEFLOW_COEFF, BASEFLOW_N, BASEFLOW_COEFF2, STORAGE_THRESHOLD
  :Units     ,     none,             none,     none,             none,              mm/d,          mm/d,            1/d,       none,             1/d,                mm,
#   [DEFAULT],      1.0,         para_x04,      0.0,         para_x16,              0.0,            0.0,            0.0,        0.0,             0.0,               0.0,
    [DEFAULT],      1.0,          par_x04,      0.0,          par_x16,              0.0,            0.0,            0.0,        0.0,             0.0,               0.0,
#    FAST_RES, _DEFAULT,         _DEFAULT,      0.0,         _DEFAULT,         _DEFAULT,            0.0,        pow__x06,        1.0,       pow__x05,          para_x10,
     FAST_RES, _DEFAULT,         _DEFAULT,      0.0,         _DEFAULT,         _DEFAULT,            0.0,        pow_x06,        1.0,         pow_x05,           par_x10,
     SLOW_RES, _DEFAULT,         _DEFAULT,      0.0,         _DEFAULT,         _DEFAULT,       _DEFAULT,           0.01,        1.0,            0.05,                 0,
:EndSoilParameterList

:VegetationParameterList
  :Parameters, SAI_HT_RATIO, MAX_CAPACITY, MAX_SNOW_CAPACITY, TFRAIN, TFSNOW
  :Units     ,         none,           mm,                mm,   frac,   frac
       FOREST,          0.0,        10000,             10000,    1.0,    1.0
:EndVegetationParameterList


:LandUseParameterList
  :Parameters ,      MELT_FACTOR,  MIN_MELT_FACTOR,   HBV_MELT_FOR_CORR,  REFREEZE_FACTOR, HBV_MELT_ASP_CORR,     DD_MELT_TEMP, FOREST_COVERAGE
  :Units      ,           mm/d/K,           mm/d/K,                none,           mm/d/K,              none,             degC,             0.1
  #           ,            KminA,            KminB,                 MRF,             CRFR,                AM,               TT,               0
  #  [DEFAULT],         para_x17,         para_x18,            para_x13,         para_x12,              0.00,         para_x01,             0.0
     [DEFAULT],          par_x17,          par_x18,             par_x13,          par_x12,              0.00,          par_x01,             0.0
:EndLandUseParameterList

:LandUseParameterList
  :Parameters,       PONDED_EXP,         PDMROF_B,             DEP_MAX, MAX_DEP_AREA_FRAC,
  :Units     ,             none,             none,                  mm,              none,
  #          ,            2/PWR,                B, SMAX=CMAX*(1/(B+1)),             MAXPA,
  # [DEFAULT],         para_x09,         para_x07,            para_x14,          para_x08,
  [DEFAULT],            par_x09,          par_x07,             par_x14,           par_x08,
:EndLandUseParameterList
