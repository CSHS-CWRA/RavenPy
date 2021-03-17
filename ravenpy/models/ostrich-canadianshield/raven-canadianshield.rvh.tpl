#########################################################################
:FileType          rvh ASCII Raven 3.0.4
:WrittenBy         Robert Chlumsky, James Craig & Juliane Mai
:CreationDate      Feb 2021
#
# Emulation of Canadian Shield simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#
#
:SubBasins
    :Attributes              NAME    DOWNSTREAM_ID   PROFILE   REACH_LENGTH    GAUGED
    :Units                   none    none            none      km              none
    1,             canadianshield,   -1,             NONE,     _AUTO,          1
:EndSubBasins


# tied parameters:
#    para_area_organic_SB1 = par_area_organic_SB1 = AREA_SB1 *    para_x34  = 4250.6 * par_x34
#    para_area_bedrock_SB1 = par_area_bedrock_SB1 = AREA_SB1 * (1-para_x34) = 4250.6 * (1-par_x34)


# split area based on parameter par_area_soilorg and par_area_soilbedrock
:HRUs
   :Attributes                 AREA    ELEVATION  LATITUDE    LONGITUDE  BASIN_ID  LAND_USE_CLASS  VEG_CLASS  SOIL_PROFILE AQUIFER_PROFILE TERRAIN_CLASS    SLOPE   ASPECT
   :Units                       km2            m       deg          deg      none            none       none          none            none          none      deg      deg
            1  par_area_organic_SB1,       843.0,  54.4848,   -123.3659,        1          FOREST     FOREST     SOILP_ORG          [NONE]        [NONE]      0.0      0.0
            2  par_area_bedrock_SB1,       843.0,  54.4848,   -123.3659,        1          FOREST     FOREST SOILP_BEDROCK          [NONE]        [NONE]      0.0      0.0
   #        1 para_area_organic_SB1,       843.0,  54.4848,   -123.3659,        1          FOREST     FOREST     SOILP_ORG          [NONE]        [NONE]      0.0      0.0
   #        2 para_area_bedrock_SB1,       843.0,  54.4848,   -123.3659,        1          FOREST     FOREST SOILP_BEDROCK          [NONE]        [NONE]      0.0      0.0
:EndHRUs
