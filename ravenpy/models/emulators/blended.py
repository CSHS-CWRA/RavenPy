from collections import defaultdict
from pathlib import Path
from typing import cast

from pydantic.dataclasses import dataclass

from ravenpy.config.commands import HRU, LU, BasinIndexCommand, HRUState, Sub
from ravenpy.config.rvs import RVI
from ravenpy.models.base import Ostrich, Raven

from .gr4jcn import GR4JCN


class BLENDED(Raven):
    @dataclass
    class Params:
        par_x01: float
        par_x02: float
        par_x03: float
        par_x04: float
        par_x05: float
        par_x06: float
        par_x07: float
        par_x08: float
        par_x09: float
        par_x10: float
        par_x11: float
        par_x12: float
        par_x13: float
        par_x14: float
        par_x15: float
        par_x16: float
        par_x17: float
        par_x18: float
        par_x19: float
        par_x20: float
        par_x21: float
        par_x22: float
        par_x23: float
        par_x24: float
        par_x25: float
        par_x26: float
        par_x27: float
        par_x28: float
        par_x29: float
        par_x30: float
        par_x31: float
        par_x32: float
        par_x33: float
        par_x34: float
        par_x35: float
        par_r01: float
        par_r02: float
        par_r03: float
        par_r04: float
        par_r05: float
        par_r06: float
        par_r07: float
        par_r08: float

    @dataclass
    class ForestHRU(HRU):
        land_use_class: str = "FOREST"
        veg_class: str = "FOREST"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "DEFAULT_T"

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "blended")
        super().__init__(*args, **kwds)

        self.config.update(
            hrus=(BLENDED.ForestHRU(),),
            subbasins=(
                Sub(
                    subbasin_id=1,
                    name="sub_001",
                    downstream_id=-1,
                    profile="None",
                    gauged=True,
                ),
            ),
            land_use_classes=(
                LU("FOREST", impermeable_frac=0.0, forest_coverage=0.02345),
            ),
            evaporation="PET_OUDIN",
            rain_snow_fraction=RVI.RainSnowFractionOptions.HBV,
        )

        #########
        # R V P #
        #########

        rvp_tmpl = """
        # tied parameters:
        # (it is important for OSTRICH to find every parameter place holder somewhere in this file)
        # (without this "5.421821E-02" and "1.596675E-01" and "2.169724E-01" wouldn't be detectable)
        #    para_1.705806E+00 = 1.705806E+00 =  para_x24 + para_x25 = 1.651588E+00 + 5.421821E-02
        #    para_2.070342E-01 = 2.070342E-01 =  para_x13 + para_x14 = 4.736668E-02 + 1.596675E-01
        #    para_2.518350E-01 = 2.518350E-01 =  para_x09 + para_x10 = 3.486263E-02 + 2.169724E-01
        #    para_2.254976E-04     = 2.254976E-04     =  10^(para_x04)       = 10^-3.646858E+00
        #    para_2.279250E-04     = 2.279250E-04     =  10^(para_x11)       = 10^-3.642208E+00

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
        {land_use_classes}

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
                 DEFAULT_P, 3, TOPSOIL,  {params.par_x29}, PHREATIC, {params.par_x30}, DEEP_GW, 1e6
        #        DEFAULT_P, 3, TOPSOIL,         x(29), PHREATIC,        x(30), DEEP_GW, 1e6
        :EndSoilProfiles

        #-----------------------------------------------------------------
        # Terrain Classes
        #-----------------------------------------------------------------
        :TerrainClasses
          :Attributes,        hillslope_len, drainage_dens,                 lambda,
               :Units,                   ??,            ??,                     ??
            DEFAULT_T,                  1.0,           1.0,           {params.par_x07}
        #                                                     TOPMODEL_LAMBDA x(7)
        :EndTerrainClasses

        #-----------------------------------------------------------------
        # Global Parameters
        #-----------------------------------------------------------------
        :GlobalParameter         SNOW_SWI_MIN {params.par_x13}  #  x(13)
        :GlobalParameter         SNOW_SWI_MAX {params.par_x14}  #  x(13)+x(14)
        :GlobalParameter     SWI_REDUCT_COEFF {params.par_x15}  #  x(15)
        :GlobalParameter             SNOW_SWI {params.par_x19}  #  x(19)
        :GlobalParameter        RAINSNOW_TEMP {params.par_x31}  #  x(31)
        :GlobalParameter       RAINSNOW_DELTA {params.par_x32}  #  x(32)
        #:GlobalParameter      TOC_MULTIPLIER 1.0               #


        #-----------------------------------------------------------------
        # Soil Parameters
        #-----------------------------------------------------------------
        :SoilParameterList
          :Parameters,        POROSITY,        PERC_COEFF,    PET_CORRECTION,     BASEFLOW_COEFF,             B_EXP,          HBV_BETA,   MAX_BASEFLOW_RATE,       BASEFLOW_N,    FIELD_CAPACITY,         SAT_WILT,
               :Units,               -,               1/d,                 -,                1/d,
              TOPSOIL,             1.0,  {params.par_x28},  {params.par_x08},    {params.par_x04}, {params.par_x02},  {params.par_x03},    {params.par_x05}, {params.par_x06},  {params.par_x10}, {params.par_x09},
             PHREATIC,             1.0,  {params.par_x35},               0.0,    {params.par_x11},              0.0,               0.0,                 0.0, {params.par_x12},               0.0,              0.0,
              DEEP_GW,             1.0,               0.0,               0.0,                 0.0,              0.0,               0.0,                 0.0,              0.0,               0.0,              0.0,
         #    TOPSOIL,             1.0,             x(28),             x(08),               x(04),            x(02),             x(03),               x(05),            x(06),       x(09)+x(10),            x(09),
         #   PHREATIC,             1.0,             x(35),               0.0,               x(11),              0.0,               0.0,                 0.0,            x(12),               0.0,              0.0,
         #    DEEP_GW,             1.0,               0.0,               0.0,                 0.0,              0.0,               0.0,                 0.0,              0.0,               0.0,              0.0,
        :EndSoilParameterList

        #-----------------------------------------------------------------
        # Land Use Parameters
        #-----------------------------------------------------------------
        :LandUseParameterList
          :Parameters,  MIN_MELT_FACTOR,     MAX_MELT_FACTOR,     DD_MELT_TEMP,   DD_AGGRADATION,  REFREEZE_FACTOR,     REFREEZE_EXP,     DD_REFREEZE_TEMP,     HMETS_RUNOFF_COEFF,
               :Units,           mm/d/C,              mm/d/C,                C,             1/mm,           mm/d/C,                -,                    C,                      -,
            [DEFAULT], {params.par_x24},    {params.par_x25}, {params.par_x26}, {params.par_x27}, {params.par_x18}, {params.par_x17},     {params.par_x16},       {params.par_x01},
        #                        x(24),         x(24)+x(25),           x(26),           x(27),               x(18),            x(17),                x(16),                  x(01),
        :EndLandUseParameterList
        :LandUseParameterList
          :Parameters,       GAMMA_SHAPE,         GAMMA_SCALE,        GAMMA_SHAPE2,        GAMMA_SCALE2,    FOREST_SPARSENESS,
               :Units,                 -,                   -,                   -,                   -,                    -,
            [DEFAULT],  {params.par_x20},    {params.par_x21},    {params.par_x22},    {params.par_x23},                  0.0,
            #                      x(20),               x(21),               x(22),               x(23),                  0.0,
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
        """
        self.config.rvp.set_tmpl(rvp_tmpl)

        #########
        # R V I #
        #########

        rvi_tmpl = """
        :PotentialMeltMethod     POTMELT_HMETS
        :RainSnowFraction        {rain_snow_fraction}
        :Evaporation             {evaporation}         # PET_OUDIN
        :CatchmentRoute          ROUTE_DUMP
        :Routing                 ROUTE_NONE
        :SoilModel               SOIL_MULTILAYER 3

        :Alias DELAYED_RUNOFF CONVOLUTION[1]

        :HydrologicProcesses
          :Precipitation   RAVEN_DEFAULT                         ATMOS_PRECIP   MULTIPLE
          :ProcessGroup #infiltration group
                        :Infiltration    INF_HMETS               PONDED_WATER   MULTIPLE
                        :Infiltration    INF_VIC_ARNO            PONDED_WATER   MULTIPLE
                        :Infiltration    INF_HBV                 PONDED_WATER   MULTIPLE
          :EndProcessGroup CALCULATE_WTS {params.par_r01} {params.par_r02}
          #                              para_r01 para_r02
                          :Overflow      OVERFLOW_RAVEN          SOIL[0]        DELAYED_RUNOFF
          :ProcessGroup #quickflow group
                        :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[0]        SURFACE_WATER   # interflow, really
                        :Baseflow        BASE_VIC                SOIL[0]        SURFACE_WATER
                        :Baseflow        BASE_TOPMODEL           SOIL[0]        SURFACE_WATER
          :EndProcessGroup CALCULATE_WTS {params.par_r03} {params.par_r04}
          #                              para_r03 para_r04
          :Percolation                   PERC_LINEAR             SOIL[0]        SOIL[1]         # recharge
            :Overflow                    OVERFLOW_RAVEN          SOIL[1]        DELAYED_RUNOFF
          :Percolation                   PERC_LINEAR             SOIL[1]        SOIL[2]         # loss to deep groundwater (simplifies to HMETS when PERC_COEFF DEEP_GW=0)
          :ProcessGroup #evaporation group
                        :SoilEvaporation SOILEVAP_ALL            SOIL[0]        ATMOSPHERE      # AET
                        :SoilEvaporation SOILEVAP_TOPMODEL       SOIL[0]        ATMOSPHERE      # AET
          :EndProcessGroup CALCULATE_WTS {params.par_r05}
          #                              para_r05
          :Convolve                      CONVOL_GAMMA            CONVOLUTION[0] SURFACE_WATER   # 'surface runoff'
          :Convolve                      CONVOL_GAMMA_2          DELAYED_RUNOFF SURFACE_WATER   # 'delayed runoff'
          :ProcessGroup #quickflow group
                        :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[1]        SURFACE_WATER
                        :Baseflow        BASE_POWER_LAW          SOIL[1]        SURFACE_WATER
          :EndProcessGroup CALCULATE_WTS {params.par_r06}
          #                              para_r06
          :ProcessGroup #snow balance group
                        :SnowBalance     SNOBAL_HMETS            MULTIPLE       MULTIPLE
                        :SnowBalance     SNOBAL_SIMPLE_MELT      SNOW           PONDED_WATER
                        :SnowBalance     SNOBAL_HBV              MULTIPLE       MULTIPLE
                        #:SnowBalance    SNOBAL_GAWSER           MULTIPLE       MULTIPLE
          :EndProcessGroup CALCULATE_WTS {params.par_r07} {params.par_r08}
          #                              para_r07 para_r08
        :EndHydrologicProcesses
        """
        self.config.rvi.set_tmpl(rvi_tmpl)

    def derived_parameters(self):
        self.config.rvi.set_extra_attributes(params=self.config.rvp.params)

        params = cast(BLENDED.Params, self.config.rvp.params)

        topsoil_hlf = params.par_x29 * 0.5 * 1000
        phreatic_hlf = params.par_x30 * 0.5 * 1000
        hru_state = HRUState(soil0=topsoil_hlf, soil1=phreatic_hlf)
        self.config.rvc.set_hru_state(hru_state)

        self.config.rvt.rain_correction = params.par_x33
        self.config.rvt.snow_correction = params.par_x34


class BLENDED_OST(Ostrich, BLENDED):

    ostrich_to_raven_param_conversion = {
        "sum_x09_x10": "par_x10",
        "sum_x13_x14": "par_x14",
        "sum_x24_x25": "par_x25",
        "pow_x04": "par_x04",
        "pow_x11": "par_x11",
    }

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "blended-ost")
        super().__init__(*args, **kwds)

        self.config.update(
            algorithm="DDS",
            max_iterations=50,
            suppress_output=True,
        )

        ####################
        # R V C (OST TMPL) #
        ####################

        rvc_tmpl = """
        # tied parameters:
        # (it is important for OSTRICH to find every parameter place holder somewhere in this file)
        # (without this "par_x29" and "par_x30" wouldn't be detectable)
        #    para_half_x29 = para_x29 * 1000. / 2. = par_x29 / 2. [m] = half_x29 [mm]
        #    para_half_x30 = para_x30 * 1000. / 2. = par_x30 / 2. [m] = half_x30 [mm]

        # initialize to 1/2 full
        #:UniformInitialConditions SOIL[0] half_x29 # x(29)*1000/2 [mm]
        #:UniformInitialConditions SOIL[1] half_x30 # x(30)*1000/2 [mm]

        :HRUStateVariableTable (formerly :IntialConditionsTable)
           :Attributes SOIL[0] SOIL[1]
           :Units mm mm
           1 half_x29 half_x30
        :EndHRUStateVariableTable
        """
        self.config.rvc.set_tmpl(rvc_tmpl, is_ostrich=True)

        ####################
        # R V I (OST TMPL) #
        ####################

        rvi_tmpl = """
        :PotentialMeltMethod     POTMELT_HMETS
        :RainSnowFraction        {rain_snow_fraction}
        :Evaporation             {evaporation}         # PET_OUDIN
        :CatchmentRoute          ROUTE_DUMP
        :Routing                 ROUTE_NONE
        :SoilModel               SOIL_MULTILAYER 3

        :Alias DELAYED_RUNOFF CONVOLUTION[1]

        :HydrologicProcesses
          :Precipitation   RAVEN_DEFAULT                         ATMOS_PRECIP   MULTIPLE
          :ProcessGroup #infiltration group
                        :Infiltration    INF_HMETS               PONDED_WATER   MULTIPLE
                        :Infiltration    INF_VIC_ARNO            PONDED_WATER   MULTIPLE
                        :Infiltration    INF_HBV                 PONDED_WATER   MULTIPLE
          :EndProcessGroup CALCULATE_WTS par_r01 par_r02
                          :Overflow      OVERFLOW_RAVEN          SOIL[0]        DELAYED_RUNOFF
          :ProcessGroup #quickflow group
                        :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[0]        SURFACE_WATER   # interflow, really
                        :Baseflow        BASE_VIC                SOIL[0]        SURFACE_WATER
                        :Baseflow        BASE_TOPMODEL           SOIL[0]        SURFACE_WATER
          :EndProcessGroup CALCULATE_WTS par_r03 par_r04
          :Percolation                   PERC_LINEAR             SOIL[0]        SOIL[1]         # recharge
            :Overflow                    OVERFLOW_RAVEN          SOIL[1]        DELAYED_RUNOFF
          :Percolation                   PERC_LINEAR             SOIL[1]        SOIL[2]         # loss to deep groundwater (simplifies to HMETS when PERC_COEFF DEEP_GW=0)
          :ProcessGroup #evaporation group
                        :SoilEvaporation SOILEVAP_ALL            SOIL[0]        ATMOSPHERE      # AET
                        :SoilEvaporation SOILEVAP_TOPMODEL       SOIL[0]        ATMOSPHERE      # AET
          :EndProcessGroup CALCULATE_WTS par_r05
          :Convolve                      CONVOL_GAMMA            CONVOLUTION[0] SURFACE_WATER   # 'surface runoff'
          :Convolve                      CONVOL_GAMMA_2          DELAYED_RUNOFF SURFACE_WATER   # 'delayed runoff'
          :ProcessGroup #quickflow group
                        :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[1]        SURFACE_WATER
                        :Baseflow        BASE_POWER_LAW          SOIL[1]        SURFACE_WATER
          :EndProcessGroup CALCULATE_WTS par_r06
          :ProcessGroup #snow balance group
                        :SnowBalance     SNOBAL_HMETS            MULTIPLE       MULTIPLE
                        :SnowBalance     SNOBAL_SIMPLE_MELT      SNOW           PONDED_WATER
                        :SnowBalance     SNOBAL_HBV              MULTIPLE       MULTIPLE
                        #:SnowBalance    SNOBAL_GAWSER           MULTIPLE       MULTIPLE
          :EndProcessGroup CALCULATE_WTS par_r07 par_r08
        :EndHydrologicProcesses
        """
        self.config.rvi.set_tmpl(rvi_tmpl, is_ostrich=True)

        ####################
        # R V P (OST TMPL) #
        ####################

        rvp_tmpl = """
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
        {land_use_classes}

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
        """
        self.config.rvp.set_tmpl(rvp_tmpl, is_ostrich=True)

        ####################
        # R V T (OST TMPL) #
        ####################

        rvt_tmpl = """
        {gauge}

        {forcing_list}

        {observed_data}
        """
        self.config.rvt.set_tmpl(rvt_tmpl, is_ostrich=True)

        ##########
        # O S T  #
        ##########

        ost_tmpl = """
        ProgramType         {algorithm}
        ObjectiveFunction   GCOP
        ModelExecutable     ./ostrich-runs-raven.sh
        PreserveBestModel   ./save_best.sh
        #OstrichWarmStart   yes

        ModelSubdir processor_

        BeginExtraDirs
        model
        #best
        EndExtraDirs

        BeginFilePairs
          {identifier}.rvi.tpl;  {identifier}.rvi
          {identifier}.rvp.tpl;  {identifier}.rvp
          {identifier}.rvc.tpl;  {identifier}.rvc
          {identifier}.rvt.tpl;  {identifier}.rvt
          #can be multiple (.rvh, .rvi)
        EndFilePairs

        #Parameter/DV Specification
        BeginParams
          #parameter init.    low               high               tx_in  tx_ost   tx_out
          par_x01    random   {lowerBounds.par_x01} {upperBounds.par_x01}  none   none     none
          par_x02    random   {lowerBounds.par_x02} {upperBounds.par_x02}  none   none     none
          par_x03    random   {lowerBounds.par_x03} {upperBounds.par_x03}  none   none     none
          par_x04    random   {lowerBounds.par_x04} {upperBounds.par_x04}  none   none     none
          par_x05    random   {lowerBounds.par_x05} {upperBounds.par_x05}  none   none     none
          par_x06    random   {lowerBounds.par_x06} {upperBounds.par_x06}  none   none     none
          par_x07    random   {lowerBounds.par_x07} {upperBounds.par_x07}  none   none     none
          par_x08    random   {lowerBounds.par_x08} {upperBounds.par_x08}  none   none     none
          par_x09    random   {lowerBounds.par_x09} {upperBounds.par_x09}  none   none     none
          par_x10    random   {lowerBounds.par_x10} {upperBounds.par_x10}  none   none     none
          par_x11    random   {lowerBounds.par_x11} {upperBounds.par_x11}  none   none     none
          par_x12    random   {lowerBounds.par_x12} {upperBounds.par_x12}  none   none     none
          par_x13    random   {lowerBounds.par_x13} {upperBounds.par_x13}  none   none     none
          par_x14    random   {lowerBounds.par_x14} {upperBounds.par_x14}  none   none     none
          par_x15    random   {lowerBounds.par_x15} {upperBounds.par_x15}  none   none     none
          par_x16    random   {lowerBounds.par_x16} {upperBounds.par_x16}  none   none     none
          par_x17    random   {lowerBounds.par_x17} {upperBounds.par_x17}  none   none     none
          par_x18    random   {lowerBounds.par_x18} {upperBounds.par_x18}  none   none     none
          par_x19    random   {lowerBounds.par_x19} {upperBounds.par_x19}  none   none     none
          par_x20    random   {lowerBounds.par_x20} {upperBounds.par_x20}  none   none     none
          par_x21    random   {lowerBounds.par_x21} {upperBounds.par_x21}  none   none     none
          par_x22    random   {lowerBounds.par_x22} {upperBounds.par_x22}  none   none     none
          par_x23    random   {lowerBounds.par_x23} {upperBounds.par_x23}  none   none     none
          par_x24    random   {lowerBounds.par_x24} {upperBounds.par_x24}  none   none     none
          par_x25    random   {lowerBounds.par_x25} {upperBounds.par_x25}  none   none     none
          par_x26    random   {lowerBounds.par_x26} {upperBounds.par_x26}  none   none     none
          par_x27    random   {lowerBounds.par_x27} {upperBounds.par_x27}  none   none     none
          par_x28    random   {lowerBounds.par_x28} {upperBounds.par_x28}  none   none     none
          par_x29    random   {lowerBounds.par_x29} {upperBounds.par_x29}  none   none     none
          par_x30    random   {lowerBounds.par_x30} {upperBounds.par_x30}  none   none     none
          par_x31    random   {lowerBounds.par_x31} {upperBounds.par_x31}  none   none     none
          par_x32    random   {lowerBounds.par_x32} {upperBounds.par_x32}  none   none     none
          par_x33    random   {lowerBounds.par_x33} {upperBounds.par_x33}  none   none     none
          par_x34    random   {lowerBounds.par_x34} {upperBounds.par_x34}  none   none     none
          par_x35    random   {lowerBounds.par_x35} {upperBounds.par_x35}  none   none     none
          par_r01    random   {lowerBounds.par_r01} {upperBounds.par_r01}  none   none     none
          par_r02    random   {lowerBounds.par_r02} {upperBounds.par_r02}  none   none     none
          par_r03    random   {lowerBounds.par_r03} {upperBounds.par_r03}  none   none     none
          par_r04    random   {lowerBounds.par_r04} {upperBounds.par_r04}  none   none     none
          par_r05    random   {lowerBounds.par_r05} {upperBounds.par_r05}  none   none     none
          par_r06    random   {lowerBounds.par_r06} {upperBounds.par_r06}  none   none     none
          par_r07    random   {lowerBounds.par_r07} {upperBounds.par_r07}  none   none     none
          par_r08    random   {lowerBounds.par_r08} {upperBounds.par_r08}  none   none     none
        EndParams

        BeginTiedParams
          # ---------------------------------------------------------------
          # MAX_MELT_FACTOR > MIN_MELT_FACTOR
          #
          # sum_x24_x25 = par_x24 + par_x25
          # Xtied =(c3 * X1 * X2) + (c2 * X2) + (c1 * X1) + c0
          # --> c0 = 0.0
          # --> c1 = 1.0
          # --> c2 = 1.0
          # --> c3 = 0.0
          #
          sum_x24_x25 2 par_x24 par_x25 linear 0.00 1.00 1.00 0.00 free
          #
          # ---------------------------------------------------------------
          # SNOW_SWI_MAX > SNOW_SWI_MIN
          #
          # sum_x13_x14 = par_x13 + par_x14
          # Xtied =(c3 * X1 * X2) + (c2 * X2) + (c1 * X1) + c0
          # --> c0 = 0.0
          # --> c1 = 1.0
          # --> c2 = 1.0
          # --> c3 = 0.0
          #
          sum_x13_x14 2 par_x13 par_x14 linear 0.00 1.00 1.00 0.00 free
          #
          # ---------------------------------------------------------------
          # FIELD_CAPACITY > SAT_WILT
          #
          # sum_x09_x10 = par_x09 + par_x10
          # Xtied =(c3 * X1 * X2) + (c2 * X2) + (c1 * X1) + c0
          # --> c0 = 0.0
          # --> c1 = 1.0
          # --> c2 = 1.0
          # --> c3 = 0.0
          #
          sum_x09_x10 2 par_x09 par_x10 linear 0.00 1.00 1.00 0.00 free
          #
          # ---------------------------------------------------------------
          # half the value but in [mm] not [m]
          #
          # half_x29 = par_x29 * 0.5 * 1000  --> half of it but in [mm] not [m]
          # Xtied = (c1 * X) + c0
          # --> c0 = 0.0
          # --> c1 = 500.
          #
          half_x29 1 par_x29 linear 500.0 0.0 free
          #
          # ---------------------------------------------------------------
          # half the value but in [mm] not [m]
          #
          # half_x30 = par_x30 * 0.5 * 1000  --> half of it but in [mm] not [m]
          # Xtied = (c1 * X) + c0
          # --> c0 = 0.0
          # --> c1 = 500.
          #
          half_x30 1 par_x30 linear 500.0 0.0 free
          #
          # ---------------------------------------------------------------
          # BASEFLOW_COEFF TOPSOIL  = 10.0^x4
          #
          # pow_x04 = 10.0**(par_x4)
          # Xtied = c2 * base ** (c1 * X) + c0
          # --> c0   = 0.0
          # --> c1   = 1.0
          # --> c2   = 1.0
          # --> base = 10.0
          #
          pow_x04 1 par_x04 exp 10.0 1.0 1.0 0.0 free
          #
          # ---------------------------------------------------------------
          # BASEFLOW_COEFF PHREATIC = 10.0^x11
          #
          # pow_x11 = 10.0**(par_x3)
          # Xtied = c2 * base ** (c1 * X) + c0
          # --> c0   = 0.0
          # --> c1   = 1.0
          # --> c2   = 1.0
          # --> base = 10.0
          #
          pow_x11 1 par_x11 exp 10.0 1.0 1.0 0.0 free
        EndTiedParams

        BeginResponseVars
          #name   filename                              keyword         line    col     token
          RawMetric  ./model/output/{run_name}-{run_index}_Diagnostics.csv;       OST_NULL        1       3       ','
        EndResponseVars

        BeginTiedRespVars
        # <name1> <np1> <pname 1,1 > <pname 1,2 > ... <pname 1,np1 > <type1> <type_data1>
          Metric 1 RawMetric wsum {evaluation_metric_multiplier}
        EndTiedRespVars

        BeginGCOP
          CostFunction Metric
          PenaltyFunction APM
        EndGCOP

        BeginConstraints
                # not needed when no constraints, but PenaltyFunction statement above is required
                # name     type     penalty    lwr   upr   resp.var
        EndConstraints

        # Randomsed control added
        {random_seed}

        #Algorithm should be last in this file:

        BeginDDSAlg
                PerturbationValue 0.20
                MaxIterations {max_iterations}
                UseRandomParamValues
                # UseInitialParamValues
                # above intializes DDS to parameter values IN the initial model input files
        EndDDSAlg
        """
        self.config.ost.set_tmpl(ost_tmpl)

    def derived_parameters(self):
        """Derived parameters are computed by Ostrich."""
        self.config.rvt.rain_correction = "par_x33"
        self.config.rvt.snow_correction = "par_x34"
