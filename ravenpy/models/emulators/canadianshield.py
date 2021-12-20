from collections import defaultdict
from pathlib import Path
from typing import cast

import xarray as xr
from pydantic.dataclasses import dataclass

from ravenpy.config import ConfigError
from ravenpy.config.commands import (
    HRU,
    LU,
    BaseDataCommand,
    BasinIndexCommand,
    HRUState,
    Sub,
)
from ravenpy.config.rvs import RVH, RVI
from ravenpy.models.base import Ostrich, Raven


class CANADIANSHIELD(Raven):
    """
    The Canadian Shield model is a useful configuration of Raven for
    Canadian shield basins characterized by shallow soils atop rock,
    with ample exposed rock and lakes.
    """

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

    @dataclass
    class HRU_ORGANIC(HRU):
        hru_id: int = 1
        land_use_class: str = "FOREST"
        veg_class: str = "FOREST"
        soil_profile: str = "SOILP_ORG"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"

    @dataclass
    class HRU_BEDROCK(HRU):
        hru_id: int = 2
        land_use_class: str = "FOREST"
        veg_class: str = "FOREST"
        soil_profile: str = "SOILP_BEDROCK"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "canadianshield")
        super().__init__(*args, **kwds)

        self.config.update(
            hrus=(CANADIANSHIELD.HRU_ORGANIC(), CANADIANSHIELD.HRU_BEDROCK()),
            subbasins=(
                Sub(
                    subbasin_id=1,
                    name="sub_001",
                    downstream_id=-1,
                    profile="None",
                    gauged=True,
                ),
            ),
        )

        #########
        # R V P #
        #########

        rvp_tmpl = """
        # tied parameters:
        # (it is important for OSTRICH to find every parameter place holder somewhere in this file)
        # (without this some parameters wouldn't be detectable)
        #    para_sum_x05_x06 = sum_x05_x06 = {params.par_x05} + {params.par_x06} = para_x05 + para_x06
        #    para_sum_x16_x17 = sum_x16_x17 = {params.par_x16} + {params.par_x17} = para_x16 + para_x17
        #    para_pow_x08     = par_pow_x08 = 10^({params.par_x08}) = 10^para_x08
        #    para_pow_x09     = par_pow_x09 = 10^({params.par_x09}) = 10^para_x09

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
          SOILP_ORG,     3, TOPSOIL, {params.par_x01},  VADOSE, {params.par_x02},  FRACBEDROCK, {params.par_x03},
          SOILP_BEDROCK, 3, TOPSOIL, {params.par_x01},  VADOSE,             0.00,  FRACBEDROCK, {params.par_x03},
        # SOILP_ORG,     3, TOPSOIL, para_x01, VADOSE, para_x02, FRACBEDROCK, para_x03,
        # SOILP_BEDROCK, 3, TOPSOIL, para_x01, VADOSE, 0.00,     FRACBEDROCK, para_x03,
        :EndSoilProfiles

        #-----------------------------------------------------------------
        # Land Use Classes
        #-----------------------------------------------------------------
        {land_use_classes}

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
        :GlobalParameter          SNOW_SWI {params.par_x15}            # para_x15
        :GlobalParameter      SNOW_SWI_MIN {params.par_x16}            # para_x16
        :GlobalParameter      SNOW_SWI_MAX {params.par_x17}            # para_x16+para_x17
        :GlobalParameter  SWI_REDUCT_COEFF {params.par_x18}            # para_x18
        :GlobalParameter     RAINSNOW_TEMP {params.par_x19}            # para_x19
        :GlobalParameter    RAINSNOW_DELTA {params.par_x20}            # para_x20
        :GlobalParameter   MAX_SWE_SURFACE {params.par_x21}            # para_x21
        :GlobalParameter    TOC_MULTIPLIER {params.par_x22}            # para_x22

        #-----------------------------------------------------------------
        # Soil Parameters
        #-----------------------------------------------------------------
        :SoilParameterList
          :Parameters,        POROSITY,         HBV_BETA,   BASEFLOW_COEFF,       BASEFLOW_N,  MAX_INTERFLOW_RATE,   FIELD_CAPACITY,         SAT_WILT,    MAX_PERC_RATE,   PET_CORRECTION,
               :Units,               -,                -,                -,                -,                   -,                -,                -,                -,                -,
              TOPSOIL,             1.0, {params.par_x07},              0.0,              0.0,    {params.par_x12}, {params.par_x06}, {params.par_x05}, {params.par_x13}, {params.par_x04},
               VADOSE,             1.0,              0.0, {params.par_x08}, {params.par_x10},                 0.0,              0.0,              0.0, {params.par_x14},              0.0,
          FRACBEDROCK,             1.0,              0.0, {params.par_x09}, {params.par_x11},                 0.0,              0.0,              0.0,              0.0,              0.0,
        #     TOPSOIL,             1.0,         para_x07,              0.0,              0.0,            para_x12, para_sum_x05_x06,         para_x05,         para_x13,         para_x04,
        #      VADOSE,             1.0,              0.0,     para_pow_x08,         para_x10,                 0.0,              0.0,              0.0,         para_x14,              0.0,
        # FRACBEDROCK,             1.0,              0.0,     para_pow_x09,         para_x11,                 0.0,              0.0,              0.0,              0.0,              0.0,
        :EndSoilParameterList

        # note: TOPSOIL FIELD_CAPACITY calculated as {params.par_x05} + {params.par_x06}
        # note: TOPSOIL BASEFLOW_COEFF calculated as  10^({params.par_x08})
        # note: PHREATIC BASEFLOW_COEFF calculated as  10^({params.par_x09})

        #-----------------------------------------------------------------
        # Land Use Parameters
        #-----------------------------------------------------------------
        :LandUseParameterList
          :Parameters, FOREST_SPARSENESS,      MELT_FACTOR,     DD_MELT_TEMP,  REFREEZE_FACTOR,          DEP_MAX,      OW_PET_CORR,
               :Units,                 -,                -,                -,                -,                -,                -,
               FOREST,               0.0, {params.par_x25}, {params.par_x24}, {params.par_x23}, {params.par_x26}, {params.par_x27},
        #      FOREST,               0.0,         para_x25,         para_x24,         para_x23,         para_x26,         para_x27,
        :EndLandUseParameterList

        #-----------------------------------------------------------------
        # Vegetation Parameters
        #-----------------------------------------------------------------
        :VegetationParameterList
          :Parameters,  SVF_EXTINCTION,   SAI_HT_RATIO,  RAIN_ICEPT_FACT,  SNOW_ICEPT_FACT,     MAX_CAPACITY, MAX_SNOW_CAPACITY,
               :Units,               -,              -,                -,                -,                -,                 -,
            FOREST,                0.5,            1.0, {params.par_x28}, {params.par_x29}, {params.par_x30},  {params.par_x31},
        #   FOREST,                0.5,            1.0,         para_x28,         para_x29,         para_x30,          para_x31,
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
        """
        self.config.rvp.set_tmpl(rvp_tmpl)

        #########
        # R V I #
        #########

        rvi_tmpl = """
        :Routing               {routing}
        :CatchmentRoute        ROUTE_TRI_CONVOLUTION
        :Evaporation           {evaporation}
        :OW_Evaporation        {ow_evaporation}
        :SWCanopyCorrect       SW_CANOPY_CORR_STATIC
        :RainSnowFraction      {rain_snow_fraction}
        :PotentialMeltMethod   POTMELT_DEGREE_DAY
        :PrecipIceptFract      PRECIP_ICEPT_LAI
        :SoilModel             SOIL_MULTILAYER 3

        :MonthlyInterpolationMethod MONTHINT_LINEAR_MID

        :Alias TOPSOIL          SOIL[0]
        :Alias VADOSE           SOIL[1]
        :Alias FRACBEDROCK      SOIL[2]

        :HydrologicProcesses
            :SnowRefreeze             FREEZE_DEGREE_DAY    SNOW_LIQ        SNOW
            :Precipitation            PRECIP_RAVEN         ATMOS_PRECIP    MULTIPLE
            :CanopyEvaporation        CANEVP_MAXIMUM       CANOPY          ATMOSPHERE
            :CanopySnowEvap           CANEVP_MAXIMUM       CANOPY_SNOW     ATMOSPHERE
            :SnowBalance              SNOBAL_TWO_LAYER     MULTIPLE        MULTIPLE
            :Abstraction              ABST_FILL            PONDED_WATER    DEPRESSION
            :OpenWaterEvaporation     OPEN_WATER_EVAP      DEPRESSION      ATMOSPHERE
            :Infiltration             INF_HBV              PONDED_WATER    MULTIPLE
            :Interflow                INTERFLOW_PRMS       TOPSOIL         SURFACE_WATER
            :Baseflow                 BASE_POWER_LAW       VADOSE          SURFACE_WATER
                  :-->Conditional HRU_TYPE IS SOIL_ORG
            :Baseflow                 BASE_POWER_LAW       FRACBEDROCK     SURFACE_WATER
                :Percolation          PERC_GAWSER          TOPSOIL         VADOSE
                  :-->Conditional HRU_TYPE IS SOIL_ORG
            :Percolation              PERC_GAWSER          VADOSE          FRACBEDROCK
                  :-->Conditional HRU_TYPE IS SOIL_ORG
                :Percolation          PERC_GAWSER          TOPSOIL         FRACBEDROCK
                  :-->Conditional HRU_TYPE IS SOIL_BEDROCK
            :SoilEvaporation          SOILEVAP_ROOT        TOPSOIL         ATMOSPHERE
        :EndHydrologicProcesses
        """

        self.config.rvi.set_tmpl(rvi_tmpl)

        self.config.rvi.routing = RVI.RoutingOptions.NONE
        self.config.rvi.rain_snow_fraction = RVI.RainSnowFractionOptions.DINGMAN
        self.config.rvi.evaporation = RVI.EvaporationOptions.PET_HARGREAVES_1985
        self.config.rvi.ow_evaporation = RVI.EvaporationOptions.PET_HARGREAVES_1985

    def derived_parameters(self):
        params = cast(CANADIANSHIELD.Params, self.config.rvp.params)

        if len(self.config.rvh.hrus) != 2:
            raise ConfigError("CANADIANSHIELD must have exactly two HRUs")

        if self.config.rvh.hrus[0].area != self.config.rvh.hrus[1].area:
            raise ConfigError("CANADIANSHIELD HRUs must have equal areas")

        self.config.rvh.hrus[0].area *= params.par_x34
        self.config.rvh.hrus[1].area *= 1.0 - params.par_x34

        soil0 = params.par_x01 * 1000.0 * 0.5
        soil1 = params.par_x02 * 1000.0 * 0.5
        soil2 = params.par_x03 * 1000.0 * 0.5
        self.config.rvc.set_hru_state(
            HRUState(index=1, soil0=soil0, soil1=soil1, soil2=soil2)
        )
        self.config.rvc.set_hru_state(
            HRUState(index=2, soil0=soil0, soil1=0, soil2=soil2)
        )

        self.config.rvt.rain_correction = params.par_x32
        self.config.rvt.snow_correction = params.par_x33


class CANADIANSHIELD_OST(Ostrich, CANADIANSHIELD):

    ostrich_to_raven_param_conversion = {
        "par_sum_x05_x06": "par_x06",
        "par_sum_x16_x17": "par_x17",
        "par_pow_x08": "par_x08",
        "par_pow_x09": "par_x09",
    }

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "canadianshield-ost")
        super().__init__(*args, **kwds)

        self.config.update(
            algorithm="DDS",
            max_iterations=50,
            suppress_output=True,
        )

        ####################
        # R V P (OST TMPL) #
        ####################

        rvp_tmpl = """
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
        #   FOREST,                0.0,   forest_frac,
        # :EndLandUseClasses
        {land_use_classes}

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
        """
        self.config.rvp.set_tmpl(rvp_tmpl, is_ostrich=True)

        ####################
        # R V T (OST TMPL) #
        ####################

        self.config.rvt.set_tmpl(is_ostrich=True)

        self.config.rvt.update("rain_correction", "par_x32")
        self.config.rvt.update("snow_correction", "par_x33")

        ####################
        # R V H (OST TMPL) #
        ####################

        self.config.rvh.set_tmpl(is_ostrich=True)

        ####################
        # R V C (OST TMPL) #
        ####################

        self.config.rvc.set_tmpl(is_ostrich=True)

        ##########
        # O S T  #
        ##########

        ost_tmpl = """
        ProgramType         DDS
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
          {identifier}.rvp.tpl;    {identifier}.rvp
          {identifier}.rvc.tpl;    {identifier}.rvc
          {identifier}.rvh.tpl;    {identifier}.rvh
          {identifier}.rvt.tpl;    {identifier}.rvt
          #can be multiple (.rvh, .rvi)
        EndFilePairs

        #Parameter/DV Specification
        BeginParams
          #parameter init.    low                   high                    tx_in  tx_ost  tx_out
          par_x01    random   {lowerBounds.par_x01} {upperBounds.par_x01}   none   none    none # SoilProfiles:            thickness TOPSOIL (m)
          par_x02    random   {lowerBounds.par_x02} {upperBounds.par_x02}   none   none    none # SoilProfiles:            thickness VADOSE (m) in organic soil class only
          par_x03    random   {lowerBounds.par_x03} {upperBounds.par_x03}   none   none    none # SoilProfiles:            thickness FRACBEDROCK (m)
          par_x04    random   {lowerBounds.par_x04} {upperBounds.par_x04}   none   none    none # SoilParameterList:       PET_CORRECTION TOPSOIL
          par_x05    random   {lowerBounds.par_x05} {upperBounds.par_x05}   none   none    none # SoilParameterList:       SAT_WILT TOPSOIL
          par_x06    random   {lowerBounds.par_x06} {upperBounds.par_x06}   none   none    none # SoilParameterList:       FIELD_CAPACITY TOPSOIL = SAT_WILT TOPSOIL + x05
          par_x07    random   {lowerBounds.par_x07} {upperBounds.par_x07}   none   none    none # SoilParameterList:       HBV_BETA TOPSOIL
          par_x08    random   {lowerBounds.par_x08} {upperBounds.par_x08}   none   none    none # SoilParameterList:       BASEFLOW_COEFF VADOSE = 10^x08
          par_x09    random   {lowerBounds.par_x09} {upperBounds.par_x09}   none   none    none # SoilParameterList:       BASEFLOW_COEFF FRACBEDROCK = 10^x09
          par_x10    random   {lowerBounds.par_x10} {upperBounds.par_x10}   none   none    none # SoilParameterList:       BASEFLOW_N VADOSE
          par_x11    random   {lowerBounds.par_x11} {upperBounds.par_x11}   none   none    none # SoilParameterList:       BASEFLOW_N FRACBEDROCK
          par_x12    random   {lowerBounds.par_x12} {upperBounds.par_x12}   none   none    none # SoilParameterList:       MAX_INTERFLOW_RATE TOPSOIL
          par_x13    random   {lowerBounds.par_x13} {upperBounds.par_x13}   none   none    none # SoilParameterList:       MAX_PERC_RATE TOPSOIL
          par_x14    random   {lowerBounds.par_x14} {upperBounds.par_x14}   none   none    none # SoilParameterList:       MAX_PERC_RATE VADOSE
          par_x15    random   {lowerBounds.par_x15} {upperBounds.par_x15}   none   none    none # GlobalParameter:         SNOW_SWI
          par_x16    random   {lowerBounds.par_x16} {upperBounds.par_x16}   none   none    none # GlobalParameter:         SNOW_SWI_MIN
          par_x17    random   {lowerBounds.par_x17} {upperBounds.par_x17}   none   none    none # GlobalParameter:         SNOW_SWI_MAX = SNOW_SWI_MIN + x16
          par_x18    random   {lowerBounds.par_x18} {upperBounds.par_x18}   none   none    none # GlobalParameter:         SWI_REDUCT_COEFF
          par_x19    random   {lowerBounds.par_x19} {upperBounds.par_x19}   none   none    none # GlobalParameter:         RAINSNOW_TEMP (dC)
          par_x20    random   {lowerBounds.par_x20} {upperBounds.par_x20}   none   none    none # GlobalParameter:         RAINSNOW_DELTA (dC)
          par_x21    random   {lowerBounds.par_x21} {upperBounds.par_x21}   none   none    none # GlobalParameter:         MAX_SWE_SURFACE (mm)
          par_x22    random   {lowerBounds.par_x22} {upperBounds.par_x22}   none   none    none # GlobalParameter:         TOC_MULTIPLIER
          par_x23    random   {lowerBounds.par_x23} {upperBounds.par_x23}   none   none    none # LandUseParameterList:    REFREEZE_FACTOR
          par_x24    random   {lowerBounds.par_x24} {upperBounds.par_x24}   none   none    none # LandUseParameterList:    DD_MELT_TEMP
          par_x25    random   {lowerBounds.par_x25} {upperBounds.par_x25}   none   none    none # LandUseParameterList:    MELT_FACTOR
          par_x26    random   {lowerBounds.par_x26} {upperBounds.par_x26}   none   none    none # LandUseParameterList:    DEP_MAX (mm)
          par_x27    random   {lowerBounds.par_x27} {upperBounds.par_x27}   none   none    none # LandUseParameterList:    OW_PET_CORR
          par_x28    random   {lowerBounds.par_x28} {upperBounds.par_x28}   none   none    none # VegetationParameterList: RAIN_ICEPT_FACT
          par_x29    random   {lowerBounds.par_x29} {upperBounds.par_x29}   none   none    none # VegetationParameterList: SNOW_ICEPT_FACT
          par_x30    random   {lowerBounds.par_x30} {upperBounds.par_x30}   none   none    none # VegetationParameterList: MAX_CAPACITY (mm)
          par_x31    random   {lowerBounds.par_x31} {upperBounds.par_x31}   none   none    none # VegetationParameterList: MAX_SNOW_CAPACITY (mm)
          par_x32    random   {lowerBounds.par_x32} {upperBounds.par_x32}   none   none    none # Gauge:                   RAINCORRECTION
          par_x33    random   {lowerBounds.par_x33} {upperBounds.par_x33}   none   none    none # Gauge:                   SNOWCORRECTION
          par_x34    random   {lowerBounds.par_x34} {upperBounds.par_x34}   none   none    none # HRUs:                    ratio to derive fraction of organic soil (x34*AREA) and bedrock ((1-x34)*AREA)
        EndParams

        BeginTiedParams
          #
          par_sum_x05_x06  2 par_x05 par_x06 linear 0.00 1.00 1.00 0.00 free
          par_sum_x16_x17  2 par_x16 par_x17 linear 0.00 1.00 1.00 0.00 free
          par_pow_x08      1 par_x08         exp    10.0 1.00 1.00 0.00 free
          par_pow_x09      1 par_x09         exp    10.0 1.00 1.00 0.00 free
          par_half_x01     1 par_x01         linear 500.0     0.00      free
          par_half_x02     1 par_x02         linear 500.0     0.00      free
          par_half_x03     1 par_x03         linear 500.0     0.00      free
          #
          # AREA HRU1=SOILP_ORG (km2)
          # soil organic area = x34 * AREA
          # area SB1 = {area}
          par_area_organic_SB1 1 par_x34         linear  {area}      0.0    free
          #
          # AREA HRU2=SOILP_BEDROCK (km2)
          # soil bedrock area = (1-x34)*AREA = AREA - x34 * AREA
          # area SB1 = {area}
          par_area_bedrock_SB1 1 par_x34         linear -{area}   {area}    free
          #
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
        RandomSeed 0

        #Algorithm should be last in this file:

        BeginDDSAlg
                PerturbationValue 0.20
                MaxIterations 10
                UseRandomParamValues
                # UseInitialParamValues
                # above intializes DDS to parameter values IN the initial model input files
        EndDDSAlg
        """
        self.config.ost.set_tmpl(ost_tmpl)

    def derived_parameters(self):

        self.config.ost.set_extra_attributes(area=self.config.rvh.hrus[0].area)

        # Here we are abusing the pydantic.dataclass type checking
        # mechanism, in virtue of the fact that it's only operating
        # when passing the args to the constructor (not when setting
        # them manually after the object creation, like we do here):
        # the following attributes are all numeric in nature, but in
        # the context of Ostrich templating we need to set them to
        # non-numeric values.

        self.config.rvh.hrus[0].area = "par_area_organic_SB1"  # type: ignore
        self.config.rvh.hrus[1].area = "par_area_bedrock_SB1"  # type: ignore

        self.config.rvc.set_hru_state(HRUState(index=1))
        self.config.rvc.hru_states[1].soil0 = "par_half_x01"  # type: ignore
        self.config.rvc.hru_states[1].soil1 = "par_half_x02"  # type: ignore
        self.config.rvc.hru_states[1].soil2 = "par_half_x03"  # type: ignore

        self.config.rvc.set_hru_state(HRUState(index=2))
        self.config.rvc.hru_states[2].soil0 = "par_half_x01"  # type: ignore
        self.config.rvc.hru_states[2].soil2 = "par_half_x03"  # type: ignore
