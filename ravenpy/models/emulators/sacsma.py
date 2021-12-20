from collections import defaultdict
from pathlib import Path
from typing import cast

from pydantic.dataclasses import dataclass

from ravenpy.config.commands import HRU, LU, BasinIndexCommand, HRUState, Sub
from ravenpy.config.rvs import RVI
from ravenpy.models.base import Ostrich, Raven


class SACSMA(Raven):
    """
    Details about SAC-SMA:
    https://wiki.ewater.org.au/display/SD41/Sacramento+Model+-+SRG
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

    @dataclass
    class HRU(HRU):  # type: ignore
        land_use_class: str = "FOREST"
        veg_class: str = "FOREST"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "sacsma")
        super().__init__(*args, **kwds)

        self.config.update(
            hrus=(SACSMA.HRU(),),
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
        # --------------------------
        # Parameters and description
        # (adopted from Table 1 under https://wiki.ewater.org.au/display/SD41/Sacramento+Model+-+SRG)
        # --------------------------
        # para_x01 0.001  0.015    # LZPK;           The ratio of water in LZFPM, which drains as base flow each day.; fraction; default=0.01
        # para_x02 0.03   0.2      # LZSK;           The ratio of water in LZFSM which drains as base flow each day.; fraction; default=0.05
        # para_x03 0.2    0.5      # UZK;            The fraction of water in UZFWM, which drains as interflow each day.; fraction; default=0.3
        # para_x04 0.025  0.125    # UZTWM;          Upper Zone Tension Water Maximum. The maximum volume of water held by the upper zone between
        #                                            field capacity and the wilting point which can be lost by direct evaporation and evapotranspiration
        #                                            from soil surface. This storage is filled before any water in the upper zone is transferred to
        #                                            other storages.; m; default=0.050
        # para_x05 0.010  0.075    # UZFWM;          Upper Zone Free Water Maximum, this storage is the source of water for interflow and the driving
        #                                            force for transferring water to deeper depths.; m; default=0.040
        # para_x06 0.075  0.300    # LZTWM;          Lower Zone Tension Water Maximum, the maximum capacity of lower zone tension water. Water from
        #                                            this store can only be removed through evapotranspiration.; m; default=0.130
        # para_x07 0.015  0.300    # LZFSM;          Lower Zone Free Water Supplemental Maximum, the maximum volume from which supplemental base flow
        #                                            can be drawn.; m; default=0.025
        # para_x08 0.040  0.600    # LZFPM;          Lower Zone Free Water Primary Maximum, the maximum capacity from which primary base flow can be
        #                                            drawn.; m; default=0.060
        # para_x09 0.0    0.5      # PFREE;          The minimum proportion of percolation from the upper zone to the lower zone directly available
        #                                            for recharging the lower zone free water stores.; percent/100; default=0.06
        # para_x10 0.0    3.0      # REXP;           An exponent determining the rate of change of the percolation rate with changing lower zone water
        #                                            storage.; none; default=1.0
        # para_x11 0.0    80       # ZPERC;          The proportional increase in Pbase that defines the maximum percolation rate.;
        #                                            none; default=40
        # para_x12 0.0    0.8      # SIDE;           The ratio of non-channel baseflow (deep recharge) to channel (visible) baseflow.; ratio;
        #                                            default=0.0
        # n/a      0.0    0.1      # SSOUT;          The volume of the flow which can be conveyed by porous material in the bed of stream.; mm;
        #                                            default=0.0
        # para_x13 0.0    0.05     # PCTIM;          The permanently impervious fraction of the basin contiguous with stream channels, which
        #                                            contributes to direct runoff.; percent/100; default=0.01
        # para_x14 0.0    0.2      # ADIMP;          The additional fraction of the catchment which develops impervious characteristics under
        #                                            soil saturation conditions.; percent/100; default=0.0
        # para_x15 0.0    0.1      # RIVA=SARVA;     A decimal fraction representing that portion of the basin normally covered by streams,
        #                                            lakes and vegetation that can deplete stream flow by evapotranspiration.; percent/100;
        #                                            default=0.0
        # para_x16 0.0    0.4      # RSERV;          Fraction of lower zone free water unavailable for transpiration; percent/100; default=0.3
        # n/a      0.0    1.0      # UH1;            The first  component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # n/a      0.0    1.0      # UH2;            The second component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # n/a      0.0    1.0      # UH3;            The third  component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # n/a      0.0    1.0      # UH4;            The fourth component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # n/a      0.0    1.0      # UH5;            The fifth  component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # para_x17 0.0    8.0      # MELT_FACTOR;    maximum snow melt factor used in degree day models (not in original SAC-SMA
        #                                            model); mm/d/C; default=1.5
        # para_x18 0.3    20.0     # GAMMA_SHAPE;    used to build unit hydrograph, LandUseParameterList; none ; default=1.0
        # para_x19 0.01   5.0      # GAMMA_SCALE;    used to build unit hydrograph, LandUseParameterList; none ; default=1.0
        # para_x20 0.8    1.2      # RAINCORRECTION; Muliplier to correct rain, Gauge; none ; default=1.0
        # para_x21 0.8    1.2      # SNOWCORRECTION; Muliplier to correct snow, Gauge; none ; default=1.0

        # tied parameters:
        # (it is important for OSTRICH to find every parameter place holder somewhere in this file)
        # (without this parameters that are used to derive parameters wouldn't be detectable by Ostrich)
        #    para_x04_mm       = par_soil0_mm     = para_x04 * 1000       =  par_x04 * 1000
        #    para_x06_mm       = par_soil0_mm     = para_x06 * 1000       =  par_x06 * 1000
        #    para_bf_loss_frac = PAR_BF_LOSS_FRAC = para_x12/(1+para_x12) = par_x12/(1+par_x12)
        #    para_pow_x01      = POW_X01          = 10^(para_x01)         = 10^par_x01
        #    para_pow_x02      = POW_X02          = 10^(para_x02)         = 10^par_x02
        #    para_pow_x03      = POW_X03          = 10^(para_x03)         = 10^par_x03

        #-----------------------------------------------------------------
        # Soil Classes
        #-----------------------------------------------------------------
        :SoilClasses
          :Attributes,
          :Units,
            UZT,
            UZF,
            LZT,
            LZFP,
            LZFS,
            ADIM,
            GW
        :EndSoilClasses

        #-----------------------------------------------------------------
        # Land Use Classes
        #-----------------------------------------------------------------
        :LandUseClasses,
          :Attributes,           IMPERM,    FOREST_COV,
               :Units,             frac,          frac,
               FOREST, {params.par_x13},           1.0,
        #      FOREST,          <PCTIM>,           1.0,
        #      FOREST,         para_x13,           1.0,
        #
        :EndLandUseClasses

        #-----------------------------------------------------------------
        # Vegetation Classes
        #-----------------------------------------------------------------
        :VegetationClasses,
          :Attributes,        MAX_HT,       MAX_LAI, MAX_LEAF_COND,
               :Units,             m,          none,      mm_per_s,
               FOREST,             4,             5,             5,
        #      FOREST,             4,             5,             5,
        :EndVegetationClasses

        #-----------------------------------------------------------------
        # Soil Profiles
        #-----------------------------------------------------------------
        :SoilProfiles
               LAKE, 0
          DEFAULT_P, 7, UZT, {params.par_x04}, UZF, {params.par_x05}, LZT,  {params.par_x06}, LZFP,  {params.par_x08}, LZFS, {params.par_x07}, ADIM, 100, GW, 100
        # DEFAULT_P, 7, UZT,          <UZTWM>, UZF,          <UZFWM>, LZT,           <LZTWM>, LZFP,           <LZFPM>, LZFS,          <LZFSM>, ADIM, 100, GW, 100
        # DEFAULT_P, 7, UZT,         para_x04, UZF,         para_x05, LZT,          para_x06, LZFP,          para_x08, LZFS,         para_x07, ADIM, 100, GW, 100
        :EndSoilProfiles

        #-----------------------------------------------------------------
        # Global Parameters
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Soil Parameters
        #-----------------------------------------------------------------
        :SoilParameterList
          :Parameters,        POROSITY,   SAC_PERC_ALPHA,   SAC_PERC_EXPON,   SAC_PERC_PFREE,
               :Units,             [-],              [-],              [-],           [0..1],
            [DEFAULT],             1.0, {params.par_x11}, {params.par_x10}, {params.par_x09},
        #   [DEFAULT],             1.0,          <ZPERC>,           <REXP>,          <PFREE>,
        #   [DEFAULT],             1.0,         para_x11,         para_x10,         para_x09,
        :EndSoilParameterList
        :SoilParameterList
          :Parameters,              BASEFLOW_COEFF,     UNAVAIL_FRAC
               :Units,                         1/d,             0..1
            [DEFAULT],                         0.0,              0.0
               UZF,               {params.par_x03},              0.0
               LZFP,              {params.par_x01}, {params.par_x16}
               LZFS,              {params.par_x02},              0.0
        #   [DEFAULT],                         0.0,              0.0
        #      UZF,                          <UZK>,              0.0
        #      LZFP,                        <LZPK>,          <RSERV>
        #      LZFS,                        <LZSK>,              0.0
        #   [DEFAULT],                         0.0,              0.0
        #      UZF,                    10^para_x03,              0.0
        #      LZFP,                   10^para_x01,         para_x16
        #      LZFS,                   10^para_x02,              0.0
        :EndSoilParameterList

        #-----------------------------------------------------------------
        # Land Use Parameters
        #-----------------------------------------------------------------
        :LandUseParameterList
          :Parameters,      GAMMA_SHAPE,      GAMMA_SCALE,      MELT_FACTOR,  STREAM_FRACTION, MAX_SAT_AREA_FRAC,      BF_LOSS_FRACTION
               :Units,                -,                -,           mm/d/C,           [0..1],            [0..1],                [0..1]
            [DEFAULT], {params.par_x18}, {params.par_x19}, {params.par_x17}, {params.par_x15},  {params.par_x14},      {params.par_x12}
        #                      para_x18,         para_x19,         para_x17,           <RIVA>,           <ADIMP>,     <SIDE>/(1+<SIDE>)
        #                      para_x18,         para_x19,         para_x17,         para_x15,          para_x14, para_x12/(1+para_x12)
        :EndLandUseParameterList

        #-----------------------------------------------------------------
        # Vegetation Parameters
        #-----------------------------------------------------------------
        :VegetationParameterList
          :Parameters,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,
               :Units,               -,               -,
            [DEFAULT],             0.0,             0.0,
        :EndVegetationParameterList
        """

        self.config.rvp.set_tmpl(rvp_tmpl)

        #########
        # R V I #
        #########

        rvi_tmpl = """
        :PotentialMeltMethod     POTMELT_DEGREE_DAY
        :RainSnowFraction        {rain_snow_fraction}
        :Evaporation             {evaporation}
        :CatchmentRoute          ROUTE_GAMMA_CONVOLUTION
        :Routing                 ROUTE_NONE

        :SoilModel               SOIL_MULTILAYER 7

        :Alias UZ_T  SOIL[0]
        :Alias UZ_F  SOIL[1]
        :Alias LZ_T  SOIL[2]
        :Alias LZ_PF SOIL[3]
        :Alias LZ_PS SOIL[4]

        :HydrologicProcesses
          :SnowBalance           SNOBAL_SIMPLE_MELT   SNOW          PONDED_WATER
          :Precipitation         RAVEN_DEFAULT        ATMOS_PRECIP  MULTIPLE
          :SoilEvaporation       SOILEVAP_SACSMA      MULTIPLE      ATMOSPHERE
          :SoilBalance           SOILBAL_SACSMA       MULTIPLE      MULTIPLE
          :OpenWaterEvaporation  OPEN_WATER_RIPARIAN  SURFACE_WATER ATMOSPHERE
        :EndHydrologicProcesses
        """
        self.config.rvi.set_tmpl(rvi_tmpl)

        self.config.rvi.rain_snow_fraction = RVI.RainSnowFractionOptions.DATA
        self.config.rvi.evaporation = RVI.EvaporationOptions.PET_OUDIN

    def derived_parameters(self):
        params = cast(SACSMA.Params, self.config.rvp.params)

        soil0 = params.par_x04 * 1000.0
        soil2 = params.par_x06 * 1000.0
        self.config.rvc.hru_states[1] = HRUState(soil0=soil0, soil2=soil2)

        self.config.rvt.rain_correction = params.par_x20
        self.config.rvt.snow_correction = params.par_x21


class SACSMA_OST(Ostrich, SACSMA):

    ostrich_to_raven_param_conversion = {
        "pow_x01": "par_x01",
        "pow_x02": "par_x02",
        "pow_x03": "par_x03",
        "par_bf_loss_frac": "par_x12",
    }

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "sacsma-ost")
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
        # initialize to full tension storage
        :UniformInitialConditions SOIL[0] par_soil0_mm #<UZTWM> = para_x04*1000
        :UniformInitialConditions SOIL[2] par_soil2_mm #<LZTWM> = para_x06*1000
        """
        self.config.rvc.set_tmpl(rvc_tmpl, is_ostrich=True)

        ####################
        # R V P (OST TMPL) #
        ####################

        rvp_tmpl = """
        # Parameters and description
        # (adopted from Table 1 under https://wiki.ewater.org.au/display/SD41/Sacramento+Model+-+SRG)
        # --------------------------
        # para_x01 0.001  0.015    # LZPK;           The ratio of water in LZFPM, which drains as base flow each day.; fraction; default=0.01
        # para_x02 0.03   0.2      # LZSK;           The ratio of water in LZFSM which drains as base flow each day.; fraction; default=0.05
        # para_x03 0.2    0.5      # UZK;            The fraction of water in UZFWM, which drains as interflow each day.; fraction; default=0.3
        # para_x04 0.025  0.125    # UZTWM;          Upper Zone Tension Water Maximum. The maximum volume of water held by the upper zone between
        #                                            field capacity and the wilting point which can be lost by direct evaporation and evapotranspiration
        #                                            from soil surface. This storage is filled before any water in the upper zone is transferred to
        #                                            other storages.; m; default=0.050
        # para_x05 0.010  0.075    # UZFWM;          Upper Zone Free Water Maximum, this storage is the source of water for interflow and the driving
        #                                            force for transferring water to deeper depths.; m; default=0.040
        # para_x06 0.075  0.300    # LZTWM;          Lower Zone Tension Water Maximum, the maximum capacity of lower zone tension water. Water from
        #                                            this store can only be removed through evapotranspiration.; m; default=0.130
        # para_x07 0.015  0.300    # LZFSM;          Lower Zone Free Water Supplemental Maximum, the maximum volume from which supplemental base flow
        #                                            can be drawn.; m; default=0.025
        # para_x08 0.040  0.600    # LZFPM;          Lower Zone Free Water Primary Maximum, the maximum capacity from which primary base flow can be
        #                                            drawn.; m; default=0.060
        # para_x09 0.0    0.5      # PFREE;          The minimum proportion of percolation from the upper zone to the lower zone directly available
        #                                            for recharging the lower zone free water stores.; percent/100; default=0.06
        # para_x10 0.0    3.0      # REXP;           An exponent determining the rate of change of the percolation rate with changing lower zone water
        #                                            storage.; none; default=1.0
        # para_x11 0.0    80       # ZPERC;          The proportional increase in Pbase that defines the maximum percolation rate.;
        #                                            none; default=40
        # para_x12 0.0    0.8      # SIDE;           The ratio of non-channel baseflow (deep recharge) to channel (visible) baseflow.; ratio;
        #                                            default=0.0
        # n/a      0.0    0.1      # SSOUT;          The volume of the flow which can be conveyed by porous material in the bed of stream.; mm;
        #                                            default=0.0
        # para_x13 0.0    0.05     # PCTIM;          The permanently impervious fraction of the basin contiguous with stream channels, which
        #                                            contributes to direct runoff.; percent/100; default=0.01
        # para_x14 0.0    0.2      # ADIMP;          The additional fraction of the catchment which develops impervious characteristics under
        #                                            soil saturation conditions.; percent/100; default=0.0
        # para_x15 0.0    0.1      # RIVA=SARVA;     A decimal fraction representing that portion of the basin normally covered by streams,
        #                                            lakes and vegetation that can deplete stream flow by evapotranspiration.; percent/100;
        #                                            default=0.0
        # para_x16 0.0    0.4      # RSERV;          Fraction of lower zone free water unavailable for transpiration; percent/100; default=0.3
        # n/a      0.0    1.0      # UH1;            The first  component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # n/a      0.0    1.0      # UH2;            The second component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # n/a      0.0    1.0      # UH3;            The third  component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # n/a      0.0    1.0      # UH4;            The fourth component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # n/a      0.0    1.0      # UH5;            The fifth  component of the unit hydrograph, i.e. the proportion of instantaneous runoff
        #                                            not lagged; percent/100; default=1.0
        # para_x17 0.0    8.0      # MELT_FACTOR;    maximum snow melt factor used in degree day models (not in original SAC-SMA
        #                                            model); mm/d/C; default=1.5
        # para_x18 0.3    20.0     # GAMMA_SHAPE;    used to build unit hydrograph, LandUseParameterList; none ; default=1.0
        # para_x19 0.01   5.0      # GAMMA_SCALE;    used to build unit hydrograph, LandUseParameterList; none ; default=1.0
        # para_x20 0.8    1.2      # RAINCORRECTION; Muliplier to correct rain, Gauge; none ; default=1.0
        # para_x21 0.8    1.2      # SNOWCORRECTION; Muliplier to correct snow, Gauge; none ; default=1.0

        # tied parameters:
        # (it is important for OSTRICH to find every parameter place holder somewhere in this file)
        # (without this parameters that are used to derive parameters wouldn't be detectable by Ostrich)
        #    para_x04_mm       = par_soil0_mm     = para_x04 * 1000       =  par_x04 * 1000
        #    para_x06_mm       = par_soil2_mm     = para_x06 * 1000       =  par_x06 * 1000
        #    para_bf_loss_frac = par_bf_loss_frac = para_x12/(1+para_x12) = par_x12/(1+par_x12)
        #    para_pow__x01     = pow_x01          = 10^(para_x01)         = 10^par_x01
        #    para_pow__x02     = pow_x02          = 10^(para_x02)         = 10^par_x02
        #    para_pow__x03     = pow_x03          = 10^(para_x03)         = 10^par_x03

        #-----------------------------------------------------------------
        # Soil Classes
        #-----------------------------------------------------------------
        :SoilClasses
          :Attributes,
          :Units,
            UZT,
            UZF,
            LZT,
            LZFP,
            LZFS,
            ADIM,
            GW
        :EndSoilClasses

        #-----------------------------------------------------------------
        # Land Use Classes
        #-----------------------------------------------------------------
        :LandUseClasses,
          :Attributes,           IMPERM,    FOREST_COV,
               :Units,             frac,          frac,
               FOREST,          par_x13,           1.0,
        #      FOREST,          <PCTIM>,           1.0,
        #      FOREST,         para_x13,           1.0,
        #
        :EndLandUseClasses

        #-----------------------------------------------------------------
        # Vegetation Classes
        #-----------------------------------------------------------------
        :VegetationClasses,
          :Attributes,        MAX_HT,       MAX_LAI, MAX_LEAF_COND,
               :Units,             m,          none,      mm_per_s,
               FOREST,             4,             5,             5,
        #      FOREST,             4,             5,             5,
        :EndVegetationClasses

        #-----------------------------------------------------------------
        # Soil Profiles
        #-----------------------------------------------------------------
        :SoilProfiles
               LAKE, 0
          DEFAULT_P, 7, UZT,          par_x04, UZF,          par_x05, LZT,           par_x06, LZFP,           par_x08, LZFS,          par_x07, ADIM, 100, GW, 100
        # DEFAULT_P, 7, UZT,          <UZTWM>, UZF,          <UZFWM>, LZT,           <LZTWM>, LZFP,           <LZFPM>, LZFS,          <LZFSM>, ADIM, 100, GW, 100
        # DEFAULT_P, 7, UZT,         para_x04, UZF,         para_x05, LZT,          para_x06, LZFP,          para_x08, LZFS,         para_x07, ADIM, 100, GW, 100
        :EndSoilProfiles

        #-----------------------------------------------------------------
        # Global Parameters
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Soil Parameters
        #-----------------------------------------------------------------
        :SoilParameterList
          :Parameters,        POROSITY,   SAC_PERC_ALPHA,   SAC_PERC_EXPON,   SAC_PERC_PFREE,
               :Units,             [-],              [-],              [-],           [0..1],
            [DEFAULT],             1.0,          par_x11,          par_x10,          par_x09,
        #   [DEFAULT],             1.0,          <ZPERC>,           <REXP>,          <PFREE>,
        #   [DEFAULT],             1.0,         para_x11,         para_x10,         para_x09,
        :EndSoilParameterList
        :SoilParameterList
          :Parameters,        BASEFLOW_COEFF,     UNAVAIL_FRAC
               :Units,                   1/d,             0..1
            [DEFAULT],                   0.0,              0.0
               UZF,                  pow_x03,              0.0
               LZFP,                 pow_x01,          par_x16
               LZFS,                 pow_x02,              0.0
        #   [DEFAULT],                   0.0,              0.0
        #      UZF,                    <UZK>,              0.0
        #      LZFP,                  <LZPK>,          <RSERV>
        #      LZFS,                  <LZSK>,              0.0
        #   [DEFAULT],                   0.0,              0.0
        #      UZF,              10^para_x03,              0.0
        #      LZFP,             10^para_x01,         para_x16
        #      LZFS,             10^para_x02,              0.0
        :EndSoilParameterList

        #-----------------------------------------------------------------
        # Land Use Parameters
        #-----------------------------------------------------------------
        :LandUseParameterList
          :Parameters,      GAMMA_SHAPE,      GAMMA_SCALE,      MELT_FACTOR,  STREAM_FRACTION, MAX_SAT_AREA_FRAC,      BF_LOSS_FRACTION
               :Units,                -,                -,           mm/d/C,           [0..1],            [0..1],                [0..1]
            [DEFAULT],          par_x18,          par_x19,          par_x17,          par_x15,           par_x14,      par_bf_loss_frac
        #                      para_x18,         para_x19,         para_x17,           <RIVA>,           <ADIMP>,     <SIDE>/(1+<SIDE>)
        #                      para_x18,         para_x19,         para_x17,         para_x15,          para_x14, para_x12/(1+para_x12)
        :EndLandUseParameterList

        #-----------------------------------------------------------------
        # Vegetation Parameters
        #-----------------------------------------------------------------
        :VegetationParameterList
          :Parameters,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,
               :Units,               -,               -,
            [DEFAULT],             0.0,             0.0,
        :EndVegetationParameterList
        """
        self.config.rvp.set_tmpl(rvp_tmpl, is_ostrich=True)

        ####################
        # R V T (OST TMPL) #
        ####################

        self.config.rvt.set_tmpl(is_ostrich=True)

        self.config.rvt.update("rain_correction", "par_x20")
        self.config.rvt.update("snow_correction", "par_x21")

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
          sacsma-ost.rvp.tpl;  sacsma-ost.rvp
          sacsma-ost.rvc.tpl;  sacsma-ost.rvc
          sacsma-ost.rvt.tpl;  sacsma-ost.rvt
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
        EndParams

        BeginTiedParams

          # ---------------------------------------------------------------
          # BASEFLOW_COEFF LZPK  = 10.0^x01
          #
          # pow_x01 = 10.0**(par_x01)
          # Xtied = c2 * base ** (c1 * X) + c0
          # --> c0   = 0.0
          # --> c1   = 1.0
          # --> c2   = 1.0
          # --> base = 10.0
          #
          pow_x01 1 par_x01 exp 10.0 1.0 1.0 0.0 free
          #
          # ---------------------------------------------------------------
          # BASEFLOW_COEFF LZSK = 10.0^x02
          #
          # pow_x02 = 10.0**(par_x02)
          # Xtied = c2 * base ** (c1 * X) + c0
          # --> c0   = 0.0
          # --> c1   = 1.0
          # --> c2   = 1.0
          # --> base = 10.0
          #
          pow_x02 1 par_x02 exp 10.0 1.0 1.0 0.0 free
          #
          # ---------------------------------------------------------------
          # BASEFLOW_COEFF UZK = 10.0^x03
          #
          # pow_x03 = 10.0**(par_x03)
          # Xtied = c2 * base ** (c1 * X) + c0
          # --> c0   = 0.0
          # --> c1   = 1.0
          # --> c2   = 1.0
          # --> base = 10.0
          #
          pow_x03 1 par_x03 exp 10.0 1.0 1.0 0.0 free
          #
          # ---------------------------------------------------------------
          # PAR_BF_LOSS_FRAC = x12/(1+x12)
          #
          # sum_x24_x25 = par_x24 + par_x25
          # Xtied =(c3*X1+c2)/(c1*X2+c0)
          # --> c0 = 1.0
          # --> c1 = 1.0
          # --> c2 = 0.0
          # --> c3 = 1.0
          #
          par_bf_loss_frac 2 par_x12 par_x12 ratio 1.00 0.00 1.00 1.00 free
          #
          # ---------------------------------------------------------------
          # SOIL[0] in [mm] not [m]
          #
          # par_soil0_mm = par_x04 * 1000
          # Xtied = (c1 * X) + c0
          # --> c0 = 0.0
          # --> c1 = 1000.
          #
          par_soil0_mm 1 par_x04 linear 1000.0 0.0 free
          #
          # ---------------------------------------------------------------
          # SOIL[2] in [mm] not [m]
          #
          # par_soil2_mm = par_x06 * 1000
          # Xtied = (c1 * X) + c0
          # --> c0 = 0.0
          # --> c1 = 1000.
          #
          par_soil2_mm 1 par_x06 linear 1000.0 0.0 free
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
        # We don't want to use the parent class in this context
        pass
