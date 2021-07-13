from collections import defaultdict
from dataclasses import fields
from pathlib import Path
from typing import cast

import xarray as xr
from pydantic.dataclasses import dataclass

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


class HYPR(Raven):
    """
    The HYdrological model for Prairie Region (HYPR) is based on the
    conceptual Hydrologiska Byrans Vattenbalansavdelning (HBV)-light
    model. HBV is modified to work in the prairies by incorporating a
    conceptual lateral-flow component to represent the pothole storage
    complexities. HYPR can be used for prairie streamflow simulation.

    References
    ----------

    Mohamed I. Ahmed, Amin Elshorbagy, and Alain Pietroniro (2020).
    Toward Simple Modeling Practices in the Complex Canadian Prairie
    Watersheds. Journal of Hydrologic Engineering, 25(6), 04020024.
    doi: 10.1061/(ASCE)HE.1943-5584.0001922.
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
        land_use_class: str = "OPEN_1"
        veg_class: str = "FOREST"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "hypr")
        super().__init__(*args, **kwds)

        self.config.update(
            hrus=(HYPR.HRU(),),
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
           DEFAULT_P, 3, TOPSOIL, {params.par_x11}, FAST_RES,1e99, SLOW_RES, 1e99
        :EndSoilProfiles

        # --- Parameter Specification ----------------------------

        :GlobalParameter RAINSNOW_TEMP   {params.par_x02}
        #                                para_x02 = TT
        :GlobalParameter RAINSNOW_DELTA  0.0
        :GlobalParameter SNOW_SWI        {params.par_x03} # para_x03 = CWH
        :AdiabaticLapseRate              {params.par_x19} # para_x19 = TCALT
        :GlobalParameter PRECIP_LAPSE    {params.par_x20} # para_x20 = PCALT

        :SoilParameterList
          :Parameters, POROSITY,   FIELD_CAPACITY, SAT_WILT,         HBV_BETA, MAX_CAP_RISE_RATE, MAX_PERC_RATE,   BASEFLOW_COEFF, BASEFLOW_N,  BASEFLOW_COEFF2, STORAGE_THRESHOLD
          :Units     ,     none,             none,     none,             none,              mm/d,          mm/d,              1/d,       none,              1/d,                mm,
        #   [DEFAULT],      1.0,         para_x04,      0.0,         para_x16,              0.0,            0.0,              0.0,        0.0,              0.0,               0.0,
            [DEFAULT],      1.0, {params.par_x04},      0.0, {params.par_x16},              0.0,            0.0,              0.0,        0.0,              0.0,               0.0,
        #    FAST_RES, _DEFAULT,         _DEFAULT,      0.0,         _DEFAULT,         _DEFAULT,            0.0,         para_x06,        1.0,         para_x05,          para_x10,
             FAST_RES, _DEFAULT,         _DEFAULT,      0.0,         _DEFAULT,         _DEFAULT,            0.0, {params.par_x06},        1.0, {params.par_x05},   {params.par_x10},
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
             [DEFAULT], {params.par_x17}, {params.par_x18},    {params.par_x13}, {params.par_x12},              0.00, {params.par_x01},             0.0
        :EndLandUseParameterList

        :LandUseParameterList
          :Parameters,       PONDED_EXP,         PDMROF_B,             DEP_MAX, MAX_DEP_AREA_FRAC,
          :Units     ,             none,             none,                  mm,              none,
          #          ,            2/PWR,                B, SMAX=CMAX*(1/(B+1)),             MAXPA,
          # [DEFAULT],         para_x09,         para_x07,            para_x14,          para_x08,
          [DEFAULT],   {params.par_x09}, {params.par_x07},    {params.par_x14},  {params.par_x08},
        :EndLandUseParameterList
        """
        self.config.rvp.set_tmpl(rvp_tmpl)

        #########
        # R V I #
        #########

        rvi_tmpl = """
        :Routing                     {routing}
        :CatchmentRoute              TRIANGULAR_UH

        :Evaporation                 {evaporation}
        :OW_Evaporation              {ow_evaporation}
        :RainSnowFraction            {rain_snow_fraction}
        :SWRadiationMethod           SW_RAD_DEFAULT
        :SWCloudCorrect              SW_CLOUD_CORR_NONE
        :SWCanopyCorrect             SW_CANOPY_CORR_NONE
        :LWRadiationMethod           LW_RAD_DEFAULT
        :PotentialMeltMethod         POTMELT_HBV
        :CloudCoverMethod            CLOUDCOV_NONE
        :PrecipIceptFract            PRECIP_ICEPT_USER
        :MonthlyInterpolationMethod  MONTHINT_LINEAR_21

        :SoilModel               SOIL_MULTILAYER 3

        #------------------------------------------------------------------------
        # Soil Layer Alias Definitions
        #
        :Alias       FAST_RESERVOIR SOIL[1]
        :Alias       SLOW_RESERVOIR SOIL[2]
        :LakeStorage SLOW_RESERVOIR

        #------------------------------------------------------------------------
        # Hydrologic process order for HYPR Emulation
        #
        :HydrologicProcesses
          :SnowRefreeze      FREEZE_DEGREE_DAY  SNOW_LIQ        SNOW
          :Precipitation     PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE
          :CanopyEvaporation CANEVP_ALL         CANOPY          ATMOSPHERE
          :CanopySnowEvap    CANEVP_ALL         CANOPY_SNOW     ATMOSPHERE
          :SnowBalance       SNOBAL_SIMPLE_MELT SNOW            PONDED_WATER
          :Infiltration      INF_HBV            PONDED_WATER    MULTIPLE
          :Flush             RAVEN_DEFAULT      SURFACE_WATER   PONDED_WATER
          :Abstraction       ABST_PDMROF        PONDED_WATER    DEPRESSION
          :Flush             RAVEN_DEFAULT      SURFACE_WATER   FAST_RESERVOIR
          :SoilEvaporation   SOILEVAP_HYPR      MULTIPLE        ATMOSPHERE    #ET from both soils and depressions
          :Baseflow          BASE_LINEAR        FAST_RESERVOIR  SURFACE_WATER
          :Baseflow          BASE_THRESH_STOR   FAST_RESERVOIR  SURFACE_WATER
        :EndHydrologicProcesses
        """
        self.config.rvi.set_tmpl(rvi_tmpl)

        self.config.rvi.routing = RVI.RoutingOptions.NONE
        self.config.rvi.rain_snow_fraction = RVI.RainSnowFractionOptions.HBV
        self.config.rvi.evaporation = RVI.EvaporationOptions.PET_FROMMONTHLY
        self.config.rvi.ow_evaporation = RVI.EvaporationOptions.PET_FROMMONTHLY

        #########
        # R V H #
        #########

        rvh_tmpl = RVH.tmpl
        rvh_tmpl += """
        :SubBasinProperties
        #                                 x_15,
        #                               MAXBAS,
          :Parameters,               TIME_CONC,
          :Units     ,  none,                d,
                     1, {par_x15},
        :EndSubBasinProperties
        """
        self.config.rvh.set_tmpl(rvh_tmpl)

        #########
        # R V C #
        #########

        self.config.rvc.hru_states[1] = HRUState(soil0=10, soil1=10, depression=20)

    def derived_parameters(self):
        params = cast(HYPR.Params, self.config.rvp.params)
        self.config.rvt.snow_correction = params.par_x21
        self.config.rvh.set_extra_attributes(par_x15=params.par_x15)
        self._monthly_average()

    # TODO: Support index specification and unit changes.
    def _monthly_average(self):

        if (
            self.config.rvi.evaporation == "PET_FROMMONTHLY"
            or self.config.rvi.ow_evaporation == "PET_FROMMONTHLY"
        ):
            # If this fails, it's likely the input data is missing some necessary variables (e.g. evap).
            tas_cmd = cast(BaseDataCommand, self.config.rvt._var_cmds["tas"])
            tasmin_cmd = cast(BaseDataCommand, self.config.rvt._var_cmds["tasmin"])
            tasmax_cmd = cast(BaseDataCommand, self.config.rvt._var_cmds["tasmax"])
            evspsbl_cmd = cast(BaseDataCommand, self.config.rvt._var_cmds["evspsbl"])

            if tas_cmd:
                tas = xr.open_dataset(tas_cmd.file_name_nc)[tas_cmd.var_name_nc]
            else:
                tasmax = xr.open_dataset(tasmax_cmd.file_name_nc)[
                    tasmax_cmd.var_name_nc
                ]
                tasmin = xr.open_dataset(tasmin_cmd.file_name_nc)[
                    tasmin_cmd.var_name_nc
                ]
                tas = (tasmax + tasmin) / 2.0

            if evspsbl_cmd:
                evap = xr.open_dataset(evspsbl_cmd.file_name_nc)[
                    evspsbl_cmd.var_name_nc
                ]

            mat = tas.groupby("time.month").mean().values
            mae = evap.groupby("time.month").mean().values

            self.config.rvt.monthly_ave_evaporation = tuple(mae)
            self.config.rvt.monthly_ave_temperature = tuple(mat)


class HYPR_OST(Ostrich, HYPR):

    # Since the `par_x05` and `par_x06` values that Ostrich have found are the base-10 logarithms of the
    # corresponding Raven values, we perform the transformation here, so that Raven receives 10^par_x05
    # and 10^par_x06 for its own `Params.par_x05` and `Params.par_x06`, respectively.
    ostrich_to_raven_param_conversion = {
        "pow_x05": "par_x05",
        "pow_x06": "par_x06",
    }

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "hypr-ost")
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
        """
        self.config.rvp.set_tmpl(rvp_tmpl, is_ostrich=True)

        ####################
        # R V T (OST TMPL) #
        ####################

        self.config.rvt.set_tmpl(is_ostrich=True)

        self.config.rvt.update("snow_correction", "par_x21")

        ####################
        # R V H (OST TMPL) #
        ####################

        rvh_tmpl = RVH.tmpl
        rvh_tmpl += """
        :SubBasinProperties
        #                                 x_15,
        #                               MAXBAS,
          :Parameters,               TIME_CONC,
          :Units     ,  none,                d,
                     1,                par_x15,
        :EndSubBasinProperties
        """
        self.config.rvh.set_tmpl(rvh_tmpl, is_ostrich=True)

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
          hypr-ost.rvp.tpl;  hypr-ost.rvp
          hypr-ost.rvh.tpl;  hypr-ost.rvh
          hypr-ost.rvt.tpl;  hypr-ost.rvt
          #can be multiple (.rvh, .rvi)
        EndFilePairs

        #Parameter/DV Specification
        BeginParams
          #parameter  init.    low    high     tx_in  tx_ost  tx_out
          par_x01     random  -1.0     1.0     none   none    none  # TT     :: DD_MELT_TEMP
          par_x02     random  -3.0     3.0     none   none    none  # TT     :: RAINSNOW_TEMP
          par_x03     random   0.0     0.8     none   none    none  # CWH    :: SNOW_SWI
          par_x04     random   0.3     1.0     none   none    none  # LP     :: FIELD_CAPACITY
          par_x05     random  -1.3     0.3     none   none    none  # K0     :: BASEFLOW_COEFF[2] = [10^-1.3, 10^0.3] = [0.05, 2.0] (Julie added logarithmic sampling)
          par_x06     random  -2.0     0.0     none   none    none  # K1     :: BASEFLOW_COEFF[1] = [10^-2.0, 10^0.0] = [0.01, 1.0] (Julie added logarithmic sampling)
          par_x07     random   0.0    30.0     none   none    none  # B      :: PDMROF_B
          par_x08     random   0.1     0.8     none   none    none  # MAXPA  :: MAX_DEP_AREA_FRAC
          par_x09     random   0.4     2.0     none   none    none  # PWR    :: 2/PONDED_EXP
          par_x10     random   0.0   100.0     none   none    none  # UZL    :: Threshhold storage = STORAGE_THRESHOLD
          par_x11     random   0.0     0.5     none   none    none  # FC     :: max storage capacity  = THICKNESS[0]*POROSITY[0]
          #                                                                     calibrate only THICKNESS TOPSOIL  (range taken from xSSA paper)
          par_x12     random   0.0     5.0     none   none    none  # CFR    :: REFREEZE_FACTOR                   (range taken from xSSA paper)
          par_x13     random   0.0     1.0     none   none    none  # MRF    :: HBV_MELT_FOR_CORR                 (range from Raven manual)
          par_x14     random   0.0  1000.0     none   none    none  # CMAX   :: DEP_MAX = CMAX/(B+1)
          par_x15     random   0.01    6.0     none   none    none  # MAXBAS :: TIME_CONC
          par_x16     random   0.0     7.0     none   none    none  # BETA   :: HBV_BETA
          par_x17     random   0.0     8.0     none   none    none  # KMIN   :: MELT_FACTOR     (was in HYPR set to C0(paper)=KMIN(code)=MELT_FACTOR=MIN_MELT_FACTOR (here separated))
          par_x18     random   1.5     3.0     none   none    none  # KMIN   :: MIN_MELT_FACTOR (was in HYPR set to C0(paper)=KMIN(code)=MELT_FACTOR=MIN_MELT_FACTOR (here separated))
          par_x19     random   0.0     5.0     none   none    none  # TCALT  :: AdiabaticLapseRate           (Julie added according to HBV-EC setup)
          par_x20     random   0.0     5.0     none   none    none  # PCALT  :: PRECIP_LAPSE                 (Julie added according to HBV-EC setup)
          par_x21     random   0.8     1.2     none   none    none  # SCF    :: :SnowCorrection
        EndParams

        BeginTiedParams
          # ---------------------------------------------------------------
          # BASEFLOW_COEFF K0  = 10.0^x5
          #
          # pow_x05 = 10.0**(par_x5)
          # Xtied = c2 * base ** (c1 * X) + c0
          # --> c0   = 0.0
          # --> c1   = 1.0
          # --> c2   = 1.0
          # --> base = 10.0
          #
          pow_x05 1 par_x05 exp 10.0 1.0 1.0 0.0 free
          #
          # ---------------------------------------------------------------
          # BASEFLOW_COEFF K1 = 10.0^x6
          #
          # pow_x06 = 10.0**(par_x6)
          # Xtied = c2 * base ** (c1 * X) + c0
          # --> c0   = 0.0
          # --> c1   = 1.0
          # --> c2   = 1.0
          # --> base = 10.0
          #
          pow_x06 1 par_x06 exp 10.0 1.0 1.0 0.0 free
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
        self._monthly_average()
