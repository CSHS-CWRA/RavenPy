from collections import defaultdict
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
from ravenpy.models.base import Ostrich, Raven

from .gr4jcn import GR4JCN


class HBVEC(Raven):
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

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "hbvec")
        super().__init__(*args, **kwds)

        self.config.update(
            hrus=(GR4JCN.LandHRU(),),
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
        #------------------------------------------------------------------------
        # Global parameters
        #
        #                             HBV_PARA_13=TCALT
        :AdiabaticLapseRate                   {params.par_x13}
        #                                   HBV_PARA_01, CONSTANT,
        :RainSnowTransition                   {params.par_x01},      2.0
        #                                   HBV_PARA_04,
        :IrreducibleSnowSaturation            {params.par_x04}
        #                             HBV_PARA_12=PCALT
        :GlobalParameter PRECIP_LAPSE         {params.par_x12}

        #---------------------------------------------------------
        # Soil classes
        :SoilClasses
         :Attributes,
         :Units,
           TOPSOIL,      1.0,    0.0,       0
           SLOW_RES,     1.0,    0.0,       0
           FAST_RES,     1.0,    0.0,       0
        :EndSoilClasses

        :SoilParameterList
          :Parameters,                POROSITY,FIELD_CAPACITY,    SAT_WILT,    HBV_BETA, MAX_CAP_RISE_RATE,MAX_PERC_RATE,BASEFLOW_COEFF,            BASEFLOW_N
          :Units     ,                    none,          none,        none,        none,              mm/d,         mm/d,           1/d,                  none
          #                        HBV_PARA_05,   HBV_PARA_06, HBV_PARA_14, HBV_PARA_07,       HBV_PARA_16,     CONSTANT,      CONSTANT,              CONSTANT,
            [DEFAULT],               {params.par_x05},     {params.par_x06},   {params.par_x14},    {params.par_x07},        {params.par_x16},          0.0,           0.0,                   0.0
          #                                                       CONSTANT,                                  HBV_PARA_08,   HBV_PARA_09, 1+HBV_PARA_15=1+ALPHA,
             FAST_RES,                _DEFAULT,      _DEFAULT,         0.0,    _DEFAULT,          _DEFAULT,    {params.par_x08},     {params.par_x09},    {one_plus_par_x15}
          #                                                       CONSTANT,                                                 HBV_PARA_10,              CONSTANT,
             SLOW_RES,                _DEFAULT,      _DEFAULT,         0.0,    _DEFAULT,          _DEFAULT,     _DEFAULT,     {params.par_x10},                   1.0
        :EndSoilParameterList

        #---------------------------------------------------------
        # Soil profiles
        # name, layers, (soilClass, thickness) x layers
        #
        :SoilProfiles
        #                        HBV_PARA_17,           CONSTANT,           CONSTANT,
           DEFAULT_P, 3, TOPSOIL,  {params.par_x17}, FAST_RES,    100.0, SLOW_RES,    100.0
        :EndSoilProfiles

        #---------------------------------------------------------
        # Vegetation classes
        #
        :VegetationClasses
         :Attributes,   MAX_HT,  MAX_LAI, MAX_LEAF_COND
         :Units,             m,     none,      mm_per_s
           VEG_ALL,         25,      6.0,           5.3
        :EndVegetationClasses

        :VegetationParameterList
          :Parameters,  MAX_CAPACITY, MAX_SNOW_CAPACITY,  TFRAIN,  TFSNOW,
          :Units,                 mm,                mm,    frac,    frac,
          VEG_ALL,             10000,             10000,    0.88,    0.88,
        :EndVegetationParameterList

        #---------------------------------------------------------
        # LandUse classes
        #
        :LandUseClasses
         :Attributes,     IMPERM, FOREST_COV
         :Units,            frac,       frac
              LU_ALL,        0.0,          1
        :EndLandUseClasses

        :LandUseParameterList
          :Parameters,   MELT_FACTOR, MIN_MELT_FACTOR,   HBV_MELT_FOR_CORR, REFREEZE_FACTOR, HBV_MELT_ASP_CORR
          :Units     ,        mm/d/K,          mm/d/K,                none,          mm/d/K,              none
          #              HBV_PARA_02,        CONSTANT,         HBV_PARA_18,     HBV_PARA_03,          CONSTANT
            [DEFAULT],     {params.par_x02},             2.2,           {params.par_x18},       {params.par_x03},              0.48
        :EndLandUseParameterList

        :LandUseParameterList
         :Parameters, HBV_MELT_GLACIER_CORR,   HBV_GLACIER_KMIN, GLAC_STORAGE_COEFF, HBV_GLACIER_AG
         :Units     ,                  none,                1/d,                1/d,           1/mm
           #                       CONSTANT,           CONSTANT,        HBV_PARA_19,       CONSTANT,
           [DEFAULT],                  1.64,               0.05,          {params.par_x19},           0.05
        :EndLandUseParameterList
        """
        self.config.rvp.set_tmpl(rvp_tmpl)

        #########
        # R V I #
        #########

        rvi_tmpl = """
        :Routing             	    ROUTE_NONE
        :CatchmentRoute      	    TRIANGULAR_UH

        :Evaporation         	    {evaporation}  # PET_FROM_MONTHLY
        :OW_Evaporation      	    {ow_evaporation}  # PET_FROM_MONTHLY
        :SWRadiationMethod   	    SW_RAD_DEFAULT
        :SWCloudCorrect      	    SW_CLOUD_CORR_NONE
        :SWCanopyCorrect     	    SW_CANOPY_CORR_NONE
        :LWRadiationMethod   	    LW_RAD_DEFAULT
        :RainSnowFraction           {rain_snow_fraction}
        :PotentialMeltMethod 	    POTMELT_HBV
        :OroTempCorrect      	    OROCORR_HBV
        :OroPrecipCorrect    	    OROCORR_HBV
        :OroPETCorrect       	    OROCORR_HBV
        :CloudCoverMethod    	    CLOUDCOV_NONE
        :PrecipIceptFract    	    PRECIP_ICEPT_USER
        :MonthlyInterpolationMethod MONTHINT_LINEAR_21

        :SoilModel                  SOIL_MULTILAYER 3

        #------------------------------------------------------------------------
        # Soil Layer Alias Definitions
        #
        :Alias       FAST_RESERVOIR SOIL[1]
        :Alias       SLOW_RESERVOIR SOIL[2]
        :LakeStorage SLOW_RESERVOIR

        #------------------------------------------------------------------------
        # Hydrologic process order for HBV-EC Emulation
        #
        :HydrologicProcesses
          :SnowRefreeze      FREEZE_DEGREE_DAY  SNOW_LIQ        SNOW
          :Precipitation     PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE
          :CanopyEvaporation CANEVP_ALL         CANOPY          ATMOSPHERE
          :CanopySnowEvap    CANEVP_ALL         CANOPY_SNOW     ATMOSPHERE
          :SnowBalance       SNOBAL_SIMPLE_MELT SNOW            SNOW_LIQ
            :-->Overflow     RAVEN_DEFAULT      SNOW_LIQ        PONDED_WATER
          :Flush             RAVEN_DEFAULT      PONDED_WATER    GLACIER
            :-->Conditional HRU_TYPE IS GLACIER
          :GlacierMelt       GMELT_HBV          GLACIER_ICE     GLACIER
          :GlacierRelease    GRELEASE_HBV_EC    GLACIER         SURFACE_WATER
          :Infiltration      INF_HBV            PONDED_WATER    MULTIPLE
          :Flush             RAVEN_DEFAULT      SURFACE_WATER   FAST_RESERVOIR
            :-->Conditional HRU_TYPE IS_NOT GLACIER
          :SoilEvaporation   SOILEVAP_HBV       SOIL[0]         ATMOSPHERE
          :CapillaryRise     RISE_HBV           FAST_RESERVOIR 	SOIL[0]
          :LakeEvaporation   LAKE_EVAP_BASIC    SLOW_RESERVOIR  ATMOSPHERE
          :Percolation       PERC_CONSTANT      FAST_RESERVOIR 	SLOW_RESERVOIR
          :Baseflow          BASE_POWER_LAW     FAST_RESERVOIR  SURFACE_WATER
          :Baseflow          BASE_LINEAR        SLOW_RESERVOIR  SURFACE_WATER
        :EndHydrologicProcesses
        """
        self.config.rvi.set_tmpl(rvi_tmpl)

        #########
        # R V H #
        #########

        rvh_tmpl = """
        {subbasins}

        {hrus}

        :SubBasinProperties
        #                       HBV_PARA_11, DERIVED FROM HBV_PARA_11,
        #                            MAXBAS,                 MAXBAS/2,
           :Parameters,           TIME_CONC,             TIME_TO_PEAK
           :Units     ,                   d,                        d,
                     1,            {par_x11},          {par_x11_half},
        :EndSubBasinProperties
        """
        self.config.rvh.set_tmpl(rvh_tmpl)

        #########
        # R V I #
        #########

        self.config.rvi.evaporation = "PET_FROMMONTHLY"
        self.config.rvi.ow_evaporation = "PET_FROMMONTHLY"
        self.config.rvi.rain_snow_fraction = "RAINSNOW_HBV"

    def derived_parameters(self):
        params = cast(HBVEC.Params, self.config.rvp.params)
        self.config.rvp.set_extra_attributes(
            one_plus_par_x15=params.par_x15 + 1.0, par_x11_half=params.par_x11 / 2.0
        )

        self.config.rvh.set_extra_attributes(
            par_x11=params.par_x11, par_x11_half=params.par_x11 / 2.0
        )

        self.config.rvt.rain_correction = params.par_x20
        self.config.rvt.snow_correction = params.par_x21

        self._monthly_average()

        # Default initial conditions if none are given
        if not self.config.rvc.hru_states:
            soil2 = 0.50657
            self.config.rvc.hru_states[1] = HRUState(soil2=soil2)

        if not self.config.rvc.basin_states:
            self.config.rvc.basin_states[1] = BasinIndexCommand()

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


class HBVEC_OST(Ostrich, HBVEC):
    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "hbvec-ost")
        super().__init__(*args, **kwds)

        self.config.update(
            algorithm="DDS",
            max_iterations=50,
            run_name="run",
            run_index=0,
            suppress_output=True,
        )

        #########
        # R V C #
        #########

        rvc_tmpl = """
        :BasinInitialConditions
         :Attributes, ID,              Q
         :Units,      none,         m3/s
         #                  HBV_PARA_???
                      1,             1.0
        :EndBasinInitialConditions

        #---------------------------------------------------------
        # Initial Lower groundwater storage - for each HRU

        :InitialConditions SOIL[2]
             # HBV_PARA_???
             0.50657
        :EndInitialConditions
        """
        self.config.rvc.set_tmpl(rvc_tmpl)

        ####################
        # R V H (OST TMPL) #
        ####################

        rvh_tmpl = """
        :SubBasins
                :Attributes     NAME    DOWNSTREAM_ID   PROFILE   REACH_LENGTH    GAUGED
                :Units          none    none            none      km              none
                1,               hbv,   -1,             NONE,     _AUTO,          1
        :EndSubBasins

        :HRUs
                :Attributes     AREA    ELEVATION  LATITUDE    LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT
                :Units           km2            m       deg          deg       none            none       none           none             none            none   ratio      deg
                     1,       4250.6,       843.0,  54.4848,   -123.3659,         1,         LU_ALL,   VEG_ALL,     DEFAULT_P,          [NONE],         [NONE], [NONE],  [NONE]
        :EndHRUs

        :SubBasinProperties
        #                       HBV_PARA_11, DERIVED FROM HBV_PARA_11,
        #                            MAXBAS,                 MAXBAS/2,
           :Parameters,           TIME_CONC,             TIME_TO_PEAK
           :Units     ,                   d,                        d,
                     1,             par_x11,             par_half_x11,
        :EndSubBasinProperties
        """
        self.config.rvh.set_tmpl(rvh_tmpl, is_ostrich=True)

        ####################
        # R V P (OST TMPL) #
        ####################

        rvp_tmpl = """
        # tied parameters:
        # (it is important for OSTRICH to find every parameter place holder somewhere in this file)
        # (without this "para_x05" and "para_x15" wouldn't be detectable)
        #    para_1_+_x15       = 1.0 + par_x15
        #    para_half_x11      = 0.5 * par_x11

        #------------------------------------------------------------------------
        # Global parameters
        #
        #                             HBV_PARA_13=TCALT
        :AdiabaticLapseRate                     par_x13
        #                                   HBV_PARA_01, CONSTANT,
        :RainSnowTransition                     par_x01       2.0
        #                                   HBV_PARA_04,
        :IrreducibleSnowSaturation              par_x04
        #                             HBV_PARA_12=PCALT
        :GlobalParameter PRECIP_LAPSE           par_x12

        #---------------------------------------------------------
        # Soil classes
        :SoilClasses
         :Attributes,
         :Units,
           TOPSOIL,      1.0,    0.0,       0
           SLOW_RES,     1.0,    0.0,       0
           FAST_RES,     1.0,    0.0,       0
        :EndSoilClasses

        :SoilParameterList
          :Parameters,                POROSITY,FIELD_CAPACITY,    SAT_WILT,    HBV_BETA, MAX_CAP_RISE_RATE,MAX_PERC_RATE,BASEFLOW_COEFF,            BASEFLOW_N
          :Units     ,                    none,          none,        none,        none,              mm/d,         mm/d,           1/d,                  none
          #                        HBV_PARA_05,   HBV_PARA_06, HBV_PARA_14, HBV_PARA_07,       HBV_PARA_16,     CONSTANT,      CONSTANT,              CONSTANT,
            [DEFAULT],                 par_x05,       par_x06,     par_x14,     par_x07,           par_x16,          0.0,           0.0,                   0.0
          #                                                       CONSTANT,                                  HBV_PARA_08,   HBV_PARA_09, 1+HBV_PARA_15=1+ALPHA,
             FAST_RES,                _DEFAULT,      _DEFAULT,         0.0,    _DEFAULT,          _DEFAULT,      par_x08,       par_x09,           par_1_+_x15
          #                                                       CONSTANT,                                                 HBV_PARA_10,              CONSTANT,
             SLOW_RES,                _DEFAULT,      _DEFAULT,         0.0,    _DEFAULT,          _DEFAULT,     _DEFAULT,       par_x10,                   1.0
        :EndSoilParameterList

        #---------------------------------------------------------
        # Soil profiles
        # name, layers, (soilClass, thickness) x layers
        #
        :SoilProfiles
        #                        HBV_PARA_17,           CONSTANT,           CONSTANT,
           DEFAULT_P, 3, TOPSOIL,    par_x17, FAST_RES,    100.0, SLOW_RES,    100.0
        :EndSoilProfiles

        #---------------------------------------------------------
        # Vegetation classes
        #
        :VegetationClasses
         :Attributes,   MAX_HT,  MAX_LAI, MAX_LEAF_COND
         :Units,             m,     none,      mm_per_s
           VEG_ALL,         25,      6.0,           5.3
        :EndVegetationClasses

        :VegetationParameterList
          :Parameters,  MAX_CAPACITY, MAX_SNOW_CAPACITY,  TFRAIN,  TFSNOW,
          :Units,                 mm,                mm,    frac,    frac,
          VEG_ALL,             10000,             10000,    0.88,    0.88,
        :EndVegetationParameterList

        #---------------------------------------------------------
        # LandUse classes
        #
        :LandUseClasses
         :Attributes,     IMPERM, FOREST_COV
         :Units,            frac,       frac
              LU_ALL,        0.0,          1
        :EndLandUseClasses

        :LandUseParameterList
          :Parameters,   MELT_FACTOR, MIN_MELT_FACTOR,   HBV_MELT_FOR_CORR, REFREEZE_FACTOR, HBV_MELT_ASP_CORR
          :Units     ,        mm/d/K,          mm/d/K,                none,          mm/d/K,              none
          #              HBV_PARA_02,        CONSTANT,         HBV_PARA_18,     HBV_PARA_03,          CONSTANT
            [DEFAULT],       par_x02,             2.2,             par_x18,         par_x03,              0.48
        :EndLandUseParameterList

        :LandUseParameterList
         :Parameters, HBV_MELT_GLACIER_CORR,   HBV_GLACIER_KMIN, GLAC_STORAGE_COEFF, HBV_GLACIER_AG
         :Units     ,                  none,                1/d,                1/d,           1/mm
           #                       CONSTANT,           CONSTANT,        HBV_PARA_19,       CONSTANT,
           [DEFAULT],                  1.64,               0.05,            par_x19,           0.05
        :EndLandUseParameterList
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

        #########
        # O S T #
        #########

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
          {identifier}.rvp.tpl;  {identifier}.rvp
          {identifier}.rvh.tpl;  {identifier}.rvh
          {identifier}.rvt.tpl;  {identifier}.rvt
          #can be multiple (.rvh, .rvi)
        EndFilePairs

        #Parameter/DV Specification
        BeginParams
          #parameter      init.    low      high    tx_in  tx_ost  tx_out   # in HBV called
          par_x01         random   {lowerBounds.par_x01}  {upperBounds.par_x01}  none   none     none    TT
          par_x02         random   {lowerBounds.par_x02}  {upperBounds.par_x02}  none   none     none    CFMAX
          par_x03         random   {lowerBounds.par_x03}  {upperBounds.par_x03}  none   none     none    CFR
          par_x04         random   {lowerBounds.par_x04}  {upperBounds.par_x04}  none   none     none    CWH
          par_x05         random   {lowerBounds.par_x05}  {upperBounds.par_x05}  none   none     none    par_x5 = FC/par_x21
          par_x06         random   {lowerBounds.par_x06}  {upperBounds.par_x06}  none   none     none    LP
          par_x07         random   {lowerBounds.par_x07}  {upperBounds.par_x07}  none   none     none    BETA
          par_x08         random   {lowerBounds.par_x08}  {upperBounds.par_x08}  none   none     none    PERC
          par_x09         random   {lowerBounds.par_x09}  {upperBounds.par_x09}  none   none     none    K1
          par_x10         random   {lowerBounds.par_x10}  {upperBounds.par_x10}  none   none     none    K2
          par_x11         random   {lowerBounds.par_x11}  {upperBounds.par_x11}  none   none     none    MAXBAS
          par_x12         random   {lowerBounds.par_x12}  {upperBounds.par_x12}  none   none     none    PCALT
          par_x13         random   {lowerBounds.par_x13}  {upperBounds.par_x13}  none   none     none    TCALT
          par_x14         random   {lowerBounds.par_x14}  {upperBounds.par_x14}  none   none     none    saturation at the wilting point
          par_x15         random   {lowerBounds.par_x15}  {upperBounds.par_x15}  none   none     none    ALPHA
          par_x16         random   {lowerBounds.par_x16}  {upperBounds.par_x16}  none   none     none    maximum interflow rate for capillary rise
          par_x17         random   {lowerBounds.par_x17}  {upperBounds.par_x17}  none   none     none    thickness of top soil layer
          par_x18         random   {lowerBounds.par_x18}  {upperBounds.par_x18}  none   none     none    melt correction factor (forest)
          par_x19         random   {lowerBounds.par_x19}  {upperBounds.par_x19}  none   none     none    release from glacier as it melts
          par_x20         random   {lowerBounds.par_x20}  {upperBounds.par_x20}  none   none     none    RFCF
          par_x21         random   {lowerBounds.par_x21}  {upperBounds.par_x21}  none   none     none    SFCF
        EndParams

        BeginTiedParams
          # par_1_+_x15 = 1.0 + par_x15
          # Xtied = (c1 * X1) + c0
          # --> c0 = 1.0
          # --> c1 = 1.0
          #
          par_1_+_x15 1 par_x15 linear 1.00 1.00 free
          #
          # par_half_x11 = par_x11 * 0.5
          # Xtied = (c1 * X) + c0
          # --> c0 = 0.0
          # --> c1 = 0.5
          #
          par_half_x11 1 par_x11 linear 0.5 0.0 free
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

    # TODO: Support index specification and unit changes.
    def derived_parameters(self):
        self.config.rvt.rain_correction = "par_x20"
        self.config.rvt.snow_correction = "par_x21"
        self._monthly_average()
