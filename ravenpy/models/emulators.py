from collections import defaultdict, namedtuple
from dataclasses import dataclass
from pathlib import Path

import xarray as xr

from ravenpy.config.commands import BasinIndexCommand
from ravenpy.config.rvs import Config

from .base import Ostrich, Raven
from .rv import HRU, LU, RV, HRUState, Ost, Sub

__all__ = [
    "GR4JCN",
    "MOHYSE",
    "HMETS",
    "HBVEC",
    "BLENDED",
    "GR4JCN_OST",
    "MOHYSE_OST",
    "HMETS_OST",
    "HBVEC_OST",
    "BLENDED_OST",
    "get_model",
    "Routing",
]


class GR4JCN(Raven):
    """GR4J + Cemaneige global hydrological model

    References
    ----------
    Perrin, C., C. Michel and V. Andréassian (2003). Improvement of a parsimonious model for streamflow simulation.
    Journal of Hydrology, 279(1-4), 275-289. doi: 10.1016/S0022-1694(03)00225-7.

    Valéry, Audrey, Vazken Andréassian, and Charles Perrin. 2014. “’As Simple as Possible but Not Simpler’: What Is
    Useful in a Temperature-Based Snow-Accounting Routine? Part 2 - Sensitivity Analysis of the Cemaneige Snow
    Accounting Routine on 380 Catchments.” Journal of Hydrology, no. 517(0): 1176–87,
    doi: 10.1016/j.jhydrol.2014.04.058.
    """

    identifier = "gr4jcn"

    @dataclass
    class Params:
        GR4J_X1: float = None
        GR4J_X2: float = None
        GR4J_X3: float = None
        GR4J_X4: float = None
        CEMANEIGE_X1: float = None
        CEMANEIGE_X2: float = None

    @dataclass
    class DerivedParams:
        one_minus_CEMANEIGE_X2: float = None
        GR4J_X1_hlf: float = None

    @dataclass
    class LandHRU(HRU):
        land_use_class: str = "LU_ALL"
        veg_class: str = "VEG_ALL"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        _hru_type: str = "land"

    @dataclass
    class LakeHRU(HRU):
        land_use_class: str = "LU_WATER"
        veg_class: str = "VEG_WATER"
        soil_profile: str = "LAKE"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        _hru_type: str = "lake"

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self.config = Config(
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
            params=GR4JCN.Params(),
            derived_params=GR4JCN.DerivedParams(),
        )

        self.config.rvp.tmpl = """
        # -Global snow parameters-------------------------------------
        :RainSnowTransition 0 1.0
        :AirSnowCoeff       {derived_params.one_minus_CEMANEIGE_X2}  # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6
        :AvgAnnualSnow      {params.CEMANEIGE_X1}            # [mm]  =       CEMANEIGE_X1 =       x5

        # -Orographic Corrections-------------------------------------
        :PrecipitationLapseRate 0.0004
        :AdiabaticLapseRate 0.0065

        # - Soil classes ---------------------------------------------
        :SoilClasses
          :Attributes
          :Units
           SOIL_PROD
           SOIL_ROUT
           SOIL_TEMP
           SOIL_GW
           AQUIFER
        :EndSoilClasses
        :SoilParameterList
         :Parameters, POROSITY ,  GR4J_X3, GR4J_X2
         :Units     ,     none ,       mm,    mm/d
           [DEFAULT],      1.0 ,     {params.GR4J_X3},     {params.GR4J_X2}
        :EndSoilParameterList

        # ----Soil Profiles--------------------------------------------
        #     name,#horizons,(soiltype,thickness)x(#horizons)
        #     GR4J_X1 is thickness of first layer (SOIL_PROD), here {params.GR4J_X1}
        :SoilProfiles
          DEFAULT_P,      4, SOIL_PROD   ,   {params.GR4J_X1}, SOIL_ROUT  ,   0.300, SOIL_TEMP  ,   1.000, SOIL_GW  ,   1.000,
          LAKE,           0
        :EndSoilProfiles

        # ----Vegetation Classes---------------------------------------
        :VegetationClasses
           :Attributes,       MAX_HT,       MAX_LAI,      MAX_LEAF_COND
                :Units,            m,          none,           mm_per_s
               VEG_ALL,           0.0,          0.0,                0.0
               VEG_WATER,         0.0,          0.0,                0.0
        :EndVegetationClasses

        # --Land Use Classes------------------------------------------
        :LandUseClasses
          :Attributes, IMPERM, FOREST_COV
          :Units     ,   frac,       frac
               LU_ALL,    0.0,        0.0
               LU_WATER,  0.0,        0.0
        :EndLandUseClasses
        :LandUseParameterList
         :Parameters, GR4J_X4, MELT_FACTOR
         :Units     ,       d,      mm/d/C
           [DEFAULT],    {params.GR4J_X4},        7.73
        :EndLandUseParameterList

        {avg_annual_runoff}

        # List of channel profiles
        {channel_profiles}
        """

        self.config.rvi.tmpl = """
        :Calendar              {calendar}
        :RunName               {run_name}-{run_index}
        :StartDate             {start_date}
        :EndDate               {end_date}
        :TimeStep              {time_step}
        :Method                ORDERED_SERIES

        :SoilModel             SOIL_MULTILAYER  4
        {routing}
        :CatchmentRoute        ROUTE_DUMP
        :Evaporation           {evaporation}  # PET_OUDIN
        :RainSnowFraction      {rain_snow_fraction}  # RAINSNOW_DINGMAN
        :PotentialMeltMethod   POTMELT_DEGREE_DAY
        :OroTempCorrect        OROCORR_SIMPLELAPSE
        :OroPrecipCorrect      OROCORR_SIMPLELAPSE

        #------------------------------------------------------------------------
        # Soil Layer Alias Definitions
        #
        :Alias PRODUCT_STORE      SOIL[0]
        :Alias ROUTING_STORE      SOIL[1]
        :Alias TEMP_STORE         SOIL[2]
        :Alias GW_STORE           SOIL[3]

        #------------------------------------------------------------------------
        # Hydrologic process order for GR4J Emulation
        #
        :HydrologicProcesses
         :Precipitation            PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE
         :SnowTempEvolve           SNOTEMP_NEWTONS    SNOW_TEMP
         :SnowBalance              SNOBAL_CEMA_NIEGE  SNOW            PONDED_WATER
         :OpenWaterEvaporation     OPEN_WATER_EVAP    PONDED_WATER    ATMOSPHERE     			 # Pn
         :Infiltration             INF_GR4J           PONDED_WATER    MULTIPLE       			 # Ps-
         :SoilEvaporation          SOILEVAP_GR4J      PRODUCT_STORE   ATMOSPHERE     			 # Es
         :Percolation              PERC_GR4J          PRODUCT_STORE   TEMP_STORE     			 # Perc
         :Flush                    RAVEN_DEFAULT      SURFACE_WATER   TEMP_STORE     			 # Pn-Ps
         :Split                    RAVEN_DEFAULT      TEMP_STORE      CONVOLUTION[0] CONVOLUTION[1] 0.9  # Split Pr
         :Convolve                 CONVOL_GR4J_1      CONVOLUTION[0]  ROUTING_STORE  			 # Q9
         :Convolve                 CONVOL_GR4J_2      CONVOLUTION[1]  TEMP_STORE     			 # Q1
         :Percolation              PERC_GR4JEXCH      ROUTING_STORE   GW_STORE       			 # F(x1)
         :Percolation              PERC_GR4JEXCH2     TEMP_STORE      GW_STORE       			 # F(x1)
         :Flush                    RAVEN_DEFAULT      TEMP_STORE      SURFACE_WATER  			 # Qd
         :Baseflow                 BASE_GR4J          ROUTING_STORE   SURFACE_WATER  			 # Qr
        :EndHydrologicProcesses
        #------------------------------------------------------------------------

        #---------------------------------------------------------
        # Output Options
        #
        :WriteForcingFunctions
        :EvaluationMetrics {evaluation_metrics}
        :WriteNetcdfFormat  yes
        #:NoisyMode
        :SilentMode
        :PavicsMode
        {suppress_output}

        :NetCDFAttribute title Simulated river discharge
        :NetCDFAttribute history Created on {now} by Raven
        :NetCDFAttribute references  Craig, J.R., and the Raven Development Team, Raven user's and developer's manual (Version 2.8), URL: http://raven.uwaterloo.ca/ (2018).
        :NetCDFAttribute comment Raven Hydrological Framework version {raven_version}

        :NetCDFAttribute model_id gr4jcn

        :NetCDFAttribute time_frequency day
        :NetCDFAttribute time_coverage_start {start_date}
        :NetCDFAttribute time_coverage_end {end_date}
        """

        self.config.rvi.rain_snow_fraction = "RAINSNOW_DINGMAN"
        self.config.rvi.evaporation = "PET_OUDIN"

        # Initialize the stores to 1/2 full. Declare the parameters that can be user-modified
        self.config.rvc.soil0 = None
        self.config.rvc.soil1 = 15

    def derived_parameters(self):
        self.config.rvp.derived_params.GR4J_X1_hlf = (
            self.config.rvp.params.GR4J_X1 * 1000.0 / 2.0
        )
        self.config.rvp.derived_params.one_minus_CEMANEIGE_X2 = (
            1.0 - self.config.rvp.params.CEMANEIGE_X2
        )

        # Default initial conditions if none are given
        if not self.config.rvc.hru_states:
            soil0 = (
                self.config.rvp.derived_params.GR4J_X1_hlf
                if self.config.rvc.soil0 is None
                else self.config.rvc.soil0
            )
            soil1 = self.config.rvc.soil1

        # subbassin_id -> has at least one LakeHRU
        sb_contains_lake = defaultdict(lambda: False)

        if not self.config.rvc.hru_states:
            # If self.rvc.hru_states is set, it means that we are using `resume()` and we don't
            # want to interfere
            for hru in self.config.rvh.hrus:
                if isinstance(hru, GR4JCN.LandHRU) or hru._hru_type == "land":
                    self.config.rvc.hru_states[hru.hru_id] = HRUState(
                        index=hru.hru_id, soil0=soil0, soil1=soil1
                    )
                elif isinstance(hru, GR4JCN.LakeHRU) or hru._hru_type == "lake":
                    self.config.rvc.hru_states[hru.hru_id] = HRUState(index=hru.hru_id)
                    sb_contains_lake[hru.subbasin_id] = True
                else:
                    raise Exception(
                        "Type of HRU must be either `GR4JCN.LandHRU` or `GR4JCN.LakeHRU` (or its `_hru_type` must be either 'land' or 'lake')"
                    )

        if not self.config.rvc.basin_states:
            # If self.rvc.basin_states is set, it means that we are using `resume()` and we don't
            # want to interfere
            for sb in self.config.rvh.subbasins:
                self.config.rvc.basin_states[sb.subbasin_id] = BasinIndexCommand(
                    index=sb.subbasin_id
                )

        self.config.rvh.lake_subbasins = tuple(
            [
                sb.subbasin_id
                for sb in self.config.rvh.subbasins
                if sb_contains_lake[sb.subbasin_id]
            ]
        )
        self.config.rvh.land_subbasins = tuple(
            [
                sb.subbasin_id
                for sb in self.config.rvh.subbasins
                if not sb_contains_lake[sb.subbasin_id]
            ]
        )


class GR4JCN_OST(Ostrich, GR4JCN):
    _p = Path(__file__).parent / "ostrich-gr4j-cemaneige"
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob("*.t??"))

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvi.suppress_output = True
        self.txt = Ost(
            algorithm="DDS",
            max_iterations=50,
            lowerBounds=GR4JCN.params(None, None, None, None, None, None),
            upperBounds=GR4JCN.params(None, None, None, None, None, None),
        )

    def derived_parameters(self):
        """Derived parameters are computed by Ostrich."""
        pass


class MOHYSE(Raven):
    identifier = "mohyse"

    @dataclass
    class Params:
        par_x01: float = None
        par_x02: float = None
        par_x03: float = None
        par_x04: float = None
        par_x05: float = None
        par_x06: float = None
        par_x07: float = None
        par_x08: float = None
        par_x09: float = None
        par_x10: float = None

    @dataclass
    class DerivedParams:
        par_rezi_x10: float = None

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self.config = Config(
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
            params=MOHYSE.Params(),
            derived_params=MOHYSE.DerivedParams(),
        )

        self.config.rvp.tmpl = """
        #-----------------------------------------------------------------
        # Soil Classes
        #-----------------------------------------------------------------
        :SoilClasses
          :Attributes,
          :Units,
          TOPSOIL
          GWSOIL
        :EndSoilClasses

        #-----------------------------------------------------------------
        # Land Use Classes
        #-----------------------------------------------------------------
        :LandUseClasses,
          :Attributes,        IMPERM,    FOREST_COV,
          :Units,             frac,      frac,
          LU_ALL,             0.0,       1.0
        :EndLandUseClasses

        #-----------------------------------------------------------------
        # Vegetation Classes
        #-----------------------------------------------------------------
        :VegetationClasses,
          :Attributes,        MAX_HT,       MAX_LAI,    MAX_LEAF_COND,
          :Units,             m,            none,       mm_per_s,
         VEG_ALL,             0.0,          0.0,        0.0
        :EndVegetationClasses

        #-----------------------------------------------------------------
        # Soil Profiles
        #-----------------------------------------------------------------
        :SoilProfiles
                 LAKE, 0
                 ROCK, 0
               # DEFAULT_P,      2, TOPSOIL, MOHYSE_PARA_5, GWSOIL, 10.0
                 DEFAULT_P,      2, TOPSOIL,     {params.par_x05}, GWSOIL, 10.0
        :EndSoilProfiles

        #-----------------------------------------------------------------
        # Global Parameters
        #-----------------------------------------------------------------
        #:GlobalParameter      RAINSNOW_TEMP              -2.0
        :GlobalParameter       TOC_MULTIPLIER              1.0
        # :GlobalParameter     MOHYSE_PET_COEFF  MOHYSE_PARA_1
        :GlobalParameter       MOHYSE_PET_COEFF      {params.par_x01}

        #-----------------------------------------------------------------
        # Soil Parameters
        #-----------------------------------------------------------------
        :SoilParameterList
          :Parameters,        POROSITY,  PET_CORRECTION,        HBV_BETA,  BASEFLOW_COEFF,      PERC_COEFF,
               :Units,               -,               -,               -,             1/d,             1/d, # (units not generated by .rvp template)
            # TOPSOIL,            1.0 ,             1.0,             1.0,   MOHYSE_PARA_7,   MOHYSE_PARA_6,
            #  GWSOIL,            1.0 ,             1.0,             1.0,   MOHYSE_PARA_8,             0.0,
              TOPSOIL,            1.0 ,             1.0,             1.0,       {params.par_x07},       {params.par_x06},
               GWSOIL,            1.0 ,             1.0,             1.0,       {params.par_x08},             0.0,
        :EndSoilParameterList

        #-----------------------------------------------------------------
        # Land Use Parameters
        #-----------------------------------------------------------------
        :LandUseParameterList
          :Parameters,     MELT_FACTOR,       AET_COEFF, FOREST_SPARSENESS, DD_MELT_TEMP,
               :Units,          mm/d/K,            mm/d,                 -,         degC,
          # [DEFAULT],   MOHYSE_PARA_3,   MOHYSE_PARA_2,               0.0,MOHYSE_PARA_4,
            [DEFAULT],       {params.par_x03},       {params.par_x02},               0.0,    {params.par_x04},
        :EndLandUseParameterList

        #-----------------------------------------------------------------
        # Vegetation Parameters
        #-----------------------------------------------------------------
        :VegetationParameterList
          :Parameters,    SAI_HT_RATIO,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,
               :Units,               -,               -,               -,
            [DEFAULT],             0.0,             0.0,             0.0,
        :EndVegetationParameterList
        """

        self.config.rvi.tmpl = """
        :Calendar              {calendar}
        :RunName               {run_name}-{run_index}
        :StartDate             {start_date}
        :EndDate               {end_date}
        :TimeStep              {time_step}
        :Method                ORDERED_SERIES

        :SoilModel             SOIL_TWO_LAYER
        :PotentialMeltMethod   POTMELT_DEGREE_DAY
        :Routing               ROUTE_NONE
        :CatchmentRoute        ROUTE_GAMMA_CONVOLUTION
        :Evaporation           {evaporation}  # PET_MOHYSE
        :DirectEvaporation
        :RainSnowFraction      {rain_snow_fraction}

        :HydrologicProcesses
             :SoilEvaporation  SOILEVAP_LINEAR    SOIL[0]            ATMOSPHERE
             :SnowBalance      SNOBAL_SIMPLE_MELT SNOW PONDED_WATER
             :Precipitation    RAVEN_DEFAULT      ATMOS_PRECIP       MULTIPLE
             :Infiltration     INF_HBV            PONDED_WATER       SOIL[0]
             :Baseflow         BASE_LINEAR        SOIL[0]            SURFACE_WATER
             :Percolation      PERC_LINEAR        SOIL[0]            SOIL[1]
             :Baseflow         BASE_LINEAR        SOIL[1]            SURFACE_WATER
        :EndHydrologicProcesses

        #:CreateRVPTemplate

        # :Alias MOHYSE_PARA_1      1.5589    # :GlobalParameter         MOHYSE_PET_COEFF
        # :Alias MOHYSE_PARA_2	    0.9991    # LandUseParameterList --> AET_COEFF
        # :Alias MOHYSE_PARA_3	    2.1511    # LandUseParameterList --> MELT_FACTOR
        # :Alias MOHYSE_PARA_4	   -1.6101    # LandUseParameterList --> DD_MELT_TEMP
        # :Alias MOHYSE_PARA_5	    0.5000    # SoilProfiles         --> thickness of TOPSOIL (in mm????? must be m!!!)
        # :Alias MOHYSE_PARA_6	    0.1050    # SoilParameterList    --> PERC_COEFF (TOPSOIL)
        # :Alias MOHYSE_PARA_7	    0.0533    # SoilParameterList    --> BASEFLOW_COEFF (TOPSOIL)
        # :Alias MOHYSE_PARA_8	    0.0132    # SoilParameterList    --> BASEFLOW_COEFF (GWSOIL)
        # :Alias MOHYSE_PARA_9	    1.0474    # :SubBasinProperties  --> GAMMA_SHAPE
        # :Alias MOHYSE_PARA_10	    7.9628    # :SubBasinProperties  --> TIME_CONC = MOHYSE_PARA_10 / 0.3 = 26.542666666

        #---------------------------------------------------------
        # Output Options
        #
        #:WriteForcingFunctions
        :EvaluationMetrics {evaluation_metrics}
        :WriteNetcdfFormat  yes
        #:NoisyMode
        :SilentMode
        :PavicsMode
        {suppress_output}

        :NetCDFAttribute title Simulated river discharge
        :NetCDFAttribute history Created on {now} by Raven
        :NetCDFAttribute references  Craig, J.R., and the Raven Development Team, Raven user's and developer's manual (Version 2.8), URL: http://raven.uwaterloo.ca/ (2018).
        :NetCDFAttribute comment Raven Hydrological Framework version {raven_version}

        :NetCDFAttribute model_id mohyse

        :NetCDFAttribute time_frequency day
        :NetCDFAttribute time_coverage_start {start_date}
        :NetCDFAttribute time_coverage_end {end_date}
        """

        self.config.rvh.tmpl = """
        {subbasins}

        {hrus}

        :SubBasinProperties
        #                  1.0 / MOHYSE_PARA_10,   MOHYSE_PARA_9
           :Parameters,             GAMMA_SCALE,     GAMMA_SHAPE,
           :Units,                          1/d,               -
                      1,         {par_rezi_x10},       {par_x09}
        :EndSubBasinProperties
        """

        self.config.rvi.rain_snow_fraction = "RAINSNOW_DATA"
        self.config.rvi.evaporation = "PET_MOHYSE"

        # This is not stricly necessary it seems
        self.config.rvc.hru_states[1] = HRUState()
        self.config.rvc.basin_states[1] = BasinIndexCommand()

    def derived_parameters(self):
        self.config.rvp.derived_params.par_rezi_x10 = (
            1.0 / self.config.rvp.params.par_x10
        )

        # These need to be injected in the RVH
        self.config.rvh.par_rezi_x10 = self.config.rvp.derived_params.par_rezi_x10
        self.config.rvh.par_x09 = self.config.rvp.params.par_x09


class MOHYSE_OST(Ostrich, MOHYSE):
    _p = Path(__file__).parent / "ostrich-mohyse"
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob("*.t??"))

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvi.suppress_output = True
        self.txt = Ost(
            algorithm="DDS",
            max_iterations=50,
            lowerBounds=MOHYSE.params(
                None, None, None, None, None, None, None, None, None, None
            ),
            upperBounds=MOHYSE.params(
                None, None, None, None, None, None, None, None, None, None
            ),
        )

    def derived_parameters(self):
        """  Derived parameters are computed by Ostrich.  """
        pass


class HMETS(Raven):
    identifier = "hmets"

    @dataclass
    class Params:
        GAMMA_SHAPE: float = None
        GAMMA_SCALE: float = None
        GAMMA_SHAPE2: float = None
        GAMMA_SCALE2: float = None
        MIN_MELT_FACTOR: float = None
        MAX_MELT_FACTOR: float = None
        DD_MELT_TEMP: float = None
        DD_AGGRADATION: float = None
        SNOW_SWI_MIN: float = None
        SNOW_SWI_MAX: float = None
        SWI_REDUCT_COEFF: float = None
        DD_REFREEZE_TEMP: float = None
        REFREEZE_FACTOR: float = None
        REFREEZE_EXP: float = None
        PET_CORRECTION: float = None
        HMETS_RUNOFF_COEFF: float = None
        PERC_COEFF: float = None
        BASEFLOW_COEFF_1: float = None
        BASEFLOW_COEFF_2: float = None
        TOPSOIL: float = None
        PHREATIC: float = None

    @dataclass
    class DerivedParams:
        TOPSOIL_m: float = None
        PHREATIC_m: float = None
        SUM_MELT_FACTOR: float = None
        SUM_SNOW_SWI: float = None
        TOPSOIL_hlf: float = None
        PHREATIC_hlf: float = None

    @dataclass
    class ForestHRU(HRU):
        land_use_class: str = "FOREST"
        veg_class: str = "FOREST"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        # _hru_type: str = "land"

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self.config = Config(
            hrus=(HMETS.ForestHRU(),),
            subbasins=(
                Sub(
                    subbasin_id=1,
                    name="sub_001",
                    downstream_id=-1,
                    profile="None",
                    gauged=True,
                ),
            ),
            params=HMETS.Params(),
            derived_params=HMETS.DerivedParams(),
        )

        self.config.rvp.tmpl = """
        #-----------------------------------------------------------------
        # Soil Classes
        #-----------------------------------------------------------------
        :SoilClasses
          :Attributes,
          :Units,
          TOPSOIL,
          PHREATIC,
        :EndSoilClasses

        #-----------------------------------------------------------------
        # Land Use Classes
        #-----------------------------------------------------------------
        :LandUseClasses,
          :Attributes,        IMPERM,    FOREST_COV,
               :Units,          frac,          frac,
               FOREST,           0.0,           1.0,
        :EndLandUseClasses

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
          DEFAULT_P, 2, TOPSOIL,  {derived_params.TOPSOIL_m}, PHREATIC, {derived_params.PHREATIC_m},
        # DEFAULT_P, 2, TOPSOIL,   x(20)/1000, PHREATIC,   x(21)/1000,
        :EndSoilProfiles

        #-----------------------------------------------------------------
        # Global Parameters
        #-----------------------------------------------------------------
        :GlobalParameter         SNOW_SWI_MIN {params.SNOW_SWI_MIN}     # x(9)
        :GlobalParameter         SNOW_SWI_MAX {derived_params.SUM_SNOW_SWI}     # x(9)+x(10)
        :GlobalParameter     SWI_REDUCT_COEFF {params.SWI_REDUCT_COEFF} # x(11)
        :GlobalParameter             SNOW_SWI 0.05         # not sure why/if needed...

        #-----------------------------------------------------------------
        # Soil Parameters
        #-----------------------------------------------------------------
        :SoilParameterList
          :Parameters,        POROSITY,      PERC_COEFF,  PET_CORRECTION,  BASEFLOW_COEFF
               :Units,               -,             1/d,               -,             1/d
              TOPSOIL,             1.0,     {params.PERC_COEFF},{params.PET_CORRECTION},{params.BASEFLOW_COEFF_1}
             PHREATIC,             1.0,             0.0,             0.0, {params.BASEFLOW_COEFF_2}
         #    TOPSOIL,             1.0,           x(17),           x(15),           x(18)
         #   PHREATIC,             1.0,             0.0,             0.0,           x(19)
        :EndSoilParameterList

        #-----------------------------------------------------------------
        # Land Use Parameters
        #-----------------------------------------------------------------
        :LandUseParameterList
          :Parameters, MIN_MELT_FACTOR,  MAX_MELT_FACTOR,    DD_MELT_TEMP,  DD_AGGRADATION,  REFREEZE_FACTOR,    REFREEZE_EXP,  DD_REFREEZE_TEMP,  HMETS_RUNOFF_COEFF,
               :Units,          mm/d/C,          mm/d/C,               C,            1/mm,          mm/d/C,               -,                C,                  -,
            [DEFAULT],{params.MIN_MELT_FACTOR},{derived_params.SUM_MELT_FACTOR},  {params.DD_MELT_TEMP},{params.DD_AGGRADATION},{params.REFREEZE_FACTOR},  {params.REFREEZE_EXP},{params.DD_REFREEZE_TEMP},{params.HMETS_RUNOFF_COEFF},
        #                         x(5),       x(5)+x(6),            x(7),            x(8),           x(13),           x(14),            x(12),              x(16),
        :EndLandUseParameterList
        :LandUseParameterList
          :Parameters,   GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2,
               :Units,             -,             1/d,               -,             1/d,
            [DEFAULT],  {params.GAMMA_SHAPE},   {params.GAMMA_SCALE},  {params.GAMMA_SHAPE2},  {params.GAMMA_SCALE2},
            #                   x(1),            x(2),            x(3),            x(4),
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

        self.config.rvi.tmpl = """
        :Calendar              {calendar}
        :RunName               {run_name}-{run_index}
        :StartDate             {start_date}
        :EndDate               {end_date}
        :TimeStep              {time_step}
        :Method                ORDERED_SERIES

        :PotentialMeltMethod     POTMELT_HMETS
        :RainSnowFraction        {rain_snow_fraction}
        :Evaporation             {evaporation}  # PET_OUDIN
        :CatchmentRoute          ROUTE_DUMP
        :Routing                 ROUTE_NONE

        :SoilModel               SOIL_TWO_LAYER

        :Alias DELAYED_RUNOFF CONVOLUTION[1]

        :HydrologicProcesses
          :SnowBalance     SNOBAL_HMETS    MULTIPLE     MULTIPLE
          :Precipitation   RAVEN_DEFAULT   ATMOS_PRECIP MULTIPLE
          :Infiltration    INF_HMETS       PONDED_WATER MULTIPLE
            :Overflow      OVERFLOW_RAVEN  SOIL[0]      DELAYED_RUNOFF
          :Baseflow        BASE_LINEAR     SOIL[0]      SURFACE_WATER   # interflow, really
          :Percolation     PERC_LINEAR     SOIL[0]      SOIL[1]         # recharge
            :Overflow      OVERFLOW_RAVEN  SOIL[1]      DELAYED_RUNOFF
          :SoilEvaporation SOILEVAP_ALL    SOIL[0]      ATMOSPHERE      # AET
          :Convolve        CONVOL_GAMMA    CONVOLUTION[0] SURFACE_WATER #'surface runoff'
          :Convolve        CONVOL_GAMMA_2  DELAYED_RUNOFF SURFACE_WATER #'delayed runoff'
          :Baseflow        BASE_LINEAR     SOIL[1]      SURFACE_WATER
        :EndHydrologicProcesses

        #:CreateRVPTemplate

        #---------------------------------------------------------
        # Output Options
        #
        #:WriteForcingFunctions
        :EvaluationMetrics {evaluation_metrics}
        :WriteNetcdfFormat  yes
        #:NoisyMode
        :SilentMode
        :PavicsMode
        {suppress_output}

        :NetCDFAttribute title Simulated river discharge
        :NetCDFAttribute history Created on {now} by Raven
        :NetCDFAttribute references  Craig, J.R., and the Raven Development Team, Raven user's and developer's manual (Version 2.8), URL: http://raven.uwaterloo.ca/ (2018).
        :NetCDFAttribute comment Raven Hydrological Framework version {raven_version}

        :NetCDFAttribute model_id hmets

        :NetCDFAttribute time_frequency day
        :NetCDFAttribute time_coverage_start {start_date}
        :NetCDFAttribute time_coverage_end {end_date}
        """

        self.config.rvi.evaporation = "PET_OUDIN"
        self.config.rvi.rain_snow_fraction = "RAINSNOW_DATA"

        self.config.rvc.soil0 = None
        self.config.rvc.soil1 = None

    def derived_parameters(self):
        self.config.rvp.derived_params.TOPSOIL_hlf = (
            self.config.rvp.params.TOPSOIL * 0.5
        )
        self.config.rvp.derived_params.PHREATIC_hlf = (
            self.config.rvp.params.PHREATIC * 0.5
        )
        self.config.rvp.derived_params.TOPSOIL_m = (
            self.config.rvp.params.TOPSOIL / 1000.0
        )
        self.config.rvp.derived_params.PHREATIC_m = (
            self.config.rvp.params.PHREATIC / 1000.0
        )

        self.config.rvp.derived_params.SUM_MELT_FACTOR = (
            self.config.rvp.params.MAX_MELT_FACTOR
        )
        self.config.rvp.derived_params.SUM_SNOW_SWI = (
            self.config.rvp.params.SNOW_SWI_MAX
        )

        # Default initial conditions if none are given
        if not self.config.rvc.hru_states:
            soil0 = (
                self.config.rvp.derived_params.TOPSOIL_hlf
                if self.config.rvc.soil0 is None
                else self.config.rvc.soil0
            )
            soil1 = (
                self.config.rvp.derived_params.PHREATIC_hlf
                if self.config.rvc.soil1 is None
                else self.config.rvc.soil1
            )
            self.config.rvc.hru_states[1] = HRUState(soil0=soil0, soil1=soil1)


class HMETS_OST(Ostrich, HMETS):
    _p = Path(__file__).parent / "ostrich-hmets"
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob("*.t??"))

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvi.suppress_output = True
        self.txt = Ost(
            algorithm="DDS",
            max_iterations=50,
            lowerBounds=HMETS.params(
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ),
            upperBounds=HMETS.params(
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ),
        )

    def derived_parameters(self):
        """Derived parameters are computed by Ostrich."""
        pass

    def ost2raven(self, ops):
        """Return a list of parameter names calibrated by Ostrich that match Raven's parameters.

        Parameters
        ----------
        ops: dict
          Optimal parameter set returned by Ostrich.

        Returns
        -------
        HMETSParams named tuple
          Parameters expected by Raven.
        """
        names = ["par_x{:02}".format(i) for i in range(1, 22)]
        names[5] = "par_sum_x05_x06"
        names[9] = "par_sum_x09_x10"

        out = [ops[n] for n in names]
        out[19] *= 1000
        out[20] *= 1000
        return self.params(*out)


class HBVEC(Raven):
    identifier = "hbvec"

    @dataclass
    class Params:
        par_x01: float = None
        par_x02: float = None
        par_x03: float = None
        par_x04: float = None
        par_x05: float = None
        par_x06: float = None
        par_x07: float = None
        par_x08: float = None
        par_x09: float = None
        par_x10: float = None
        par_x11: float = None
        par_x12: float = None
        par_x13: float = None
        par_x14: float = None
        par_x15: float = None
        par_x16: float = None
        par_x17: float = None
        par_x18: float = None
        par_x19: float = None
        par_x20: float = None
        par_x21: float = None

    @dataclass
    class DerivedParams:
        one_plus_par_x15: float = None
        par_x11_half: float = None

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self.config = Config(
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
            params=HBVEC.Params(),
            derived_params=HBVEC.DerivedParams(),
        )

        self.config.rvp.tmpl = """
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
             FAST_RES,                _DEFAULT,      _DEFAULT,         0.0,    _DEFAULT,          _DEFAULT,    {params.par_x08},     {params.par_x09},    {derived_params.one_plus_par_x15}
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

        self.config.rvi.tmpl = """
        :Calendar              {calendar}
        :RunName               {run_name}-{run_index}
        :StartDate             {start_date}
        :EndDate               {end_date}
        :TimeStep              {time_step}
        :Method                ORDERED_SERIES

        #------------------------------------------------------------------------
        # Model options
        #
        :Method              	    ORDERED_SERIES
        #:Interpolation      	    INTERP_NEAREST_NEIGHBOR

        :Routing             	    ROUTE_NONE
        :CatchmentRoute      	    TRIANGULAR_UH

        :Evaporation         	    {evaporation}  # PET_FROM_MONTHLY
        :OW_Evaporation      	    {ow_evaporation}  # PET_FROM_MONTHLY
        :SWRadiationMethod   	    SW_RAD_DEFAULT
        :SWCloudCorrect      	    SW_CLOUD_CORR_NONE
        :SWCanopyCorrect     	    SW_CANOPY_CORR_NONE
        :LWRadiationMethod   	    LW_RAD_DEFAULT
        :RainSnowFraction    	    {rain_snow_fraction}  # RAINSNOW_HBV
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

        #---------------------------------------------------------
        # Output Options
        #
        #:WriteForcingFunctions
        :EvaluationMetrics {evaluation_metrics}
        :WriteNetcdfFormat  yes
        #:NoisyMode
        :SilentMode
        :PavicsMode
        {suppress_output}

        :NetCDFAttribute title Simulated river discharge
        :NetCDFAttribute history Created on {now} by Raven
        :NetCDFAttribute references  Craig, J.R., and the Raven Development Team, Raven user's and developer's manual (Version 2.8), URL: http://raven.uwaterloo.ca/ (2018).
        :NetCDFAttribute comment Raven Hydrological Framework version {raven_version}

        :NetCDFAttribute model_id hbvec

        :NetCDFAttribute time_frequency day
        :NetCDFAttribute time_coverage_start {start_date}
        :NetCDFAttribute time_coverage_end {end_date}
        """

        self.config.rvh.tmpl = """
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

        self.config.rvi.evaporation = "PET_FROMMONTHLY"
        self.config.rvi.ow_evaporation = "PET_FROMMONTHLY"
        self.config.rvi.rain_snow_fraction = "RAINSNOW_HBV"

        self.config.rvc.soil2 = 0.50657

    def derived_parameters(self):
        self.config.rvp.derived_params.one_plus_par_x15 = (
            self.config.rvp.params.par_x15 + 1.0
        )
        self.config.rvp.derived_params.par_x11_half = (
            self.config.rvp.params.par_x11 / 2.0
        )

        # These need to be injected in the RVH
        self.config.rvh.par_x11 = self.config.rvp.params.par_x11
        self.config.rvh.par_x11_half = self.config.rvp.derived_params.par_x11_half

        self.config.rvt.rain_correction = self.config.rvp.params.par_x20
        self.config.rvt.snow_correction = self.config.rvp.params.par_x21

        self._monthly_average()

        # Default initial conditions if none are given
        if not self.config.rvc.hru_states:
            self.config.rvc.hru_states[1] = HRUState(soil2=self.config.rvc.soil2)
        if not self.config.rvc.basin_states:
            self.config.rvc.basin_states[1] = BasinIndexCommand()

    # TODO: Support index specification and unit changes.
    def _monthly_average(self):

        if (
            self.config.rvi.evaporation == "PET_FROMMONTHLY"
            or self.config.rvi.ow_evaporation == "PET_FROMMONTHLY"
        ):
            # If this fails, it's likely the input data is missing some necessary variables (e.g. evap).
            tas_cmd = self.config.rvt._var_cmds["tas"]
            tasmin_cmd = self.config.rvt._var_cmds["tasmin"]
            tasmax_cmd = self.config.rvt._var_cmds["tasmax"]
            evspsbl_cmd = self.config.rvt._var_cmds["evspsbl"]

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
    _p = Path(__file__).parent / "ostrich-hbv-ec"
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob("*.t??"))

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvi.suppress_output = True
        self.low = HBVEC.params
        self.high = HBVEC.params
        self.txt = Ost(
            algorithm="DDS",
            max_iterations=50,
            lowerBounds=self.low(
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ),
            upperBounds=self.high(
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ),
        )

    # TODO: Support index specification and unit changes.
    def derived_parameters(self):
        self.rvt.raincorrection = "par_x20"
        self.rvt.snowcorrection = "par_x21"
        self._monthly_average()


class BLENDED(Raven):
    identifier = "blended"

    @dataclass
    class Params:
        par_x01: float = None
        par_x02: float = None
        par_x03: float = None
        par_x04: float = None
        par_x05: float = None
        par_x06: float = None
        par_x07: float = None
        par_x08: float = None
        par_x09: float = None
        par_x10: float = None
        par_x11: float = None
        par_x12: float = None
        par_x13: float = None
        par_x14: float = None
        par_x15: float = None
        par_x16: float = None
        par_x17: float = None
        par_x18: float = None
        par_x19: float = None
        par_x20: float = None
        par_x21: float = None
        par_x22: float = None
        par_x23: float = None
        par_x24: float = None
        par_x25: float = None
        par_x26: float = None
        par_x27: float = None
        par_x28: float = None
        par_x29: float = None
        par_x30: float = None
        par_x31: float = None
        par_x32: float = None
        par_x33: float = None
        par_x34: float = None
        par_x35: float = None
        par_r01: float = None
        par_r02: float = None
        par_r03: float = None
        par_r04: float = None
        par_r05: float = None
        par_r06: float = None
        par_r07: float = None
        par_r08: float = None

    @dataclass
    class DerivedParams:
        TOPSOIL_mm: float = None
        PHREATIC_mm: float = None
        TOPSOIL_hlf: float = None
        PHREATIC_hlf: float = None
        POW_X04: float = None
        POW_X11: float = None
        SUM_X09_X10: float = None
        SUM_X13_X14: float = None
        SUM_X24_X25: float = None

    @dataclass
    class ForestHRU(HRU):
        land_use_class: str = "FOREST"
        veg_class: str = "FOREST"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "DEFAULT_T"

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self.config = Config(
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
            params=BLENDED.Params(),
            derived_params=BLENDED.DerivedParams(),
            land_use_classes=(
                LU("FOREST", impermeable_frac=0.0, forest_coverage=0.02345),
            ),
            evaporation="PET_OUDIN",
            rain_snow_fraction="RAINSNOW_HBV",
        )

        self.config.rvp.tmpl = """
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
        :GlobalParameter         SNOW_SWI_MAX {derived_params.SUM_X13_X14}     #  x(13)+x(14)
        :GlobalParameter     SWI_REDUCT_COEFF {params.par_x15}  #  x(15)
        :GlobalParameter             SNOW_SWI {params.par_x19}  #  x(19)
        :GlobalParameter        RAINSNOW_TEMP {params.par_x31}  #  x(31)
        :GlobalParameter       RAINSNOW_DELTA {params.par_x32}  #  x(32)
        #:GlobalParameter      TOC_MULTIPLIER 1.0               #


        #-----------------------------------------------------------------
        # Soil Parameters
        #-----------------------------------------------------------------
        :SoilParameterList
          :Parameters,        POROSITY,           PERC_COEFF,  PET_CORRECTION,  BASEFLOW_COEFF,        B_EXP,     HBV_BETA, MAX_BASEFLOW_RATE,   BASEFLOW_N, FIELD_CAPACITY,     SAT_WILT,
               :Units,               -,                  1/d,               -,             1/d
              TOPSOIL,             1.0,         {params.par_x28},    {params.par_x08},       {derived_params.POW_X04}, {params.par_x02}, {params.par_x03}, 	 {params.par_x05}, {params.par_x06},  {derived_params.SUM_X09_X10}, {params.par_x09},
             PHREATIC,             1.0,         {params.par_x35},             0.0,       {derived_params.POW_X11},          0.0,          0.0, 	          0.0, {params.par_x12},            0.0,          0.0,
              DEEP_GW,             1.0,             0.0,                  0.0,             0.0,          0.0,          0.0, 	          0.0,          0.0,            0.0,          0.0,
         #    TOPSOIL,             1.0,           x(28),                x(08),           x(04),        x(02),        x(03), 	        x(05),        x(06),    x(09)+x(10),        x(09),
         #   PHREATIC,             1.0,           x(35),                  0.0,           x(11),          0.0,          0.0, 	          0.0,        x(12),            0.0,          0.0,
         #    DEEP_GW,             1.0,             0.0,                  0.0,             0.0,          0.0,          0.0, 	          0.0,          0.0,            0.0,          0.0,
        :EndSoilParameterList

        #-----------------------------------------------------------------
        # Land Use Parameters
        #-----------------------------------------------------------------
        :LandUseParameterList
          :Parameters, MIN_MELT_FACTOR,     MAX_MELT_FACTOR,    DD_MELT_TEMP,  DD_AGGRADATION, REFREEZE_FACTOR, REFREEZE_EXP, DD_REFREEZE_TEMP, HMETS_RUNOFF_COEFF,
               :Units,          mm/d/C,              mm/d/C,               C,            1/mm,          mm/d/C,            -,                C,                  -,
            [DEFAULT],    {params.par_x24},       {derived_params.SUM_X24_X25},    {params.par_x26},    {params.par_x27},    {params.par_x18}, {params.par_x17},     {params.par_x16},       {params.par_x01},
        #                        x(24),         x(24)+x(25),           x(26),           x(27),           x(18),        x(17),            x(16),              x(01),
        :EndLandUseParameterList
        :LandUseParameterList
          :Parameters,   GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2,    FOREST_SPARSENESS,
               :Units,             -,               -,               -,               -,                    -,
            [DEFAULT],  {params.par_x20},    {params.par_x21},    {params.par_x22},    {params.par_x23},                  0.0,
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

        self.config.rvi.tmpl = """
        :Calendar              {calendar}
        :RunName               {run_name}-{run_index}
        :StartDate             {start_date}
        :EndDate               {end_date}
        :TimeStep              {time_step}
        :Method                ORDERED_SERIES

        :PotentialMeltMethod     POTMELT_HMETS
        :RainSnowFraction        {rain_snow_fraction}  # RAINSNOW_HBV
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

        #:CreateRVPTemplate

        #---------------------------------------------------------
        # Output Options
        #
        #:WriteForcingFunctions
        :EvaluationMetrics {evaluation_metrics}
        :WriteNetcdfFormat  yes
        #:NoisyMode
        :SilentMode
        :PavicsMode
        {suppress_output}

        :NetCDFAttribute title Simulated river discharge
        :NetCDFAttribute history Created on {now} by Raven
        :NetCDFAttribute references  Craig, J.R., and the Raven Development Team, Raven user's and developer's manual (Version 2.8), URL: http://raven.uwaterloo.ca/ (2018).
        :NetCDFAttribute comment Raven Hydrological Framework version {raven_version}

        :NetCDFAttribute model_id hmets

        :NetCDFAttribute time_frequency day
        :NetCDFAttribute time_coverage_start {start_date}
        :NetCDFAttribute time_coverage_end {end_date}
        """

        self.config.rvi.params = self.config.rvp.params

        self.config.rvc.tmpl = """
        # initialize to 1/2 full
        # x(29)*1000/2
        :UniformInitialConditions SOIL[0] {TOPSOIL_hlf}
        # x(30)*1000/2
        :UniformInitialConditions SOIL[1] {PHREATIC_hlf}

        :HRUStateVariableTable (formerly :InitialConditionsTable)
           :Attributes SOIL[0] SOIL[1]
           :Units mm mm
           1 {TOPSOIL_hlf} {PHREATIC_hlf}
           # x(29)*1000/2       x(30)*1000/2
        :EndHRUStateVariableTable
        """

        self.config.rvc.soil0 = None
        self.config.rvc.soil1 = None
        self.config.rvc.basin_states[1] = BasinIndexCommand()

    def derived_parameters(self):
        self.config.rvp.derived_params.TOPSOIL_hlf = (
            self.config.rvp.params.par_x29 * 0.5 * 1000.0
        )
        self.config.rvp.derived_params.PHREATIC_hlf = (
            self.config.rvp.params.par_x30 * 0.5 * 1000.0
        )
        self.config.rvp.derived_params.TOPSOIL_mm = (
            self.config.rvp.params.par_x29 * 1000.0
        )
        self.config.rvp.derived_params.PHREATIC_mm = (
            self.config.rvp.params.par_x30 * 1000.0
        )
        self.config.rvp.derived_params.SUM_X09_X10 = (
            self.config.rvp.params.par_x10
        )  # + self.config.rvp.params.par_x09
        self.config.rvp.derived_params.SUM_X13_X14 = (
            self.config.rvp.params.par_x14
        )  # + self.config.rvp.params.par_x13
        self.config.rvp.derived_params.SUM_X24_X25 = (
            self.config.rvp.params.par_x25
        )  # + self.config.rvp.params.par_x24
        # 10.0**self.config.rvp.params.par_x04  #
        self.config.rvp.derived_params.POW_X04 = self.config.rvp.params.par_x04
        # 10.0**self.config.rvp.params.par_x11  #
        self.config.rvp.derived_params.POW_X11 = self.config.rvp.params.par_x11

        self.config.rvc.TOPSOIL_hlf = self.config.rvp.derived_params.TOPSOIL_hlf
        self.config.rvc.PHREATIC_hlf = self.config.rvp.derived_params.PHREATIC_hlf

        # Default initial conditions if none are given
        # if not self.config.rvc.hru_states:
        #     soil0 = (
        #         self.config.rvp.derived_params.TOPSOIL_hlf if self.config.rvc.soil0 is None else self.config.rvc.soil0
        #     )
        #     soil1 = (
        #         self.config.rvp.derived_params.PHREATIC_hlf if self.config.rvc.soil1 is None else self.config.rvc.soil1
        #     )
        #     self.config.rvc.hru_states[1] = HRUState(soil0=soil0, soil1=soil1)

        self.config.rvt.rain_correction = self.config.rvp.params.par_x33
        self.config.rvt.snow_correction = self.config.rvp.params.par_x34


class BLENDED_OST(Ostrich, BLENDED):
    _p = Path(__file__).parent / "ostrich-blended"
    templates = tuple(_p.glob("model/*.rv?")) + tuple(_p.glob("*.t??"))

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvi.suppress_output = True
        self.txt = Ost(
            algorithm="DDS",
            max_iterations=50,
            lowerBounds=BLENDED.params(
                None,  # par_x01
                None,  # par_x02
                None,  # par_x03
                None,  # 10**par_x04
                None,  # par_x05
                None,  # par_x06
                None,  # par_x07
                None,  # par_x08
                None,  # par_x09
                None,  # par_x09+par_x10
                None,  # 10**par_x11
                None,  # par_x12
                None,  # par_x13
                None,  # par_x13+par_x14
                None,  # par_x15
                None,  # par_x16
                None,  # par_x17
                None,  # par_x18
                None,  # par_x19
                None,  # par_x20
                None,  # par_x21
                None,  # par_x22
                None,  # par_x23
                None,  # par_x24
                None,  # par_x24+par_x25
                None,  # par_x26
                None,  # par_x27
                None,  # par_x28
                None,  # par_x29
                None,  # par_x30
                None,  # par_x31
                None,  # par_x32
                None,  # par_x33
                None,  # par_x34
                None,  # par_x35
                None,  # par_r01
                None,  # par_r02
                None,  # par_r03
                None,  # par_r04
                None,  # par_r05
                None,  # par_r06
                None,  # par_r07
                None,  # par_r08
            ),
            upperBounds=BLENDED.params(
                None,  # par_x01
                None,  # par_x02
                None,  # par_x03
                None,  # 10**par_x04
                None,  # par_x05
                None,  # par_x06
                None,  # par_x07
                None,  # par_x08
                None,  # par_x09
                None,  # par_x09+par_x10
                None,  # 10**par_x11
                None,  # par_x12
                None,  # par_x13
                None,  # par_x13+par_x14
                None,  # par_x15
                None,  # par_x16
                None,  # par_x17
                None,  # par_x18
                None,  # par_x19
                None,  # par_x20
                None,  # par_x21
                None,  # par_x22
                None,  # par_x23
                None,  # par_x24
                None,  # par_x24+par_x25
                None,  # par_x26
                None,  # par_x27
                None,  # par_x28
                None,  # par_x29
                None,  # par_x30
                None,  # par_x31
                None,  # par_x32
                None,  # par_x33
                None,  # par_x34
                None,  # par_x35
                None,  # par_r01
                None,  # par_r02
                None,  # par_r03
                None,  # par_r04
                None,  # par_r05
                None,  # par_r06
                None,  # par_r07
                None,  # par_r08
            ),
        )

    def derived_parameters(self):
        """Derived parameters are computed by Ostrich."""
        self.rvt.raincorrection = "par_x33"
        self.rvt.snowcorrection = "par_x34"

    def ost2raven(self, ops):
        """Return a list of parameter names calibrated by Ostrich that match Raven's parameters.

        Parameters
        ----------
        ops: dict
          Optimal parameter set returned by Ostrich.

        Returns
        -------
        BLENDEDParams named tuple
          Parameters expected by Raven.
        """
        names = ["par_x{:02}".format(i) for i in range(1, 36)] + [
            "par_r{:02}".format(i) for i in range(1, 9)
        ]
        names[3] = "pow_x04"
        names[9] = "sum_x09_x10"
        names[10] = "pow_x11"
        names[13] = "sum_x13_x14"
        names[24] = "sum_x24_x25"

        out = [ops[n] for n in names]
        return self.params(*out)


class Routing(Raven):
    """Routing model - no hydrological modeling"""

    identifier = "routing"

    # params = namedtuple("RoutingParams", ())

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        # Declare the parameters that can be user-modified.

        self.config = Config()
        self.config.rvp.tmpl = """
        {soil_classes}

        {soil_profiles}

        {vegetation_classes}

        {land_use_classes}

        {avg_annual_runoff}

        {channel_profiles}
        """

        self.config.rvi.tmpl = """
        :Calendar              {calendar}
        :RunName               {run_name}-{run_index}
        :StartDate             {start_date}
        :EndDate               {end_date}
        :TimeStep              {time_step}
        :Method                ORDERED_SERIES                # Numerical method used for simulation

        :CatchmentRoute        ROUTE_DUMP                    # Catchment routing method, used to convey water from the catchment tributaries and rivulets to the subbasin outlets. DEFAULT ROUTE_DUMP, which instantly ‘dumps’ all water in the subbasin stream reach.
        :Routing               ROUTE_DIFFUSIVE_WAVE          # Channel routing method which is used to transport water from upstream to downstream within the main subbasin channels. DEFAULT ROUTE_DIFFUSIVE_WAVE, which analytically solves the diffusive wave equation along the reach using a constant reference celerity.
        :PrecipIceptFract      PRECIP_ICEPT_NONE             # Estimation of the precipitation interception fraction. In this routing model, stream input(s) are "pretending" to be precipitation going into Raven, thus using DEFAULT PRECIP_ICEPT_NONE to indicate no interception processes are adopted.
        :PotentialMeltMethod   POTMELT_NONE                  # Estimation of the potential snow melt. In this routing model, snow melt processes are not relevant, thus using DEFAULT POTMELT_NONE method.
        :SoilModel             SOIL_ONE_LAYER                # In this routing model, use DEFAULT SOIL_ONE_LAYER to define single soil layer structure.

        :HydrologicProcesses
          :Precipitation     PRECIP_RAVEN             ATMOS_PRECIP     PONDED_WATER          # Moves stream input(s) from ATMOS_PRECIP to PONDED_WATER storage (waiting for runoff). Use DEFAULT PRECIP_RAVEN method.
          :Flush             RAVEN_DEFAULT            PONDED_WATER     SURFACE_WATER         # Moves water from PONDED_WATER to SURFACE_WATER (routed to outlet). Use DEFAULT RAVEN_DEFAULT method.
        :EndHydrologicProcesses


        # Output Options
        #
        #:WriteForcingFunctions
        # Defines the hydrograph performance metrics output by Raven. Either one or multiple is acceptable.
        :EvaluationMetrics {evaluation_metrics}
        :WriteNetcdfFormat  yes
        #:NoisyMode
        :SilentMode
        :PavicsMode
        {suppress_output}

        :NetCDFAttribute title Simulated river discharge
        :NetCDFAttribute history Created on {now} by Raven
        :NetCDFAttribute references  Craig, J.R., and the Raven Development Team, Raven user's and developer's manual (Version 2.8), URL: http://raven.uwaterloo.ca/ (2018).
        :NetCDFAttribute comment Raven Hydrological Framework version {raven_version}

        :NetCDFAttribute model_id routing

        :NetCDFAttribute time_frequency day
        :NetCDFAttribute time_coverage_start {start_date}
        :NetCDFAttribute time_coverage_end {end_date}
        """

    def derived_parameters(self):
        pass


def get_model(name):
    """Return the corresponding Raven emulated model instance.

    Parameters
    ----------
    name : str
      Model class name or model identifier.

    Returns
    -------
    Raven model instance
    """
    from ravenpy.models import emulators

    model_cls = getattr(emulators, name, None)

    if model_cls is None:
        for m in [GR4JCN, MOHYSE, HMETS, HBVEC, BLENDED]:
            if m.identifier == name:
                model_cls = m

    if model_cls is None:
        raise ValueError("Model {} is not recognized.".format(name))

    return model_cls


def used_storage_variables(fn):
    """Identify variables that are used by the model."""
    import xarray as xr

    ds = xr.open_dataset(fn)
    return [
        (key, da.isel(time=-1).values.tolist(), da.units)
        for key, da in ds.data_vars.items()
        if any(ds[key] != 0)
    ]
