from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ravenpy.config.commands import BasinIndexCommand
from ravenpy.models.base import Ostrich, Raven
from ravenpy.models.rv import HRU, LU, RV, HRUState, Sub

from .gr4jcn import GR4JCN


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
