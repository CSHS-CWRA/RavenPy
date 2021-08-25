from collections import defaultdict
from pathlib import Path
from typing import cast

from pydantic.dataclasses import dataclass

from ravenpy.config.commands import HRU, LU, BasinIndexCommand, HRUState, Sub
from ravenpy.config.rvs import RVI
from ravenpy.models.base import Ostrich, Raven


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

    @dataclass
    class Params:
        GR4J_X1: float
        GR4J_X2: float
        GR4J_X3: float
        GR4J_X4: float
        CEMANEIGE_X1: float
        CEMANEIGE_X2: float

    @dataclass
    class LandHRU(HRU):
        land_use_class: str = "LU_ALL"
        veg_class: str = "VEG_ALL"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        hru_type: str = "land"

    @dataclass
    class LakeHRU(HRU):
        land_use_class: str = "LU_WATER"
        veg_class: str = "VEG_WATER"
        soil_profile: str = "LAKE"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        hru_type: str = "lake"

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "gr4jcn")
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
        # -Global snow parameters-------------------------------------
        :RainSnowTransition 0 1.0
        :AirSnowCoeff       {one_minus_CEMANEIGE_X2}  # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6
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
        self.config.rvp.set_tmpl(rvp_tmpl)

        #########
        # R V I #
        #########

        rvi_tmpl = """
        :SoilModel             SOIL_MULTILAYER  4
        :Routing               {routing}
        :CatchmentRoute        ROUTE_DUMP
        :Evaporation           {evaporation}
        :RainSnowFraction      {rain_snow_fraction}

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
        """
        self.config.rvi.set_tmpl(rvi_tmpl)

        self.config.rvi.rain_snow_fraction = RVI.RainSnowFractionOptions.DINGMAN
        self.config.rvi.evaporation = "PET_OUDIN"

    def derived_parameters(self):

        params = cast(GR4JCN.Params, self.config.rvp.params)

        self.config.rvp.set_extra_attributes(
            one_minus_CEMANEIGE_X2=(1.0 - params.CEMANEIGE_X2)
        )

        # subbassin_id -> has at least one LakeHRU
        sb_contains_lake = defaultdict(lambda: False)

        if not self.config.rvc.hru_states:
            # If self.rvc.hru_states is set, it means that we are using `resume()` and we don't
            # want to interfere

            soil0 = params.GR4J_X1 * 1000.0 / 2.0
            soil1 = 15

            for hru in self.config.rvh.hrus:
                if isinstance(hru, GR4JCN.LandHRU) or hru.hru_type == "land":
                    self.config.rvc.hru_states[hru.hru_id] = HRUState(
                        index=hru.hru_id, soil0=soil0, soil1=soil1
                    )
                elif isinstance(hru, GR4JCN.LakeHRU) or hru.hru_type == "lake":
                    self.config.rvc.hru_states[hru.hru_id] = HRUState(index=hru.hru_id)
                    sb_contains_lake[hru.subbasin_id] = True
                else:
                    raise Exception(
                        "Type of HRU must be either `GR4JCN.LandHRU` or `GR4JCN.LakeHRU`"
                    )

        if not self.config.rvc.basin_states:
            # If self.rvc.basin_states is set, it means that we are using `resume()` and we don't
            # want to interfere
            for sb in self.config.rvh.subbasins:
                self.config.rvc.basin_states[sb.subbasin_id] = BasinIndexCommand(
                    index=sb.subbasin_id
                )


class GR4JCN_OST(Ostrich, GR4JCN):

    ostrich_to_raven_param_conversion = {
        "par_x1": "GR4J_X1",
        "par_x2": "GR4J_X2",
        "par_x3": "GR4J_X3",
        "par_x4": "GR4J_X4",
        "par_x5": "CEMANEIGE_X1",
        "par_x6": "CEMANEIGE_X2",
    }

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "gr4jcn-ost")
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
        # Tied parameters:
        # (it is important for OSTRICH to find every parameter place holder somewhere in this file)
        # (without this "para_x1" wouldn't be detectable)
        #    para_half_x1 = para_x1 / 2. = par_x1 / 2. [m] = par_half_x1 [mm]

        # initialize to 1/2 full
        # GR4J_X1 * 1000. / 2.0
        :UniformInitialConditions SOIL[0] par_half_x1

        # Fixed because SOIL_ROUT layer thickness is fixed to be 0.3m
        :UniformInitialConditions SOIL[1] 15.0

        :HRUStateVariableTable (formerly :InitialConditionsTable)
           :Attributes SOIL[0]       SOIL[1]
           :Units      mm            mm
           1           par_half_x1   15.0
        :EndHRUStateVariableTable
        """
        self.config.rvc.set_tmpl(rvc_tmpl, is_ostrich=True)

        ####################
        # R V P (OST TMPL) #
        ####################

        rvp_tmpl = """
        #########################################################################
        :FileType          rvp ASCII Raven 2.8.2
        :WrittenBy         Juliane Mai & James Craig
        :CreationDate      Sep 2018
        #
        # Emulation of GR4J simulation of Salmon River near Prince George
        #------------------------------------------------------------------------

        # tied parameters:
        # (it is important for OSTRICH to find every parameter place holder somewhere in this file)
        # (without this "para_x6" wouldn't be detectable)
        #    para_1_minus_x6 = par_x6 - 1.0 = par_1_minus_x6

        # -Global snow parameters-------------------------------------
        :RainSnowTransition 0 1.0
        :AirSnowCoeff       par_1_minus_x6  # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6
        :AvgAnnualSnow      par_x5          # [mm]  =       CEMANEIGE_X1 =       x5

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
        :EndSoilClasses
        :SoilParameterList
         :Parameters, POROSITY ,  GR4J_X3, GR4J_X2
         :Units     ,     none ,       mm,    mm/d
           [DEFAULT],      1.0 ,   par_x3,  par_x2
        :EndSoilParameterList

        # ----Soil Profiles--------------------------------------------
        #     name, #horizons, (soiltype, thickness) x #horizons
        #     GR4J_X1 is thickness of first layer (SOIL_PROD), here 0.696
        :SoilProfiles
          DEFAULT_P,      4, SOIL_PROD   , par_x1, SOIL_ROUT  ,   0.300, SOIL_TEMP  ,   1.000, SOIL_GW  ,   1.000,
        :EndSoilProfiles

        # ----Vegetation Classes---------------------------------------
        :VegetationClasses
           :Attributes,       MAX_HT,       MAX_LAI,      MAX_LEAF_COND
                :Units,            m,          none,           mm_per_s
               VEG_ALL,           0.0,          0.0,                0.0
        :EndVegetationClasses

        # --Land Use Classes------------------------------------------
        :LandUseClasses
          :Attributes, IMPERM, FOREST_COV
          :Units     ,   frac,       frac
               LU_ALL,    0.0,        0.0
        :EndLandUseClasses
        :LandUseParameterList
         :Parameters, GR4J_X4, MELT_FACTOR
         :Units     ,       d,      mm/d/C
           [DEFAULT],  par_x4,        7.73
        :EndLandUseParameterList
        """
        self.config.rvp.set_tmpl(rvp_tmpl, is_ostrich=True)

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
          {identifier}.rvp.tpl;    {identifier}.rvp
          {identifier}.rvc.tpl;    {identifier}.rvc
          #can be multiple (.rvh, .rvi)
        EndFilePairs

        #Parameter/DV Specification
        BeginParams
          #parameter    init.    low                 high                 tx_in  tx_ost   tx_out
          par_x1        random   {lowerBounds.GR4J_X1}       {upperBounds.GR4J_X1}        none   none     none
          par_x2        random   {lowerBounds.GR4J_X2}       {upperBounds.GR4J_X2}        none   none     none
          par_x3        random   {lowerBounds.GR4J_X3}       {upperBounds.GR4J_X3}        none   none     none
          par_x4        random   {lowerBounds.GR4J_X4}       {upperBounds.GR4J_X4}        none   none     none
          par_x5        random   {lowerBounds.CEMANEIGE_X1}  {upperBounds.CEMANEIGE_X1}   none   none     none
          par_x6        random   {lowerBounds.CEMANEIGE_X2}  {upperBounds.CEMANEIGE_X2}   none   none     none
        EndParams

        BeginTiedParams
          # par_half_x1 = par_x1 * 0.5 * 1000  --> half of it but in [mm] not [m]
          # Xtied = (c1 * X) + c0
          # --> c0 = 0.0
          # --> c1 = 500.
          #
          par_half_x1 1 par_x1 linear 500.0 0.0 free
          #
          # par_1_minus_x6 = - par_x6 + 1.0
          # Xtied = (c1 * X1) + c0
          # --> c0 =  1.0
          # --> c1 = -1.0
          #
          par_1_minus_x6 1 par_x6 linear -1.00 1.00 free
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
        pass
