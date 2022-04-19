from collections import defaultdict
from dataclasses import asdict, field, make_dataclass
from pathlib import Path
from typing import Union, cast

import xarray as xr
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

from ravenpy.config.commands import (
    HRU,
    LU,
    PL,
    SOIL,
    BaseDataCommand,
    BasinIndexCommand,
    Config,
    HRUState,
    LandUseParameterListCommand,
    SoilClassesCommand,
    SoilParameterListCommand,
    Sub,
    VegetationClassesCommand,
    VegetationParameterListCommand,
)
from ravenpy.models.base import Ostrich, Raven

from .gr4jcn import GR4JCN

Sym = Union[Variable, float]


class HBVEC(Raven):
    """Hydrologiska Byråns Vattenbalansavdelning – Environment Canada (HBV-EC)

    References
    ----------
    Lindström, G., et al., 1997. Development and test of the distributed HBV-96 hydrological model.
    Journal of Hydrology, 201 (1–4), 272–288. doi:10.1016/S0022-1694(97)00041-3

    Hamilton, A.S., Hutchinson, D.G., and Moore, R.D., 2000. Estimating winter streamflow using
    conceptual streamflow model. Journal of Cold Regions Engineering, 14 (4), 158–175.
    doi:10.1061/(ASCE)0887-381X(2000)14:4(158)

    Canadian Hydraulics Centre, 2010. Green kenue reference manual.
    Ottawa, Ontario: National Research Council.
    """

    Params = dataclass(
        make_dataclass(
            "Params",
            [
                (f"par_x{i:02}", Sym, field(default=Variable(f"par_x{i:02}")))
                for i in range(1, 22)
            ],
        ),
        config=Config,
    )

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "hbvec")
        super().__init__(*args, **kwds)

        # Initialize symbolic parameters
        P = self.Params()

        self.config.update(
            routing="ROUTE_NONE",
            hrus=(GR4JCN.LandHRU(),),
            soil_classes=[
                SoilClassesCommand.Record(n)
                for n in ["TOPSOIL", "FAST_RES", "SLOW_RES"]
            ],
            vegetation_classes=[
                VegetationClassesCommand.Record("VEG_ALL", 25, 6.0, 5.3)
            ],
            land_use_classes=[LU("LU_ALL", 0, 1)],
            soil_profiles=[
                SOIL(
                    "DEFAULT_P",
                    ["TOPSOIL", "FAST_RES", "SLOW_RES"],
                    [P.par_x17, 100.0, 100.0],
                )
            ],
            soil_parameter_list=[
                SoilParameterListCommand(
                    names=[
                        "POROSITY",
                        "FIELD_CAPACITY",
                        "SAT_WILT",
                        "HBV_BETA",
                        "MAX_CAP_RISE_RATE",
                        "MAX_PERC_RATE",
                        "BASEFLOW_COEFF",
                        "BASEFLOW_N",
                    ],
                    records=[
                        PL(
                            name="[DEFAULT]",
                            vals=[
                                P.par_x05,
                                P.par_x06,
                                P.par_x14,
                                P.par_x07,
                                P.par_x16,
                                0,
                                0,
                                0,
                            ],
                        ),
                        PL(
                            name="FAST_RES",
                            vals=[
                                None,
                                None,
                                0,
                                None,
                                None,
                                P.par_x08,
                                P.par_x09,
                                1 + P.par_x15,
                            ],
                        ),
                        PL(
                            name="SLOW_RES",
                            vals=[None, None, 0, None, None, None, P.par_x10, 1],
                        ),
                    ],
                )
            ],
            vegetation_parameter_list=[
                VegetationParameterListCommand(
                    names=[
                        "MAX_CAPACITY",
                        "MAX_SNOW_CAPACITY",
                        "RAIN_ICEPT_PCT",
                        "SNOW_ICEPT_PCT",
                    ],
                    records=[PL(name="VEG_ALL", vals=[10000, 10000, 0.12, 0.12])],
                )
            ],
            land_use_parameter_list=[
                LandUseParameterListCommand(
                    names=[
                        "MELT_FACTOR",
                        "MIN_MELT_FACTOR",
                        "HBV_MELT_FOR_CORR",
                        "REFREEZE_FACTOR",
                        "HBV_MELT_ASP_CORR",
                        "HBV_MELT_GLACIER_CORR",
                        "HBV_GLACIER_KMIN",
                        "GLAC_STORAGE_COEFF",
                        "HBV_GLACIER_AG",
                    ],
                    records=[
                        PL(
                            name="[DEFAULT]",
                            vals=[
                                P.par_x02,
                                2.2,
                                P.par_x18,
                                P.par_x03,
                                0.48,
                                1.64,
                                0.05,
                                P.par_x19,
                                0.05,
                            ],
                        )
                    ],
                )
            ],
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
        :AdiabaticLapseRate                   {par_x13}
        #                                   HBV_PARA_01, CONSTANT,
        :RainSnowTransition                   {par_x01},      2.0
        #                                   HBV_PARA_04,
        :IrreducibleSnowSaturation            {par_x04}
        #                             HBV_PARA_12=PCALT
        :GlobalParameter PRECIP_LAPSE         {par_x12}

        #---------------------------------------------------------
        # Soil classes

        {soil_classes}

        {soil_parameter_list}

        #---------------------------------------------------------
        # Soil profiles
        # name, layers, (soilClass, thickness) x layers
        #

        {soil_profiles}

        #---------------------------------------------------------
        # Vegetation classes

        {vegetation_classes}

        {vegetation_parameter_list}

        #---------------------------------------------------------
        # LandUse classes

        {land_use_classes}

        {land_use_parameter_list}
        """
        self.config.rvp.set_tmpl(rvp_tmpl)

        #########
        # R V I #
        #########

        rvi_tmpl = """
        :Routing                    {routing}
        :CatchmentRoute      	    ROUTE_TRI_CONVOLUTION

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
        self.config.rvp.set_extra_attributes(par_x11_half=params.par_x11 / 2.0)

        # DH: Why do we need this copy of par_x11 ?
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
        {subbasins}

        {hrus}

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
        self.config.rvp.is_ostrich_tmpl = True

        ####################
        # R V T (OST TMPL) #
        ####################
        self.config.rvt.is_ostrich = True

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
