from dataclasses import asdict, field, make_dataclass
from typing import Any, cast

import xarray as xr
from pydantic.dataclasses import dataclass
from pymbolic import var

from ravenpy.config.commands import (
    HRU,
    LU,
    PL,
    SOIL,
    VEG,
    BaseDataCommand,
    BasinIndexCommand,
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


class HBVECMOD(Raven):
    """
    # para dist       lower     upper     default   informative(0)_or_noninformative(1)_or_ignored(-1)    # comment
    #                 mean      stddev
    X01    uniform    1.0        4.0         1.0        1   # MAX_LAI              | rvp | :VegetationClasses        | MAX_LAI for forest
    X02    uniform    1.0        2.0         1.0        1   # MAX_LAI              | rvp | :VegetationClasses        | MAX_LAI for agriculture, wetland, peatland
    X03    uniform    0.05       0.25        0.21941    1   # Topsoil thickness    | rvp | :SoilProfiles             | thickness of 1st layer = 0.21941 [m]
    X04    uniform    0.05       0.75        0.15725    1   # FAST_RES             | rvp | :SoilProfiles             | thickness of the 2nd layer = 0.15725 [m]
    X05    uniform    0.05       100.0       2.65       1   # SLOW_RES             | rvp | :SoilProfiles             | thickness of the 3rd layer = 2.65 [m]
    X06    uniform   -1.0        1.0         0.0        1   # RAINSNOW_TEMP        | rvp | :GlobalParameter          | rain/snow halfway transition temperature
    X07    uniform    0.0        4.0         1.0        1   # RAINSNOW_DELTA       | rvp | :GlobalParameter          | range of rain-snow transition zone
    X08    uniform    0.0        7.0         4.0        1   # ADIABATIC_LAPSE      | rvp | :GlobalParameter          | adiabatic temperature lapse rate
    X09    uniform    0.04       0.07        0.0464     1   # SNOW_SWI             | rvp | :GlobalParameter          | water saturation fraction of snow
    X10    uniform    0.0        5.0         1.0        1   # PRECIP_LAPSE         | rvp | :GlobalParameter          | PRECIP_LAPSE
    X11    uniform    0.5        2.0         1.0        1   # POROSITY             | rvp | :SoilParameterList        | effective porosity of the soil
    X12    uniform    1.0        7.0         1.0        1   # HBV_BETA             | rvp | :SoilParameterList        | HBV_BETA
    X13    uniform    0.5        2.0         1.0        1   # FIELD_CAPACITY       | rvp | :SoilParameterList        | field capacity saturation of the soil
    X14    uniform    0.01       1.0         0.01       1   # SAT_WILT             | rvp | :SoilParameterList        | SAT_WILT
    X15    uniform    0.01       1.0         0.01       1   # MAX_CAP_RISE_RATE    | rvp | :SoilParameterList        | HBV max capillary rise rate
    X16    uniform    0.01       10.0        1.0        1   # MAX_PERC_RATE        | rvp | :SoilParameterList        | VIC/ARNO/GAWSER percolation rate
    X17    uniform    0.01       1.0         0.03       1   # BASEFLOW_COEFF       | rvp | :SoilParameterList
    | linear baseflow storage/routing coefficient for the FAST_RES
    X18    uniform    0.05       0.1         0.03       1   # BASEFLOW_COEFF       | rvp | :SoilParameterList        | linear baseflow storage/routing coefficient for the SLOW_RES
    X19    uniform    0.5        2.0         1.1        1   # BASEFLOW_N           | rvp | :SoilParameterList        | VIC/ARNO baseflow exponent
    X20    uniform    0.02       0.2         0.02       1   # RAIN_ICEPT_PCT       | rvp | :VegetationParameterList  | relates percentage of throughfall of rain/snow to LAI+SAI
    X21    uniform    0.01       100.0       100.0      1   # MAX_CAPACITY         | rvp | :VegetationParameterList  | maximum canopy storage capacity
    X22    uniform    0.0        5.0         0.01       1   # MELT_FACTOR          | rvp | :LandUseParameterList     | maximum snow melt factor used in degree day models
    X23    uniform    0.01       1.0         0.01       1   # MIN_MELT_FACTOR      | rvp | :LandUseParameterList     | minimum snow melt factor used in degree day models
    X24    uniform    0.001      1.0         0.1        1   # HBV_MELT_FOR_CORR    | rvp | :LandUseParameterList     | HBV snowmelt forest correction (MRF in HBV-EC)
    X25    uniform    1.0        3.0         1.0        1   # REFREEZE_FACTOR      | rvp | :LandUseParameterList     | maximum refreeze factor used in degree day models
    X26    uniform    0.001      1.0         0.1        1   # HBV_MELT_ASP_CORR    | rvp | :LandUseParameterList     | HBV snowmelt aspect correction (AM in HBV-EC)
    X27    uniform    0.001      5.0         0.01       1   # GLAC_STORAGE_COEF    | rvp | :LandUseParameterList     | maximum linear storage coefficient for glacial melt
    """

    # TODO: enforce upper case for param names - needed for ostrich tied params.
    Params = dataclass(
        make_dataclass(
            "Params",
            [(f"X{i:02}", Any, field(default=var(f"X{i:02}"))) for i in range(1, 28)],
        )
    )

    @dataclass
    class LandHRU(HRU):
        land_use_class: str = "CONIFEROUS_FOREST"
        veg_class: str = "CONIFEROUS_FOREST"
        soil_profile: str = "sand"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        hru_type: str = "land"

    @dataclass
    class LakeHRU(HRU):
        land_use_class: str = "WATER"
        veg_class: str = "WATER"
        soil_profile: str = "LAKE"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "[NONE]"
        hru_type: str = "lake"

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "hbvecmod")
        super().__init__(*args, **kwds)

        P = HBVECMOD.Params

        self.config.update(hrus=(HBVECMOD.LandHRU(),))

        #########
        # R V P #
        #########

        rvp_tmpl = """
        # These parameters appear in tied parameters only, so we print them here for debugging.
        # x01: {X01}
        # x02: {X02}
        # x11: {X11}
        # x13: {X13}
        # x14: {X14}
        # x23: {X23}

        #------------------------------------------------------------------------
        # Global parameters
        #
        # --- Parameter Specification ----------------------------
        #                                  X6=[-1.0,1.0];  default = 0.0
        :GlobalParameter RAINSNOW_TEMP     {X06}
        #                                  X7 =[0.0,4.0];  default = 1.0
        :GlobalParameter RAINSNOW_DELTA    {X07}
        #                                  X8=[0.0,7.0];   default = 4.0
        :GlobalParameter ADIABATIC_LAPSE   {X08}
        #                                  X9=[0.04,0.07]; default = 0.04
        :GlobalParameter SNOW_SWI          {X09}

        # The PRECIP_LAPSE is the parameter for the simple linear precipitation lapse rate [mm/d/km], as used in the
        # OROCORR_SIMPLELAPSE orographic correction algorithm
        #                                  X10=[0.0,5.0] ; default = 1.002140; X10 NOT calibrated in Famine
        :GlobalParameter PRECIP_LAPSE      {X10}
        #
        {avg_annual_runoff}
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

        # List of channel profiles
        {channel_profiles}
        """
        p = self.config.rvp
        p.set_tmpl(rvp_tmpl)
        p.soil_classes = [
            SoilClassesCommand.Record(n) for n in ["TOPSOIL", "FAST_RES", "SLOW_RES"]
        ]

        p.vegetation_classes = [
            VEG("AGRICULTURE", 5, 1 + P.X02, 5),
            VEG("BARE_SOIL", 0, 1, 5),
            VEG("CONIFEROUS_FOREST", 5, 1 + P.X01, 5),
            VEG("DECIDUOUS_FOREST", 5, 1 + P.X01, 5),
            VEG("IMPERMEABLE_SURFACE", 0, 1 + P.X01, 5),
            VEG("PEATLAND", 2, 1 + P.X02, 5),
            VEG("WATER", 0, 1, 5),
            VEG("WETLAND", 4, 1 + P.X02, 5),
        ]

        p.land_use_classes = [
            LU("AGRICULTURE", 0, 0.5),
            LU("BARE_SOIL", 0, 0),
            LU("CONIFEROUS_FOREST", 0, 1),
            LU("DECIDUOUS_FOREST", 0, 1),
            LU("IMPERMEABLE_SURFACE", 1, 0),
            LU("PEATLAND", 0, 0),
            LU("WATER", 0, 0),
            LU("WETLAND", 0, 0),
        ]

        p.soil_profiles = [
            SOIL(
                "sand",
                ["TOPSOIL", "FAST_RES", "SLOW_RES"],
                [P.X03, P.X04, P.X05],
            ),
            SOIL(
                "loamy_sand",
                ["TOPSOIL", "FAST_RES", "SLOW_RES"],
                [P.X03, P.X04, P.X05],
            ),
            SOIL(
                "loam",
                ["TOPSOIL", "FAST_RES", "SLOW_RES"],
                [P.X03, P.X04, P.X05],
            ),
            SOIL(
                "sandy_loam",
                ["TOPSOIL", "FAST_RES", "SLOW_RES"],
                [P.X03, P.X04, P.X05],
            ),
            SOIL("LAKE", [], []),
        ]

        p.soil_parameter_list = [
            SoilParameterListCommand(
                names=[
                    "POROSITY",
                    "HBV_BETA",
                    "FIELD_CAPACITY",
                    "SAT_WILT",
                    "MAX_CAP_RISE_RATE",
                    "MAX_PERC_RATE",
                    "BASEFLOW_COEFF",
                    "BASEFLOW_N",
                ],
                records=[
                    PL(
                        name="TOPSOIL",
                        vals=[
                            P.X11 * 0.2,
                            P.X12,
                            P.X13 * 0.5,
                            P.X14 * 0.01,
                            0,
                            0,
                            0,
                            0,
                        ],
                    ),
                    PL(
                        name="FAST_RES",
                        vals=[
                            P.X11 * 0.4,
                            P.X12,
                            P.X13 * 0.4,
                            P.X14 * 0.01,
                            P.X15,
                            P.X16,
                            P.X17,
                            P.X19,
                        ],
                    ),
                    PL(
                        name="SLOW_RES",
                        vals=[
                            P.X11 * 0.3,
                            P.X12,
                            P.X13 * 0.2,
                            P.X14 * 0.01,
                            P.X15,
                            P.X16,
                            P.X18,
                            P.X19,
                        ],
                    ),
                ],
            )
        ]

        p.vegetation_parameter_list = [
            VegetationParameterListCommand(
                names=[
                    "SAI_HT_RATIO",
                    "RAIN_ICEPT_PCT",
                    "SNOW_ICEPT_PCT",
                    "MAX_CAPACITY",
                    "MAX_SNOW_CAPACITY",
                ],
                records=[
                    PL(name="AGRICULTURE", vals=[2, 0, 0, 0, 0]),
                    PL("BARE_SOIL", vals=[1, 0, 0, 0, 0]),
                    PL("CONIFEROUS_FOREST", vals=[3, P.X20, P.X20, P.X21, P.X21]),
                    PL("DECIDUOUS_FOREST", vals=[3, P.X20, P.X20, P.X21, P.X21]),
                    PL("IMPERMEABLE_SURFACE", vals=[1, 0, 0, 0, 0]),
                    PL("PEATLAND", vals=[1, 0, 0, 0, 0]),
                    PL("WATER", vals=[1, 0, 0, 0, 0]),
                    PL("WETLAND", vals=[2, 0, 0, 0, 0]),
                ],
            )
        ]

        p.land_use_parameter_list = [
            LandUseParameterListCommand(
                names=[
                    "MIN_MELT_FACTOR",
                    "MELT_FACTOR",
                    "HBV_MELT_FOR_CORR",
                    "REFREEZE_FACTOR",
                    "HBV_MELT_ASP_CORR",
                ],
                records=[
                    PL(
                        "AGRICULTURE",
                        vals=[0, P.X22, P.X24, P.X25, 0.5 * P.X26],
                    ),
                    PL(
                        "BARE_SOIL",
                        vals=[0, P.X22, 0.5, P.X25, 0.4 * P.X26],
                    ),
                    PL(
                        "CONIFEROUS_FOREST",
                        vals=[
                            3 * P.X23,
                            P.X22 + 3 * P.X23,
                            5 * P.X24,
                            2 * P.X25,
                            P.X26,
                        ],
                    ),
                    PL(
                        "DECIDUOUS_FOREST",
                        vals=[
                            3 * P.X23,
                            P.X22 + 3 * P.X23,
                            5 * P.X24,
                            2 * P.X25,
                            P.X26,
                        ],
                    ),
                    PL(
                        "IMPERMEABLE_SURFACE",
                        vals=[0, P.X22, 0.2, P.X25, 0.4 * P.X26],
                    ),
                    PL("PEATLAND", vals=[0, P.X22, 0.1, P.X25, 0.4 * P.X26]),
                    PL("WATER", vals=[0, P.X22, 0.2, P.X25, 0.4 * P.X26]),
                    PL(
                        "WETLAND",
                        vals=[3 * P.X23, P.X22 + 3 * P.X23, 0.4, P.X25, P.X26],
                    ),
                ],
            ),
            LandUseParameterListCommand(
                names=[
                    "HBV_MELT_GLACIER_CORR",
                    "HBV_GLACIER_KMIN",
                    "HBV_GLACIER_AG",
                    "GLAC_STORAGE_COEFF",
                ],
                records=[PL("[DEFAULT]", vals=[0, 0, 0, P.X27])],
            ),
        ]

        #########
        # R V I #
        #########

        rvi_tmpl = """
        :Interpolation              INTERP_NEAREST_NEIGHBOR

        :Routing                    {routing}
        :CatchmentRoute      	    {catchment_route}

        :Evaporation         	    {evaporation}  # PET_FROM_MONTHLY
        :OW_Evaporation      	    {ow_evaporation}  # PET_FROM_MONTHLY
        :SWRadiationMethod   	    SW_RAD_DEFAULT
        :SWCloudCorrect      	    SW_CLOUD_CORR_NONE
        :SWCanopyCorrect     	    SW_CANOPY_CORR_NONE
        :LWRadiationMethod   	    LW_RAD_DEFAULT
        :RainSnowFraction           {rain_snow_fraction}
        :PotentialMeltMethod 	    POTMELT_HBV
        :OroTempCorrect      	    OROCORR_NONE
        :OroPrecipCorrect    	    OROCORR_NONE
        :OroPETCorrect       	    OROCORR_NONE
        :CloudCoverMethod    	    CLOUDCOV_NONE
        :PrecipIceptFract    	    PRECIP_ICEPT_USER

        :SoilModel                  SOIL_MULTILAYER 3

        #------------------------------------------------------------------------
        # Soil Layer Alias Definitions
        #
        :Alias       FAST_RESERVOIR SOIL[1]
        :Alias       SLOW_RESERVOIR SOIL[2]
        :LakeStorage SLOW_RESERVOIR

        #------------------------------------------------------------------------
        # Hydrologic process order for HBV-EC Mod Emulation
        #
        :HydrologicProcesses
          :SnowRefreeze      FREEZE_DEGREE_DAY  SNOW_LIQ        SNOW
          :Precipitation     PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE
          :CanopyEvaporation CANEVP_ALL         CANOPY          ATMOSPHERE
          :CanopySnowEvap    CANEVP_ALL         CANOPY_SNOW     ATMOSPHERE
          :SnowBalance       SNOBAL_SIMPLE_MELT SNOW            SNOW_LIQ
            :-->Overflow     RAVEN_DEFAULT      SNOW_LIQ        PONDED_WATER
          :Flush             RAVEN_DEFAULT      SURFACE_WATER   FAST_RESERVOIR
            :-->Conditional HRU_TYPE IS_NOT GLACIER
          :Infiltration      INF_HBV            PONDED_WATER    MULTIPLE
          :SoilEvaporation   SOILEVAP_HBV       SOIL[0]         ATMOSPHERE
          :CapillaryRise     RISE_HBV           FAST_RESERVOIR 	SOIL[0]
          :LakeEvaporation   LAKE_EVAP_BASIC    SLOW_RESERVOIR  ATMOSPHERE
          :Percolation       PERC_CONSTANT      FAST_RESERVOIR 	SLOW_RESERVOIR
          :Baseflow          BASE_POWER_LAW     FAST_RESERVOIR  SURFACE_WATER
          :Baseflow          BASE_LINEAR        SLOW_RESERVOIR  SURFACE_WATER
        :EndHydrologicProcesses
        """
        self.config.rvi.set_tmpl(rvi_tmpl)

        self.config.rvi.routing = "ROUTE_DIFFUSIVE_WAVE"
        self.config.rvi.catchment_route = "ROUTE_DUMP"
        self.config.rvi.evaporation = "PET_HARGREAVES_1985"
        self.config.rvi.ow_evaporation = "PET_HARGREAVES_1985"
        self.config.rvi.rain_snow_fraction = "RAINSNOW_HBV"

        #########
        # R V H #
        #########
        self.config.rvh.subbasins = (
            Sub(
                subbasin_id=1,
                name="sub_001",
                downstream_id=-1,
                profile="None",
                gauged=True,
            ),
        )

    def derived_parameters(self):
        self._monthly_average()

        if not self.config.rvc.hru_states:
            # If self.rvc.hru_states is set, it means that we are using `resume()` and we don't
            # want to interfere
            for hru in self.config.rvh.hrus:
                if isinstance(hru, HBVECMOD.LandHRU) or hru.soil_profile in [
                    "SAND",
                    "LOAMY_SAND",
                    "LOAM",
                    "SANDY_LOAM",
                ]:
                    self.config.rvc.hru_states[hru.hru_id] = HRUState(index=hru.hru_id)
                elif isinstance(hru, HBVECMOD.LakeHRU) or hru.soil_profile in ["LAKE"]:
                    self.config.rvc.hru_states[hru.hru_id] = HRUState(index=hru.hru_id)

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


class HBVECMOD_OST(Ostrich, HBVECMOD):
    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "hbvecmod-ost")
        super().__init__(*args, **kwds)
        P = HBVECMOD.Params

        self.config.update(
            algorithm="DDS",
            max_iterations=50,
            run_name="run",
            run_index=0,
            suppress_output=True,
        )
        p = self.config.rvp
        p.soil_parameter_list = [
            SoilParameterListCommand(
                names=[
                    "POROSITY",
                    "HBV_BETA",
                    "FIELD_CAPACITY",
                    "SAT_WILT",
                    "MAX_CAP_RISE_RATE",
                    "MAX_PERC_RATE",
                    "BASEFLOW_COEFF",
                    "BASEFLOW_N",
                ],
                records=[
                    PL(
                        name="TOPSOIL",
                        vals=[
                            P.X11 * 0.2,
                            P.X12,
                            P.X13 * 0.5,
                            P.X14 * 0.01,
                            0,
                            0,
                            0,
                            0,
                        ],
                    ),
                    PL(
                        name="FAST_RES",
                        vals=[
                            P.X11 * 0.4,
                            P.X12,
                            P.X13 * 0.4,
                            P.X14 * 0.01,
                            P.X15,
                            P.X16,
                            P.X17,
                            P.X19,
                        ],
                    ),
                    PL(
                        name="SLOW_RES",
                        vals=[
                            P.X11 * 0.3,
                            P.X12,
                            P.X13 * 0.2,
                            P.X14 * 0.01,
                            P.X15,
                            P.X16,
                            P.X18,
                            P.X19,
                        ],
                    ),
                ],
            )
        ]

        #########
        # R V C #
        #########

        ####################
        # R V H (OST TMPL) #
        ####################
        self.config.rvh.is_ostrich = True

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
          #can be multiple (.rvh, .rvi)
        EndFilePairs

        #Parameter/DV Specification
        BeginParams
          #parameter      init.    low      high    tx_in  tx_ost  tx_out   # in HBV called
          X01         random   {lowerBounds.X01}  {upperBounds.X01}  none   none     none
          X02         random   {lowerBounds.X02}  {upperBounds.X02}  none   none     none
          X03         random   {lowerBounds.X03}  {upperBounds.X03}  none   none     none
          X04         random   {lowerBounds.X04}  {upperBounds.X04}  none   none     none
          X05         random   {lowerBounds.X05}  {upperBounds.X05}  none   none     none
          X06         random   {lowerBounds.X06}  {upperBounds.X06}  none   none     none
          X07         random   {lowerBounds.X07}  {upperBounds.X07}  none   none     none
          X08         random   {lowerBounds.X08}  {upperBounds.X08}  none   none     none
          X09         random   {lowerBounds.X09}  {upperBounds.X09}  none   none     none
          X10         random   {lowerBounds.X10}  {upperBounds.X10}  none   none     none
          X11         random   {lowerBounds.X11}  {upperBounds.X11}  none   none     none
          X12         random   {lowerBounds.X12}  {upperBounds.X12}  none   none     none
          X13         random   {lowerBounds.X13}  {upperBounds.X13}  none   none     none
          X14         random   {lowerBounds.X14}  {upperBounds.X14}  none   none     none
          X15         random   {lowerBounds.X15}  {upperBounds.X15}  none   none     none
          X16         random   {lowerBounds.X16}  {upperBounds.X16}  none   none     none
          X17         random   {lowerBounds.X17}  {upperBounds.X17}  none   none     none
          X18         random   {lowerBounds.X18}  {upperBounds.X18}  none   none     none
          X19         random   {lowerBounds.X19}  {upperBounds.X19}  none   none     none
          X20         random   {lowerBounds.X20}  {upperBounds.X20}  none   none     none
          X21         random   {lowerBounds.X21}  {upperBounds.X21}  none   none     none
          X22         random   {lowerBounds.X22}  {upperBounds.X22}  none   none     none
          X23         random   {lowerBounds.X23}  {upperBounds.X23}  none   none     none
          X24         random   {lowerBounds.X24}  {upperBounds.X24}  none   none     none
          X25         random   {lowerBounds.X25}  {upperBounds.X25}  none   none     none
          X26         random   {lowerBounds.X26}  {upperBounds.X26}  none   none     none
          X27         random   {lowerBounds.X27}  {upperBounds.X27}  none   none        none
        EndParams

        {tied_params}

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
