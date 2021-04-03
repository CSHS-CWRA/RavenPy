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
        # self.rvc = RVC(soil0=None, soil1=15)
        self.config.rvc.soil0 = None
        self.config.rvc.soil1 = 15

        # self.rvd = RV(one_minus_CEMANEIGE_X2=None, GR4J_X1_hlf=None)

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

    def run(self, ts, overwrite=False, **kwds):
        """
        This is a hook into `Raven.run` for this particular subclass, which
        allows the support of legacy HRU-related keywords in the model.__call__
        interface.
        """
        hru_attrs = {}
        for k in ["area", "latitude", "longitude", "elevation"]:
            v = kwds.pop(k, None)
            if v:
                # It seems that `v` is a list when running via a WPS interface
                hru_attrs[k] = v[0] if isinstance(v, list) else v
        if hru_attrs:
            self.config.rvh.hrus = (GR4JCN.LandHRU(**hru_attrs),)

        return super().run(ts, overwrite=overwrite, **kwds)


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
    templates = tuple((Path(__file__).parent / "raven-mohyse").glob("*.rv?"))

    params = namedtuple(
        "MOHYSEParams", ", ".join(["par_x{:02}".format(i) for i in range(1, 11)])
    )

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvp = RV(params=MOHYSE.params(*((None,) * 10)))
        self.rvh = RV(
            name=None, area=None, elevation=None, latitude=None, longitude=None
        )
        self.rvt = RVT()
        self.rvi = RVI(evaporation="PET_MOHYSE", rain_snow_fraction="RAINSNOW_DATA")
        self.rvc = RVC(
            hru_state=HRUState(),
            basin_state=BasinIndexCommand(),
        )
        self.rvd = RV(par_rezi_x10=None)

    def derived_parameters(self):
        self.rvd["par_rezi_x10"] = 1.0 / self.rvp.params.par_x10


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
    templates = tuple((Path(__file__).parent / "raven-hmets").glob("*.rv?"))

    params = namedtuple(
        "HMETSParams",
        (
            "GAMMA_SHAPE",
            "GAMMA_SCALE",
            "GAMMA_SHAPE2",
            "GAMMA_SCALE2",
            "MIN_MELT_FACTOR",
            "MAX_MELT_FACTOR",
            "DD_MELT_TEMP",
            "DD_AGGRADATION",
            "SNOW_SWI_MIN",
            "SNOW_SWI_MAX",
            "SWI_REDUCT_COEFF",
            "DD_REFREEZE_TEMP",
            "REFREEZE_FACTOR",
            "REFREEZE_EXP",
            "PET_CORRECTION",
            "HMETS_RUNOFF_COEFF",
            "PERC_COEFF",
            "BASEFLOW_COEFF_1",
            "BASEFLOW_COEFF_2",
            "TOPSOIL",
            "PHREATIC",
        ),
    )

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvp = RV(params=HMETS.params(*((None,) * len(HMETS.params._fields))))
        self.rvt = RVT()
        self.rvh = RV(
            name=None, area=None, elevation=None, latitude=None, longitude=None
        )
        self.rvi = RVI(evaporation="PET_OUDIN", rain_snow_fraction="RAINSNOW_DATA")
        self.rvc = RVC(soil0=None, soil1=None, basin_state=BasinIndexCommand())
        self.rvd = RV(
            TOPSOIL_m=None,
            PHREATIC_m=None,
            SUM_MELT_FACTOR=None,
            SUM_SNOW_SWI=None,
            TOPSOIL_hlf=None,
            PHREATIC_hlf=None,
        )

    def derived_parameters(self):
        self.rvd["TOPSOIL_hlf"] = self.rvp.params.TOPSOIL * 0.5
        self.rvd["PHREATIC_hlf"] = self.rvp.params.PHREATIC * 0.5
        self.rvd["TOPSOIL_m"] = self.rvp.params.TOPSOIL / 1000.0
        self.rvd["PHREATIC_m"] = self.rvp.params.PHREATIC / 1000.0
        self.rvd[
            "SUM_MELT_FACTOR"
        ] = self.rvp.params.MAX_MELT_FACTOR  # self.rvp.params.MIN_MELT_FACTOR +
        self.rvd[
            "SUM_SNOW_SWI"
        ] = self.rvp.params.SNOW_SWI_MAX  # self.rvp.params.SNOW_SWI_MIN +

        # Default initial conditions if none are given
        if self.rvc.hru_state is None:
            soil0 = (
                self.rvd["TOPSOIL_hlf"] if self.rvc.soil0 is None else self.rvc.soil0
            )
            soil1 = (
                self.rvd["PHREATIC_hlf"] if self.rvc.soil1 is None else self.rvc.soil1
            )
            self.rvc.hru_state = HRUState(soil0=soil0, soil1=soil1)


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
    templates = tuple((Path(__file__).parent / "raven-hbv-ec").glob("*.rv?"))

    params = namedtuple("HBVECParams", ("par_x{:02}".format(i) for i in range(1, 22)))

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvp = RV(params=HBVEC.params(*((None,) * len(HBVEC.params._fields))))
        self.rvd = RV(
            one_plus_par_x15=None,
            par_x11_half=None,
        )
        self.rvt = RVT()
        self.rvh = RV(
            name=None, area=None, elevation=None, latitude=None, longitude=None
        )
        self.rvi = RVI(
            evaporation="PET_FROMMONTHLY",
            ow_evaporation="PET_FROMMONTHLY",
            rain_snow_fraction="RAINSNOW_HBV",
        )
        self.rvc = RVC(soil2=0.50657, qout=1)

    def derived_parameters(self):
        self.rvd["one_plus_par_x15"] = self.rvp.params.par_x15 + 1.0
        self.rvd["par_x11_half"] = self.rvp.params.par_x11 / 2.0

        self.rvt["raincorrection"] = self.rvp.params.par_x20
        self.rvt["snowcorrection"] = self.rvp.params.par_x21

        self._monthly_average()

        # Default initial conditions if none are given
        if self.rvc.hru_state is None:
            self.rvc.hru_state = HRUState(soil2=self.rvc.soil2)
        if self.rvc.basin_state is None:
            self.rvc.basin_state = BasinIndexCommand(qout=(self.rvc.qout,))

    # TODO: Support index specification and unit changes.
    def _monthly_average(self):

        if (
            self.rvi.evaporation == "PET_FROMMONTHLY"
            or self.rvi.ow_evaporation == "PET_FROMMONTHLY"
        ):
            # If this fails, it's likely the input data is missing some necessary variables (e.g. evap).
            tas_cmd = self.rvt.var_cmds.get("tas")
            tasmin_cmd = self.rvt.var_cmds.get("tasmin")
            tasmax_cmd = self.rvt.var_cmds.get("tasmax")
            evspsbl_cmd = self.rvt.var_cmds.get("evspsbl")

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

            self.rvt["monthly_ave_evaporation"] = tuple(mae)
            self.rvt["monthly_ave_temperature"] = tuple(mat)


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
    templates = tuple((Path(__file__).parent / "raven-blended").glob("*.rv?"))

    params = namedtuple(
        "BLENDEDParams",
        ", ".join(
            ["par_x{:02}".format(i) for i in range(1, 36)]
            + ["par_r{:02}".format(i) for i in range(1, 9)]
        ),
    )

    @dataclass
    class HRU(HRU):
        land_use_class: str = "FOREST"
        veg_class: str = "FOREST"
        soil_profile: str = "DEFAULT_P"
        aquifer_profile: str = "[NONE]"
        terrain_class: str = "DEFAULT_T"

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.rvp = RVP(
            params=BLENDED.params(*((None,) * len(BLENDED.params._fields))),
            land_use_classes=(
                LU("FOREST", impermeable_frac=0.0, forest_coverage=0.02345),
            ),
        )
        self.rvh = RVH(hrus=(BLENDED.HRU(),))
        self.rvt = RVT()
        self.rvi = RVI(evaporation="PET_OUDIN", rain_snow_fraction="RAINSNOW_HBV")
        self.rvc = RVC(soil0=None, soil1=None, basin_state=BasinIndexCommand())
        self.rvd = RV(
            TOPSOIL_mm=None,
            PHREATIC_mm=None,
            TOPSOIL_hlf=None,
            PHREATIC_hlf=None,
            POW_X04=None,
            POW_X11=None,
            SUM_X09_X10=None,
            SUM_X13_X14=None,
            SUM_X24_X25=None,
        )

    def derived_parameters(self):
        self.rvd["TOPSOIL_hlf"] = self.rvp.params.par_x29 * 0.5 * 1000.0
        self.rvd["PHREATIC_hlf"] = self.rvp.params.par_x30 * 0.5 * 1000.0
        self.rvd["TOPSOIL_mm"] = self.rvp.params.par_x29 * 1000.0
        self.rvd["PHREATIC_mm"] = self.rvp.params.par_x30 * 1000.0
        self.rvd["SUM_X09_X10"] = self.rvp.params.par_x10  # + self.rvp.params.par_x09
        self.rvd["SUM_X13_X14"] = self.rvp.params.par_x14  # + self.rvp.params.par_x13
        self.rvd["SUM_X24_X25"] = self.rvp.params.par_x25  # + self.rvp.params.par_x24
        # 10.0**self.rvp.params.par_x04  #
        self.rvd["POW_X04"] = self.rvp.params.par_x04
        # 10.0**self.rvp.params.par_x11  #
        self.rvd["POW_X11"] = self.rvp.params.par_x11

        # Default initial conditions if none are given
        if self.rvc.hru_state is None:
            soil0 = (
                self.rvd["TOPSOIL_hlf"] if self.rvc.soil0 is None else self.rvc.soil0
            )
            soil1 = (
                self.rvd["PHREATIC_hlf"] if self.rvc.soil1 is None else self.rvc.soil1
            )
            self.rvc.hru_state = HRUState(soil0=soil0, soil1=soil1)

        self.rvt.raincorrection = self.rvp.params.par_x33
        self.rvt.snowcorrection = self.rvp.params.par_x34


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
